using Plots
using Printf

include("level.jl")

function advance(grid::Grid)

    # Definition of the local variables 
    new_u::Vector{Float64} = Vector{Float64}(undef, grid.n_points + 2 * grid.n_gzones)

    # Update the inside points
    for ind in grid.box.inside_ind
        new_u[ind] = grid.u[ind] - grid.dt * (compute_flux(grid.u, ind)) / grid.dx
    end
    return new_u
    
end


function initialize!(level::Union{Level, Nothing})

    if !isnothing(level)
        initialize_level!(level)

        if level.depth < max_depth

            if level.depth > 0
                set_initial_local_error!(level)
            end

            regrid_level!(level)

            initialize!(level.child)
        end
    end

    return nothing
end


function test(level::Union{Level, Nothing})

    if !isnothing(level)
        println("depth: ", level.depth, ", number of grid: ", length(level.grid_list), ", x0 = ", level.grid_list[1].x0)
        test(level.child)
    end

    return true
end


function vizulisation(base_level::Level, iteration::Int64, dir_out::String)
    
    # Declaration of the local variables
    position_list::Vector{Float64} = []
    solution_list::Vector{Float64} = []
    depth_list::Vector{Int64} = []

    # Update of the lists
    output!(base_level.grid_list[1], 1, position_list, solution_list, depth_list)

    # Plot 
    fig = plot(
        position_list, solution_list,
        title = @sprintf("time = %f ", base_level.grid_list[1].t),
        marker = :xcross,
        label = "solution",
        legend = :topright,
        seriestype = :line,
        linecolor = :orange,
        markercolor = :orange,
        xlims = [0.0, 1.0], 
        ylims = [0.0, 1.0]
    )

    plot!(
        fig,
        [NaN], [NaN],
        marker = :xcross,
        label = "level",
        seriestype = :line,
        linecolor = :blue,
        markercolor = :blue
    )

    plot!(
        twinx(fig),
        position_list, depth_list,
        marker = :xcross,
        label = nothing,
        legend = :topright,
        seriestype = :line,
        linecolor = :blue,
        markercolor = :blue,
        xlims = [0.0, 1.0],
        ylims = [1, max_depth],
        ylabel = "depth"
    )

    file_out::String = dir_out * @sprintf("sol_%d.png", iteration) 

    savefig(fig, file_out)

    position_list = []
    solution_list =[]
    depth_list = []

    return nothing
end






# Declaration of the variables 

dir_out::String = "../image/criteria3/test1/"

out::Bool = false

n_points::Int64 = 64
n_gzones::Int64 = 1
x0::Float64 = 0.0
l::Float64 = 1.0
dx::Float64 = (l - x0) / n_points
dt::Float64 =   cfl * dx / (2.0^(max_depth))

tfinal::Float64 = 0.5


# Initialization
level_0::Level = Level(n_points, n_gzones, x0, dt, dx)

initialize!(level_0)

if out
    vizulisation(level_0.child, 0, dir_out)
end


# Compute the solution 
iteration::Int64 = 0

while level_0.grid_list[1].t < tfinal

    # Update the iteration
    global iteration 
    iteration += 1

    # Update the time step value
    update_dt!(level_0)

    # Advance the solution 
    evolve_level!(level_0)



    if (out ) && (iteration % 10 == 0)
        vizulisation(level_0.child, iteration, dir_out)
    end
end


# Compute the global error 

position_list::Vector{Float64} = []
solution_list::Vector{Float64} = []
depth_list::Vector{Int64} = []

output!(level_0.child.grid_list[1], 1, position_list, solution_list, depth_list)

err_glob::Float64 = 0.0

for ind in eachindex(position_list)
    global err_glob
    err_glob += (dx / 2.0^(depth_list[ind])) * (solution_list[ind] - exact_solution(level_0.grid_list[1].t, position_list[ind]))^2
end



println("max_depth: ", max_depth, ", global error: ", sqrt(err_glob))










