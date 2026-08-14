using Plots
using Printf

include("level.jl")

function advance(grid::Grid)

    # Compute the boundary conditions 
    compute_boundary_condition_on_grid!(grid)

    # Definition of the locql variables 
    new_u::Vector{Float64} = Vector{Float64}(undef, grid.n_points + 2 * grid.n_gzones)

    for ind in grid.box.inside_ind
        new_u[ind] = grid.u[ind] - grid.dt * (compute_flux(grid.u, ind)) / grid.dx
    end

    return new_u
    
end

function lohner_err_estimator!(grid::Grid, epsilon::Float64)

    for ind in grid.box.inside_ind

        if ind == grid.box.inside_ind[1]
            grid.local_error[ind] = 0.0
        elseif ind == grid.box.inside_ind[end]
            grid.local_error[ind] = 0.0
        else
            
            num::Float64 = abs(grid.u[ind + 1] - 2.0 * grid.u[ind] + grid.u[ind - 1])
            denom::Float64 = abs(grid.u[ind + 1] - grid.u[ind]) + abs(grid.u[ind] - grid.u[ind - 1]) + epsilon * (abs(grid.u[ind + 1]) + 2.0 * abs(grid.u[ind]) + abs(grid.u[ind - 1]))

            grid.local_error[ind] = num / denom 

        end
    end

    return nothing
end

function initialize!(level::Union{Level, Nothing})

    if !isnothing(level)
        initialize_level!(level)

        if level.depth < max_depth
            
            compute_local_error!(level)

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
    title = @sprintf("iteration = %d", iteration),
    marker = :xcross,
    label = nothing,
    seriestype = :line,
    linecolor = :orange,
    markercolor = :orange,
    xlims = [0.0, 1.0],
    ylims = [0.0, 1.0]
    )

    plot!(
        twinx(fig),
        position_list, depth_list,
        marker = :xcross,
        label = nothing,
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

dir_out::String = "../image/criteria2/test3/"

out::Bool = true

n_points::Int64 = 128
n_gzones::Int64 = 1
x0::Float64 = 0.0
l::Float64 = 1.0
dx::Float64 = (l - x0) / n_points
dt::Float64 =   dx #/ (2.0^(max_depth - 1))

tfinal::Float64 = 0.5

base_level::Level = Level(n_points, n_gzones, x0, dt, dx)



# Initialization 
initialize!(base_level)

if out 
    vizulisation(base_level, 0, dir_out)
end

base_grid::Grid = base_level.child.grid_list[1]

for ind in base_grid.box.inside_ind
    if base_grid.local_error[ind] > 0.0
        println(base_grid.x[ind])
    end
end







