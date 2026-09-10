using Plots
using Printf

using DataFrames, GLM

include("level.jl")

function advance!(grid::Grid)

    # Update the inside points
    for ind in grid.box.inside_ind
        grid.u[ind] = grid.u_p[ind] - grid.dt * (compute_flux(grid.u_p, ind + 1) - compute_flux(grid.u_p, ind)) / grid.dx
    end
    
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

            regrid_level!(level)

            initialize!(level.child)
        end
    end

    return nothing
end


function test(level::Union{Level, Nothing})

    if !isnothing(level)
        println("depth: ", level.depth, ", number of grid: ", length(level.grid_list), ", x0 = ", level.grid_list[1].x0)
        println(level.grid_list[1].u)
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

function compute_order(n_points::Int64, n_iteration::Int64, n_gzones::Int64, x0::Float64, l::Float64, tfinal::Float64)

    # Definition of the local variables
    base_level::Level = Level()
    err_glob::Float64 = 0.0
    err_glb_list::Vector{Float64} = []
    dx_list::Vector{Float64} = []
    dt::Float64 = 0.0

    
    # Loop on the iteration

    for ite in 1:n_iteration
        new_n_points::Int64 = n_points * (2 ^(ite - 1))
        dx::Float64 = (l - x0 ) / new_n_points
        base_level = Level()
        base_level = Level(new_n_points, n_gzones, x0, dt, dx)

        # Initialization
        initialize!(base_level)
        
        # Compute the solution
        iteration = 0

        while base_level.grid_list[1].t < tfinal 

            # Update the iteration count 
            iteration += 1

            # Update the time step value 
            update_dt!(base_level)    
            
            # Evolve the solution 
            evolve_level!(base_level)
        end

        # Compute global error 
        position_list::Vector{Float64} = []
        solution_list::Vector{Float64} = []
        depth_list::Vector{Int64} = []

        output!(base_level.grid_list[1], 1, position_list, solution_list, depth_list)

        err_glob = 0.0

        for ind in eachindex(position_list)
            err_glob += (dx / 2.0^(depth_list[ind] - 1)) * (solution_list[ind] - exact_solution(base_level.grid_list[1].t, position_list[ind]))^2
        end

        err_glob = sqrt(err_glob)

        # Save 
        push!(err_glb_list, err_glob)
        push!(dx_list, dx)

    end

    return dx_list, err_glb_list
end





# Declaration of the variables 

dir_out::String = "../image/criteria4/test3/"

out::Bool = false

n_points::Int64 = 2048
n_gzones::Int64 = 1
x0::Float64 = 0.0
l::Float64 = 1.0
dx::Float64 = (l - x0) / n_points
dt::Float64 =   0.0

tfinal::Float64 = 0.5



color = [:blue, :orange, :green, :purple]



dx_list, err_list = compute_order(n_points, 5, n_gzones, x0, l, tfinal)


file_out::String = "test.png"

# Compute the order

data =  DataFrame(dx = dx_list, err = err_list)

data.logdx = log.(data.dx)
data.logerr = log.(data.err)

model = lm(@formula(logerr ~ logdx), data)

a = coef(model)[2]
b = coef(model)[1]

data.fitted_logerr = predict(model)
data.fitted_err = exp.(data.fitted_logerr)

# Plot 

fig = plot(
    data.dx, data.err,
    seriestype = :scatter,
    xscale = :log10, yscale = :log10,
    xlims = [10^(-5), 10^(-2)],
    xlabel = "dx", ylabel = "global error", 
    marker = (:xcross, 6),
    markercolor = color[max_depth],
    label = @sprintf("%d level", max_depth),
    legend = :topleft
)

plot!(
    fig, 
    data.dx, data.fitted_err,
    xscale = :log10, yscale = :log10,
    linestyle = :dash,
    linecolor = color[max_depth],
    label = @sprintf("slope %.3f, r2 %.3f", a, r2(model)),
    legend = :topleft
)

while max_depth < 4
    global max_depth += 1
    global n_points = n_points / 2

    println(max_depth)

    # Compute the order 
    global dx_list, err_list = compute_order(n_points, 5, n_gzones, x0, l, tfinal)

    global data =  DataFrame(dx = dx_list, err = err_list)

    data.logdx = log.(data.dx)
    data.logerr = log.(data.err)

    global model = lm(@formula(logerr ~ logdx), data)

    global a = coef(model)[2]
    global b = coef(model)[1]

    data.fitted_logerr = predict(model)
    data.fitted_err = exp.(data.fitted_logerr)

    # Plot 

    plot!(
        fig,
        data.dx, data.err,
        seriestype = :scatter,
        xscale = :log10, yscale = :log10,
        xlabel = "dx", ylabel = "global error", 
        marker = (:xcross, 6),
        markercolor = color[max_depth],
        label = @sprintf("%d level", max_depth),
        legend = :topleft
    )

    plot!(
        fig, 
        data.dx, data.fitted_err,
        xscale = :log10, yscale = :log10,
        linestyle = :dash,
        linecolor = color[max_depth],
        label = @sprintf("slope %.3f, r2 %.3f", a, r2(model)),
        legend = :topleft
    )

end



savefig(fig, file_out)








