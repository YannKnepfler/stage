
using Plots

include("level.jl")

function test(level::Union{Level, Nothing})

    if !isnothing(level)
        println("depth: ", level.depth, ", number of grid: ", length(level.grid_list), ", dx = ", level.grid_list[1].dx)
        test(level.child)
    end

    return true
end


function move!(level::Union{Level, Nothing}, gap::Float64)
    
    # Check if the there is a level 
    if !isnothing(level)

        for grid in level.grid_list

            for i in 1 : grid.N_points + 2 * grid.N_gzones
                grid.u[i] = init_cond(grid.x[i] - gap )
            end

            move!(level.child, gap)

        end

    end


end

function vizulisation(base_level::Level, N_ite::Int64, gap::Float64, d_max::Int64, dir_out::String)
    
    # Initialization 
    initialize!(base_level, d_max)

    depth_list::Vector{Int64} = []
    position_list::Vector{Float64} = []
    solution_list::Vector{Float64} = []

    output!(base_level.grid_list[1], 1, depth_list, position_list, solution_list)

    # Plot 
    fig = plot(
        position_list, solution_list,
        title = "iteration: 0",
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
        ylims = [0, d_max],
        ylabel = "depth"
    )

    file_out::String = dir_out * "test_ite_0.png" 

    savefig(fig, file_out)

    for ite in 1:N_ite

        # Move and adapt
        move!(base_level, ite * gap)

        adapt!(base_level, d_max)

        if ite % 10 == 0

            depth_list = []
            position_list = []
            solution_list = []

            output!(base_level.grid_list[1], 1, depth_list, position_list, solution_list)

            # Plot
            fig = plot(
                position_list, solution_list,
                title = "iteration: " * string(ite),
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
                ylims = [0, d_max],
                ylabel = "depth"
            )

            file_out = dir_out * "test_ite_" * string(ite) * ".png"
            
            savefig(fig, file_out)
        end



    end

    return true
end


N_points::Int64 = 16
N_gzones::Int64 = 2

x0::Float64 = 0.0
l::Float64 = 1.0
dx::Float64 = (l - x0) / N_points
dt::Float64 = dx 

# test

base_level::Union{Level, Nothing} = Level(N_points, N_gzones, x0, dt, dx)

N_ite::Int64 = 10
gap::Float64 = 0.005
d_max::Int64 = 8
dir_out::String = "../image/criteria1/test5/"

vizulisation(base_level, N_ite, gap, d_max, dir_out)

# initialize!(base_level, 8)

# test(base_level)

# println("")
# println("======================== MOVED 1 ========================")
# println("")

# move!(base_level, 0.01 / 2.0)

# adapt!(base_level, 8)
# test(base_level)

# println("")
# println("======================== MOVED 2 ========================")
# println("")

# move!(base_level, 0.01 / 2.0)

# adapt!(base_level, 8)
# test(base_level)





# depth_list::Vector{Int64} = []
# position_list::Vector{Float64} = []
# solution_list::Vector{Float64} = []

# output!(base_level.grid_list[1], 1, depth_list, position_list, solution_list)


# ind = 0.0:0.001:1.0
# exp_pos = collect(ind)
# exp_sol = [init_cond(x - 0.25) for x in ind]


# fig = plot( 
#     position_list, solution_list,
#     seriestype = :line,
#     marker = :xcross,
#     label = nothing,
#     title = "Basr number of cell: " * string(N_points),
#     linecolor = :orange,
#     markercolor = :orange,
#     x_lims = [0, 1],
#     y_lims = [0, 1],
#     legend = :topright
# )

# plot!( 
#     twinx(fig),
#     position_list, depth_list,
#     seriestype = :line,
#     marker = :xcross,
#     label = nothing,
#     x_lims = [0, 1],
#     y_lims = [0, 8],
#     ylabel = "depth"
# )


# savefig(fig, "../image/criteria1/test2/coef0.0/test_crit1_" * string(N_points) * ".png")

