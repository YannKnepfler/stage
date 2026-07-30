
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


N_points::Int64 = 16
N_gzones::Int64 = 0

x0::Float64 = 0.0
l::Float64 = 1.0
dx::Float64 = (l - x0) / N_points
dt::Float64 = dx 

# test

base_level::Union{Level, Nothing} = Level(N_points, N_gzones, x0, dt, dx)

initialize!(base_level, 5)

test(base_level)

move!(base_level, 0.05)

adapt!(base_level, 5)

println("")
println("======================== MOVED ========================")
println("")
# test(base_level)

println( length(base_level.child.child.child.child.grid_list), " ", base_level.child.child.child.child.depth )



depth_list::Vector{Int64} = []
position_list::Vector{Float64} = []
solution_list::Vector{Float64} = []

output!(base_level.grid_list[1], 1, depth_list, position_list, solution_list)


# ind = 0.0:0.001:1.0
# exp_pos = collect(ind)
# exp_sol = [init_cond(x - 0.25) for x in ind]


fig = plot( 
    position_list, solution_list,
    seriestype = :line,
    marker = :xcross,
    label = nothing,
    title = "test",
    linecolor = :orange,
    markercolor = :orange,
    x_lims = [0, 1],
    y_lims = [0, 1],
    legend = :topright
)

plot!( 
    twinx(fig),
    position_list, depth_list,
    seriestype = :line,
    marker = :xcross,
    label = nothing,
    x_lims = [0, 1],
    y_lims = [0, 8],
    ylabel = "depth"
)


savefig(fig, "../image/criteria1/test1/test_mv_5.png")

