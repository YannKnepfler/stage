
using Plots

include("level.jl")

function test(level::Union{Level, Nothing})

    if !isnothing(level)
        println("depth: ", level.depth, ", number of grid: ", length(level.grid_list), ", dx = ", level.grid_list[1].dx, ", x0: ", level.grid_list[1].x0)
        test(level.child)
    end

    return true
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

depth_list::Vector{Int64} = []
position_list::Vector{Float64} = []

output(base_level.grid_list[1], 1, depth_list, position_list)

solution::Vector{Vector{Any}} = [depth_list]
position::Vector{Vector{Any}} = [position_list]


ind = 0.0:0.001:1.0
exp_pos = collect(ind)
exp_sol = [init_cond(x) / 200.0 for x in ind]

push!(solution, exp_sol)
push!(position, exp_pos)


fig = plot(position, solution)

savefig(fig, "./test.png")
