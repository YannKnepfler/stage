using Plots


include("level.jl")
# Programme qui test les fonctionnalité des modules de amr
x0 = 0.0
L = 1.0
N_points_base = 32
N_gzones = 1
dx = (L - x0) / N_points_base
dt = dx



function test(level::Union{Level, Nothing})

    if !isnothing(level)
        println("depth: ", level.depth, ", number of grid: ", length(level.grid_list))
        test(level.child)
    end

    return true
end


base::Union{Level, Nothing}  = Level(N_points_base, N_gzones, x0, dt, dx)

initialize_level!(base, 2)


test(base)






