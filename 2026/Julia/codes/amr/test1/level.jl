# Here the structure level and all of the functions related to this structur are defined

include("grid.jl")


mutable struct Level
    # depth of the level
    depth::Int64 

    # list of the grid of the level
    grid_list::Vector{Grid}    

    # Child level 
    child::Union{Level, Nothing}

    # Undefined constructor
    function Level()
        grid_list::Vector{Grid} = Vector{Grid}(undef, 0) 
        new(-1, grid_list, nothing)
    end
    
    # Constructor of the base level 
    function Level(N_points::Int64, N_gzones::Int64, x0::Float64, dt::Float64, dx::Float64)
        
        # The base level is at depth 1
        depth::Int64 = 1

        # Definition of the base grid 
        grid_list::Vector{Grid} = [Grid(N_points, N_gzones, x0, dt, dx)]

        # For now the base grid is the only grid 
        child::Union{Level, Nothing} = nothing 

        new(depth, grid_list, child)

    end

    # Constructor of a child level 
    function Level(parent::Level)

        # The depth is 1 plus the the depth of the parent level 
        depth = 1 + parent.depth

        # Defintiion of the grid_list
        grid_list::Vector{Grid} = []

        for grid in parent.grid_list
            append!(grid_list, grid.children)
        end

        # For now the level has no child 
        child::Union{Level, Nothing} = nothing 

        new(depth, grid_list, child)
    end
end


# Function that initialize the the value of the grids in the level 
function initialize_level!(base::Level, max_depth::Int64)
    for grid in base.grid_list
        initialize_grid!(grid)
    end

    level::Union{Level, Nothing} = base 
    eta::Float64 = 0.0
    
    while !isnothing(level) && (level.depth < max_depth)
        
        regrid_level!(level, eta, base.grid_list[1].N_points)

        if !isnothing(level.child)
            for grid in level.child.grid_list 
                initialize_grid!(grid)
            end
        end

        level = level.child
        
    end
end

# Function that regrid a level 
function regrid_level!(level::Level, eta::Float64, N_points_base::Int64)
    
    # Compute N and Update of eta
    N::Int64 = 0
    for grid in level.grid_list
        N += grid.N_points

        u::Vector{Float64} = grid.u[grid.N_gzones + 1 : grid.N_gzones + grid.N_points]

        gu_m, gu_p = compute_mean_grad(u)

        eta += gu_m + gu_p
    end

    # Compute de tolerance 
    tol = sqrt(eta / 2 * N)

    # Refine the grid 
    for grid in level.grid_list
        children::Union{Vector{Grid}, Nothing} = refine_grid(grid, tol, N_points_base)

        if !isnothing(children)
            grid.children = children

            if isnothing(level.child)
                level.child = Level(level)

            else
                append!(level.child.grid_list, children)
            end           
        end
    end

end

