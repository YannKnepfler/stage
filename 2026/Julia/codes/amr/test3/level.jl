# Here the structur level and the function related to this structur are defined 

include("grid.jl")

mutable struct Level


    # Refinment depth of the level 
    depth::Int64

    # List of t=the grid at this level 
    grid_list::Vector{Grid}

    # Parent level 
    parent::Union{Level, Nothing}

    # Child level 
    child::Union{Level, Nothing}

    # Constructor 

    # Default constructor 
    function Level()
        new()
    end

    # Base level constructor 
    function Level(n_points::Int64, n_gzones::Int64, x0::Float64, dt::Float64, dx::Float64)

        # The base level is at level 1
        depth::Int64 = 1

        # The base level only contains the base grid
        grid_list::Vector{Grid} =  [Grid(n_points, n_gzones, x0, dt, dx)]

        # The base level does not have a parent level 
        parent::Union{Level, Nothing} = nothing 

        # For now the level does not have a child level 
        child::Union{Level, Nothing} = nothing 

        new(depth, grid_list, parent, child)
    end

    # Child level constructor 
    function Level(parent::Level)

        # Definition of the depth 
        depth::Int64 = parent.depth + 1

        # Defintion of the list of grid 
        grid_list::Vector{Grid} =Vector{Grid}(undef, 0)

        # for grid in parent.grid_list
        #     append!(grid_list, grid.children)
        # end

        # For now the level does not have a child level 
        child::Union{Level, Nothing} = nothing 

        new(depth, grid_list, parent, child)
    end

end

function initialize_base_level!(base::Level)
    # Initialization of the base grid 
    initialize_grid!(base.grid_list[1])
end