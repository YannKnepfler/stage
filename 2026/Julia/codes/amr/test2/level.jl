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
    function Level(N_points::Int64, N_gzones::Int64, x0::Float64, dt::Float64, dx::Float64)

        # The base level is at level 1
        depth::Int64 = 1

        # The base level only contains the base grid
        grid_list::Vector{Grid} =  [Grid(N_points, N_gzones, x0, dt, dx)]

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

# Function that initialize the grids 
function initialize!(base::Level, max_depth::Int64)

    # Initialization of the grid at this level 
    for grid in base.grid_list
        initialize_grid!(grid) 
    end

    level::Union{Level, Nothing} = base 
    eta::Float64 = 0.0

    while !isnothing(level) && (level.depth < max_depth)
        
        child_level::Union{Level, Nothing} = initialize_level(level, eta, base.grid_list[1].N_points)

        if !isnothing(child_level)
            for grid in child_level.grid_list
                initialize_grid!(grid)
            end
        end

        level.child = child_level
        level = child_level
        
    end
end

# Function that initialize a level 
function initialize_level(level::Level, eta::Float64, N_points_base::Int64)
    println("depth = ", level.depth)
    # Defintion of the local variables 
    child_level::Union{Level, Nothing} = nothing
    N::Int64 = 0
    tol::Float64 = 0.0

    # Compute N and update eta
    for grid in level.grid_list
        # Compute N
        N += grid.N_points

        # Compute eta
        u::Vector{Float64} = grid.u[grid.box.without_bz]

        # println("compute tol")
        gu_m, gu_p = compute_mean_grad(u)
        
        eta += gu_m + gu_p 
    end


    # Compute the tolerance 
    tol = sqrt(eta / (2 * N))

    # println("tol = ", tol)

    if level.depth == 1
        tol = 1.0 * tol
    end


    # Refine the grids 
    for grid in level.grid_list
        children::Union{Vector{Grid}, Nothing} = refine_grid(grid, tol, N_points_base)

        if !isnothing(children)

            grid.children = children

            if isnothing(child_level)
                child_level = Level(level)

            end 
            
            append!(child_level.grid_list, children)
        end
    end

    return child_level
end

# Function that adapt the level structur
function adapt!(base::Level, max_depth::Int64) 

    # Definition of the local varaibles
    level::Union{Level, Nothing} = base
    eta::Float64 = 0.0

    # Loop on the child level 
    while !isnothing(level) && (level.depth < max_depth)

        # println("depth = ", level.depth)
        adapt_level!(level, eta, base.grid_list[1].N_points)

        
        level = level.child

    end
end

# Function that adpat a level 
function adapt_level!(level::Level, eta::Float64, N_points_base::Int64)
    

    println("depth = ", level.depth)
    # Defintion of the local variables 
    N::Int64 = 0
    tol::Float64 = 0.0

    # Compute N and update eta
    for grid in level.grid_list
        # Compute N
        N += grid.N_points

        # Compute eta
        u::Vector{Float64} = grid.u[grid.box.without_bz]

        # println("compute tol")
        gu_m, gu_p = compute_mean_grad(u)

        eta += gu_m + gu_p 
        
    end


    # Compute the tolerance 
    tol = sqrt(eta / (2 * N))

    if level.depth == 1
        tol = 1.0 * tol
    end


    # println("tol = ", tol)


    # Clear the old children list of the child grid 
    if !isnothing(level.child)
        level.child.grid_list = Vector{Grid}(undef, 0)
    end

    # Refine the grids
    for grid in level.grid_list
        
        # Compute the new children of the grid
        new_children::Union{Vector{Grid}, Nothing} = refine_grid(grid, tol, N_points_base)

        if isnothing(new_children)
            grid.children = Vector{Grid}(undef, 0)
        else

            # Prolongate the new child 
            prolongate_children!(grid, new_children)


            # Add the new children to the parent grid
            grid.children = new_children


            # Add the new children list to the child level 
            if isnothing(level.child)
                level.child = Level(level)

            end

            append!(level.child.grid_list, new_children)           

        end


    end

    # Delete the child level if it has an empty grid list 
    if !isnothing(level.child)
        if length(level.child.grid_list) == 0 
            level.child = nothing
        end
    end


end



function output!(grid::Grid, depth::Int64, depth_list::Vector{Int64}, position_list::Vector{Float64}, solution_list::Vector{Float64})



    if length(grid.children) == 0

        i_start = grid.box.inside_ind[1]
        i_end = grid.box.inside_ind[end]

        append!(depth_list, [depth for i in i_start : i_end])
        append!(position_list, grid.x[i_start : i_end])
        append!(solution_list, grid.u[i_start : i_end])

    else

        # Definition of the boundary list 
        boundaries::Vector{Tuple{Int64, Int64}} = []

        for child in grid.children
            push!(boundaries, (child.box.lbnd, child.box.ubnd))
        end

        # Sorting the lists
        sort!(boundaries, by = x -> x[1])
        children = sort(grid.children, by = x -> x.box.lbnd)


        # Initialization 
        if boundaries[1][1] != grid.box.inside_ind[1]
            i_start = grid.box.inside_ind[1]
            i_end = boundaries[1][1] - 1

            append!(depth_list, [depth for i in i_start : i_end])
            append!(position_list, grid.x[ i_start : i_end ])
            append!(solution_list, grid.u[i_start : i_end])
        end

    
        # Loop on the grid child       

        for i in eachindex(boundaries)
            child::Grid = children[i]
            
            output!(child, depth + 1, depth_list, position_list, solution_list)
            if i < length(boundaries)
                i_start = boundaries[i][2] + 1
                i_end = boundaries[i+1][1] - 1

                append!(depth_list, [depth for i in i_start : i_end])
                append!(position_list, grid.x[ i_start : i_end ])
                append!(solution_list, grid.u[i_start : i_end])
            end
            
        

        
        end

        # Finalization 
        if boundaries[end][2] != grid.box.inside_ind[end]
            i_start = boundaries[end][2] + 1
            i_end = grid.box.inside_ind[end]

            append!(depth_list, [depth for i in i_start : i_end])
            append!(position_list, grid.x[ i_start : i_end ])
            append!(solution_list, grid.u[i_start : i_end])
        end
    end
    

    return depth_list, position_list, solution_list
end