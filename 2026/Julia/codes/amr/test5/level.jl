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

function initialize_level!(level::Level)
    
    # Initialize all grid in the level
    for grid in level.grid_list
        initialize_grid!(grid)
    end

end

function regrid_level!(level::Level)

    # Check for the depth 
    if level.depth < max_depth

        # Definition of the local variables 
        new_child_list::Vector{Grid} = Vector{Grid}(undef, 0)
        
        # Definition of the new child list 
        for grid in level.grid_list

            # Update the local error 
            compute_local_error!(level)

            # Find the new children
            append!(new_child_list, regrid_one_grid!(grid))

        end

        # Set a new child level if needed 
        if !isnothing(level.child) || length(new_child_list) > 0
            if isnothing(level.child)
                # Create a new child level 
                level.child = Level(level)
            end
        end


        # Replace the child grid_list
        if length(new_child_list) > 0

            # Copy over the data from old grids
            copy_data_from_grid_list!(new_child_list, level.child.grid_list)

            level.child.grid_list = new_child_list

            # Reconnect the child grid list to the level grid list
            reconnect_child_list!(level.grid_list, level.child.grid_list)

            # Reconnect the child grid list to the child child grid list
            if !isnothing(level.child.child)
                reconnect_child_list!(level.child.grid_list, level.child.child.grid_list)
            end

        else
            if !isnothing(level.child)
                if !isnothing(level.child.child)
                    error("A level is about to be orphaned")
                end            
            end

            level.child = nothing 
        end
    end

    return nothing 

end

function evolve_level!(level::Level)

    # Advance the level 
    for grid in level.grid_list

        # Advance the solution
        advance_grid!(grid, advance)

        # Update the time 
        update_time_on_grid(grid::Grid)

    end

    # Advance the child level
    if !isnothing(level.child)
        evolve_level!(level.child)
    end

    # Update the overlapping control volumes 
    for grid in level.grid_list
        
        # Check if the grid has a parent grid to update
        if !isnothing(grid.parent)
            restrict_grid!(grid)
        end
    end

    # Update the ghost zones 
    compute_boundary_condition!(level)

    # Update the previous solution u_p 
    for grid in level.grid_list
        update_previous_solution!(grid)
    end

    # Adapt the child level 
    regrid_level!(level)

    return nothing
end

function compute_local_error!(level::Level)

    # Compute the local errror on each grid 
    for grid in level.grid_list
        compute_local_error_on_grid!(grid)
    end
    
    return nothing

end

function compute_boundary_condition!(level::Level)

    for grid in level.grid_list
        compute_boundary_condition_on_grid!(grid)
    end
    return nothing 
end

function output!(grid::Grid, depth::Int64, position_list::Vector{Float64}, solution_list::Vector{Float64}, depth_list::Vector{Int64})
    


    if length(grid.children) == 0 # No children 

        append!(position_list, grid.x[grid.box.inside_ind])
        append!(solution_list, grid.u[grid.box.inside_ind])
        append!(depth_list, [depth for ind in grid.box.inside_ind])

    else

        # Declaration of the local variables 
        istart = -1
        iend  = -1

        # Sort the child list 
        child_list::Vector{Grid} = []

        child_list = sort(grid.children, by = x -> x.box.lbnd)

        # println(length(grid.children))
        # println(length(child_list))

        # Initialization 
        istart = grid.box.inside_ind[1]
        iend = child_list[1].box.lbnd - 1

        # Check istart:iend
        if (istart in grid.box.inside_ind) && (iend in grid.box.inside_ind) && (istart <= iend)
            append!(position_list, grid.x[istart:iend])
            append!(solution_list, grid.u[istart:iend])
            append!(depth_list, [depth for ind in istart:iend])
        end


        # Loop on the child grid 
        for ind in eachindex(child_list)


            child::Grid = child_list[ind]

            output!(child, depth + 1, position_list, solution_list, depth_list)

            # Update the index 
            istart = child.box.ubnd + 1
            
            if ind == length(child_list)
                iend = grid.box.inside_ind[end]

            else
                iend = child_list[ind + 1].box.lbnd - 1
            end 

            # Check istart:iend
            if (istart in grid.box.inside_ind) && (iend in grid.box.inside_ind) && ( istart <= iend)
                append!(position_list, grid.x[istart:iend])
                append!(solution_list, grid.u[istart:iend])
                append!(depth_list, [depth for ind in istart:iend])
            end

  

        end
    end    

    return nothing
end

function update_dt!(base_level::Level)

    # Declaration of the local variables 
    level::Union{Level, Nothing} = base_level

    # Look for finest level 
    while !isnothing(level.child)
        level = level.child
    end

    # Compute the time step value 
    dt::Float64 = cfl * level.grid_list[1].dx

    # Update the time step value on every level 
    level = base_level
    while !isnothing(level)
        update_dt_on_level!(dt, level)
        level = level.child
    end

    return nothing 
end

function update_dt_on_level!(dt::Float64, level::Level)

    # Loop on the grid 
    for grid in level.grid_list
        update_dt_on_grid!(dt, grid)
    end

    return nothing 
end