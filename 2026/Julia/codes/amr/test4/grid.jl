# Here the structur Grid and the function related to this structur are defined 

include("box.jl")
include("functions.jl")
include("parameter.jl")

mutable struct Grid 

    # Number of points and ghost zones 
    n_points::Int64 
    n_gzones::Int64

    # Position of the left boundary of the grid
    x0::Float64 

    # Time step and cell size
    dt::Float64 
    dx::Float64

    # Solution 
    u::Vector{Float64}

    # Previous solution
    u_p::Vector{Float64}

    # Local error
    local_error::Vector{Float64}

    # Position 
    x::Vector{Float64}

    # Time of the solution 
    t::Float64

    # Box of the grid 
    box::Box 

    # Parent grid of this grid 
    parent::Union{Grid, Nothing}

    # List of child grids 
    children::Vector{Grid}


    # Constructor

    # Default constructor
    function Grid()
        new()
    end

    # Base grid constructor
    function Grid(n_points::Int64, n_gzones::Int64, x0::Float64, dt::Float64, dx::Float64)

        # Definition of the solution list 
        u::Vector{Float64} = Vector{Float64}(undef, n_points + 2*n_gzones)

        # Definition of the previous solution list 
        u_p::Vector{Float64} = Vector{Float64}(undef, n_points + 2*n_gzones)

        # Definition of the local error list 
        local_error::Vector{Float64} = fill(1.0 + 10.0 * amr_parameter, n_points + 2*n_gzones)

        # Definition of the position list 
        x::Vector{Float64} = Vector{Float64}(undef, n_points + 2*n_gzones)

        for i in 1 : n_points + 2*n_gzones
            ind = Float64(i - n_gzones)
            x[i] = x0 + (ind - 0.5) * dx
        end

        # Definition of the initial time
        t::Float64 = 0.0

        # Definition of the Box
        box::Box = Box(n_points, n_gzones)

        # The base grid does not have parent 
        parent::Union{Grid, Nothing} = nothing 

        # For nox the base grid does not have child grid
        children::Vector{Grid} = Vector{Grid}(undef, 0)

        new(n_points, n_gzones, x0, dt, dx, u, u_p, local_error, x, t, box, parent, children)
    end

    # Child grid constructor
    function Grid(parent::Grid, lbnd::Int64, ubnd::Int64)

        # Definition of the child box 
        box::Box = Box(parent.box, lbnd, ubnd)

        # Number of points and ghost zones 
        n_points::Int64 = box.n_points 
        n_gzones::Int64 = box.n_gzones

        # Definition of the left border position 
        x0::Float64 = parent.x[lbnd] - 0.5 * parent.dx

        # Definition of the time step and the cell size 
        dt::Float64 = parent.dt
        dx::Float64 = parent.dx * 0.5

        # Definition of the solution list 
        u::Vector{Float64} = Vector{Float64}(undef, n_points + 2*n_gzones)

        # Definition of the previous solution list 
        u_p::Vector{Float64} = Vector{Float64}(undef, n_points + 2*n_gzones)

        # Definition of the local error list 
        local_error::Vector{Float64} = Vector{Float64}(undef, n_points + 2*n_gzones)

        # Definition of the position list 
        x::Vector{Float64} = Vector{Float64}(undef, n_points + 2*n_gzones)

        for i in 1 : n_points + 2*n_gzones
            ind = Float64(i - n_gzones)
            x[i] = x0 + (ind - 0.5) * dx
        end

        # Defintion of the starting time of the grid 
        t::Float64 = parent.t

        # For now the child grid does not have child grid itself
        children::Vector{Grid} = Vector{Grid}(undef, 0)

        new(n_points, n_gzones, x0, dt, dx, u, u_p, local_error, x, t, box, parent, children)
    end

    # Child grid constructor
    function Grid(parent::Grid, child_box::Box)

        # Number of points and ghost zones 
        n_points::Int64 = child_box.n_points 
        n_gzones::Int64 = child_box.n_gzones

        # Definition of x0 
        x0::Float64 = parent.x[child_box.lbnd] - 0.5 * parent.dx

        # Definition of the time step and the cell size 
        dt::Float64 = parent.dt
        dx::Float64 = parent.dx * 0.5

        # Definition of the solution list 
        u::Vector{Float64} = Vector{Float64}(undef, n_points + 2*n_gzones)

        # Definition of the previous solution list 
        u_p::Vector{Float64} = Vector{Float64}(undef, n_points + 2*n_gzones)

        # Definition of the local error list 
        local_error::Vector{Float64} = Vector{Float64}(undef, n_points + 2*n_gzones)

        # Defintion of the position list 
        x::Vector{Float64} = Vector{Float64}(undef, n_points + 2*n_gzones)

        for i in 1 : n_points + 2*n_gzones
            ind = Float64(i - n_gzones)
            x[i] = x0 + (ind - 0.5) * dx
        end

        # Defintion of the starting time of the grid 
        t::Float64 = parent.t

        # For now the child grid does not have child grid itself
        children::Vector{Grid} = Vector{Grid}(undef, 0)

        new(n_points, n_gzones, x0, dt, dx, u, u_p, local_error, x, t, child_box, parent, children)

    end
end

function initialize_grid!(grid::Grid)

    # Initialization of the solution list and the previous solution list
    for ind in 1 : grid.n_points + 2 * grid.n_gzones
        grid.u[ind]   = exact_solution(0.0, grid.x[ind])
        grid.u_p[ind] = exact_solution(0.0, grid.x[ind])
    end
    
    return nothing
end

function advance_grid!(grid::Grid, advance::Function)

    # Definition of the local variable
    new_u::Vector{Float64} = advance(grid)

    # Update of the solution 
    grid.u = new_u

    return nothing
end

function apply_error_estimator!(grid::Grid, err_estimator!::Function)

    # Update the local error 
    err_estimator!(grid, epsilon)

    return nothing
end

function add_children_to_errors!(grid::Grid)

    # Update the local error where the grid has grandchildren
    for child in grid.children
        for greatchild in child.children 

            child_start = greatchild.box.lbnd 
            grid_start = child.box.lbnd + div(child_start - (child.n_gzones + 1), 2)

            child_end = greatchild.box.ubnd 
            grid_end = child.box.lbnd + div(child_end - (child.n_gzones + 1), 2)

            grid.local_error[grid_start : grid_end] .= 1.0 + 10.0 * amr_parameter

        end
    end

    return nothing
end

function swap_current_and_previous_solution!(grid::Grid)

    for ind in 1 : grid.n_points  + 2 * grid.n_gzones
        grid.u[ind] = grid.u_p[ind]
    end

    return nothing
end

function regrid_one_grid(grid::Grid)

    # Declaration of the local variables 
    new_children::Vector{Grid} = Vector{Grid}(undef, 0)

    new_box_list::Union{Vector{Box}, Nothing} = find_error_boxes(grid)

    # Check if there are new children
    if isnothing(new_box_list) 
    else 
        # Set up the new children 
        for child_box in new_box_list
            new_child ::Grid = Grid(grid, child_box)

            # Prolongate the new child 
            prolongate_grid!(new_child)

            # Compute the boundary conditions
            compute_boundary_condition_on_grid!(new_child)

            # Add the new child to the new children list 
            push!(new_children, new_child)
        end
    end

    return new_children
end

function find_error_boxes(grid::Grid)

    # Declaration of the local variables 
    child_box_list::Union{Vector{Box}, Nothing} = nothing

    error_tmp::Vector{Bool} = fill(false, grid.n_points + 2 * grid.n_gzones)
    in_error::Vector{Bool}  = fill(false, grid.n_points + 2 * grid.n_gzones)

    # Initialization of error_tmp
    for ind in grid.box.inside_ind
        if grid.local_error[ind] > amr_parameter
            error_tmp[ind] = true
        end
    end

    # Initialization of in_error 
    for ind in grid.box.inside_ind
        istart::Int64 = max(grid.box.inside_ind[1], ind - amr_buffer_width)
        iend::Int64   = min(grid.box.inside_ind[end], ind + amr_buffer_width)

        if any(error_tmp[istart : iend])
            in_error[ind] = true
        end
    end

    # Definition of the new child box list 
    istart = 0
    iend = -1

    for ind in grid.box.inside_ind        
        if istart < 1 # If the start of the box has not been found
            if in_error[ind]
                istart = ind 

            end
        else # The start of the box has been found but the not the finish 
            if !(in_error[ind]) 
                iend = ind - 1

                # Compute the box and add the new box to the child_box_list
                if isnothing(child_box_list) 
                    child_box_list = [Box(grid.box, istart, iend)]
                else
                    push!(child_box_list, Box(grid.box, istart, iend))
                end

                # Update istart and iend
                istart = 0
                iend = -1

            end
        end
    end

    # Check if there is not a box to the right boundary
    if istart > 0
        iend = grid.box.inside_ind[end]

        # Compute the box and add the new box to the child_box_list
        if isnothing(child_box_list) 
            child_box_list = [Box(grid.box, istart, iend)]
        else
            push!(child_box_list, Box(grid.box, istart, iend))
        end

    end



    return child_box_list

end

function prolongate_grid!(child::Grid)

    # Declaration of the local variables 

    istart = child.box.lbnd - div(child.n_gzones, 2)
    iend = child.box.ubnd + div(child.n_gzones, 2)

    # Check if the number of ghost zones is odd
    if child.n_gzones % 2 == 1

        # Left most ghost zone 
        parent_i = istart - 1

        # Compute the slope 
        slope = compute_slope(child.parent, parent_i)

        # Update the ghost zone
        child.u[1] = child.parent.u[parent_i] + 0.25 * slope

        # Right most ghost zone 
        parent_i = iend + 1

        # Compute the slope 
        slope = compute_slope(child.parent, parent_i)

        # Update the ghost zone
        child.u[end] = child.parent.u[parent_i] - 0.25 * slope


    end

    # Update the rest
    for parent_i in istart : iend

        # Identify the child index 
        child_i = (parent_i - child.box.lbnd) * 2 + child.n_gzones + 1
        
        # Compute the slope 
        slope = compute_slope(child.parent, parent_i)
        

        # Update the value 
        child.u[child_i]     = child.parent.u[parent_i] - 0.25 * slope
        child.u[child_i + 1] = child.parent.u[parent_i] + 0.25 * slope

    end


    return nothing 
end

function compute_slope(grid::Grid, ind::Int64)

    # Declaration of the local variables 
    slope::Float64 = 0.0

    # Check for border index 
    if ind == grid.box.lgzone_ind[1]
        slope = grid.u[ind + 1] - grid.u[ind]

    elseif ind == grid.box.rgzone_ind[end]
        slope = grid.u[ind] - grid.u[ind - 1]

    else
        slope_l::Float64 = grid.u[ind] - grid.u[ind - 1]
        slope_r::Float64 = grid.u[ind + 1] - grid.u[ind]

        # Apply minmod slope limiter 
        if slope_l * slope_r > 0.0
            if abs(slope_l) < abs(slope_r)
                slope = slope_l

            else
                slope = slope_r

            end
        end
    end

    return slope
end

function copy_data_from_grid_list!(new_grid_list::Vector{Grid}, old_grid_list::Vector{Grid})

    for new_grid in new_grid_list
        for old_grid in old_grid_list
            intersection = grid_intersection(new_grid, old_grid)

            # Check if the grid intersect 
            if !isnothing(intersection)

                istart, iend, offset = intersection 

                # Copy the data from the old grid on the intersection 
                new_grid.u[istart : iend] = old_grid.u[(istart + offset) : (iend + offset)]
            end
        end
    end
end

function grid_intersection(new_grid::Grid, old_grid::Grid)

    # Definition of the local variables 
    dx = new_grid.dx
    new_start = new_grid.x0
    new_end   = new_grid.x0 + new_grid.n_points * dx
    old_start = old_grid.x0
    old_end   = old_grid.x0 + old_grid.n_points * dx

    if (new_start + 0.5 * dx < old_end) && (new_end > old_start + 0.5 * dx )

        intersection_start = max(new_start, old_start)
        intersection_end = min(new_end, old_end)
        istart = new_grid.box.inside_ind[round(Int64, (intersection_start - new_start) / dx) + 1]
        iend = new_grid.box.inside_ind[round(Int64, (intersection_end - new_start) / dx)]
        offset = round(Int64, (new_start - old_start) / dx)

        return istart, iend, offset 
    else
        return nothing 
    end

end

function reconnect_child_list!(parent_list::Vector{Grid}, child_list::Vector{Grid})

    # Declaration of the local variables 
    parent_of_child::Union{Grid, Nothing} = nothing 

    # Clear the parent_children list 
    for parent in parent_list
        parent.children = Vector{Grid}(undef, 0)
    end

    # Loop on the child grid 
    for child in child_list

        # Identify the parent grid
        parent_of_child = nothing 

        for parent in parent_list
            if grid_child_of_parent(child, parent)
                parent_of_child = parent 
            end
        end

        if isnothing(parent_of_child)
            error("Failed to reconnect a child grid")
        end

        # Update the child grid 
        child.parent = parent_of_child
        correct_box(child)

        # Update the parent grid 
        push!(parent_of_child.children, child)
    end
end

function grid_child_of_parent(child::Grid, parent::Grid)

    # Declaration of the local variables 
    child_start = child.x0 - 0.5 * child.dx 
    child_end = child.x0 + child.n_points * child.dx  + 0.5 * child.dx
    parent_start = parent.x0 - 0.5 * parent.dx 
    parent_end = parent.x0 + parent.n_points * parent.dx  + 0.5 * parent.dx

    return (child_start > parent_start) && (child_end < parent_end)
end

function correct_box(child::Grid)

    # Compute the offset 
    offset = round(Int64, (child.x0 - child.parent.x0) / child.parent.dx)
    
    # Update the child box 
    child.box.lbnd = child.parent.box.inside_ind[offset + 1]
    child.box.ubnd = child.box.lbnd + div(child.n_points, 2) - 1

    if !((child.box.lbnd in child.parent.box.inside_ind) && (child.box.ubnd  in child.parent.box.inside_ind))
        println("lbnd = ", child.box.lbnd, ", ubnd = ", child.box.ubnd)
    end
end

function restrict_grid!(child::Grid)


    # Delclaration of the local variables 
    istart = child.box.lbnd
    iend = child.box.ubnd

    # Correction of the overlapping control volumes
    for parent_i in istart : iend

        # Identify the first child control volume 
        child_i = (parent_i - child.box.lbnd) * 2 + child.n_gzones + 1

        # Update child local error 
        child.local_error[child_i] = abs((child.parent.u[parent_i] - 0.5 * (child.u[child_i] + child.u[child_i + 1])))
        child.local_error[child_i + 1] = child.local_error[child_i]
        # Update the parent control volume 
        child.parent.u[parent_i] = 0.5 * (child.u[child_i] + child.u[child_i + 1])

    end

    return nothing 
end

function update_previous_solution!(grid::Grid)

    # Copy u on u_p
    for ind in 1 : grid.n_points + 2 * grid.n_gzones
        grid.u_p[ind] = grid.u[ind]
    end

    return nothing 
end

function update_time_on_grid(grid::Grid)
    grid.t += grid.dt
end

function compute_boundary_condition_on_grid!(grid::Grid)

    # Left boundary
    if grid.box.bbox[1]
        for ind in grid.box.lgzone_ind
            grid.u[ind] = exact_solution(grid.t, grid.x[ind])
        end

    else # Prolongate the left ghost zone 
        
        # Declaration of the local variables 
        istart = grid.box.lbnd - div(grid.n_gzones, 2)
        iend = grid.box.lbnd - 1

        # Check if the number of ghostzone is odd 
        if grid.n_gzones % 2 == 1
            parent_i = istart - 1

            # Compute the slope 
            slope = compute_slope(grid.parent, parent_i)

            # Update the ghost zone 
            grid.u[1] = grid.parent.u[parent_i] + 0.25 * slope
        end

        # Update the rest 
        for parent_i in istart : iend

            # Identify the child control volume index 
            child_i = (parent_i - child.box.lbnd) * 2 + child.n_gzones + 1

            # Compute the slope 
            slope = compute_slope(grid.parent, parent_i)

            # Update the ghost zone 
            grid.u[child_i]     = grid.parent[parent_i] - 0.25 * slope
            grid.u[child_i + 1] = grid.parent[parent_i] + 0.25 * slope
        end
        
    end

    # Right boundary 
    if grid.box.bbox[2]
        for ind in grid.box.rgzone_ind
            grid.u[ind] = exact_solution(grid.t, grid.x[ind])
        end

    else # Prolongate the right ghost zone 

        # Declaration of the local variables
        istart = grid.box.ubnd + 1
        iend = grid.box.ubnd + div(grid.n_gzones, 2)

        # Check if the number of ghost zoe is odd 
        if grid.n_gzones % 2 == 1
            parent_i = iend + 1 

            # Compute the slope 
            slope = compute_slope(grid.parent, parent_i)
            
            # Update the ghost zone 
            grid.u[end] = grid.parent.u[parent_i] - 0.25 * slope
        end

        # Update the rest 
        for parent_i in istart : iend

            # Identify the child control volume index 
            child_i = (parent_i - child.box.lbnd) * 2 + child.n_gzones + 1

            # Compute the slope 
            slope = compute_slope(grid.parent, parent_i)

            # Update the ghost zone 
            grid.u[child_i]     = grid.parent[parent_i] - 0.25 * slope
            grid.u[child_i + 1] = grid.parent[parent_i] + 0.25 * slope
        end


    end

    return nothing
end

function set_initial_local_error_on_grid!(grid::Grid)
    
    for ind in grid.box.inside_ind
        grid.local_error[ind] = abs(grid.u[ind] - grid.u[ind - 1])
    end

    return nothing

end

function update_dt_on_grid!(dt::Float64, grid::Grid)
    grid.dt = dt
    return nothing 
end

