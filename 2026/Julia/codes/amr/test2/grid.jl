# Here the structur Grid and the function related to this structur are defined 

include("box.jl")
include("functions.jl")

mutable struct Grid 

    # Number of points and ghost zones 
    N_points::Int64 
    N_gzones::Int64

    # Position of the left boundary of the grid
    x0::Float64 

    # Time step and cell size
    dt::Float64 
    dx::Float64

    # Solution 
    u::Vector{Float64}

    # Position 
    x::Vector{Float64}

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
    function Grid(N_points::Int64, N_gzones::Int64, x0::Float64, dt::Float64, dx::Float64)

        # Definition of the solution list 
        u::Vector{Float64} = Vector{Float64}(undef, N_points + 2*N_gzones)

        # Definition of the position list 
        x::Vector{Float64} = Vector{Float64}(undef, N_points + 2*N_gzones)

        for i in 1 : N_points + 2*N_gzones
            ind = Float64(i - N_gzones)
            x[i] = x0 + (ind - 0.5) * dx
        end

        # Definition of the Box
        box::Box = Box(N_points, N_gzones)

        # The base grid does not have parent 
        parent::Union{Grid, Nothing} = nothing 

        # For nox the base grid does not have child grid
        children::Vector{Grid} = Vector{Grid}(undef, 0)

        new(N_points, N_gzones, x0, dt, dx, u, x, box, parent, children)
    end

    # Child grid Constructor
    function Grid(parent::Grid, lbnd::Int64, ubnd::Int64)

        # Definition of the child box 
        box::Box = Box(parent.box, lbnd, ubnd)

        # Number of points and ghost zones 
        N_points::Int64 = box.N_points 
        N_gzones::Int64 = box.N_gzones

        # Definition of the left border position 
        x0::Float64 = parent.x[box.lbnd] - 0.5 * parent.dx

        # Definition of the time step and the cell size 
        dt::Float64 = parent.dt
        dx::Float64 = parent.dx * 0.5

        # Definition of the solution list 
        u::Vector{Float64} = Vector{Float64}(undef, N_points + 2*N_gzones)

        # Definition of the position list 
        x::Vector{Float64} = Vector{Float64}(undef, N_points + 2*N_gzones)

        for i in 1 : N_points + 2*N_gzones
            ind = Float64(i - N_gzones)
            x[i] = x0 + (ind - 0.5) * dx
        end

        # For now the child grid does not have child grid itself
        children::Vector{Grid} = Vector{Grid}(undef, 0)

        new(N_points, N_gzones, x0, dt, dx, u, x, box, parent, children)
    end 
end

# Functions that initialize the grid 
function initialize_grid!(grid::Grid)
    
    for i in 1 : grid.N_points + 2 * grid.N_gzones
        grid.u[i] = init_cond(grid.x[i])
    end
end


# Function that compute the mean grad of each half of a vector 
function compute_mean_grad(u::Vector{Float64})

    # Definition of the mean gradient of the two half of the grid 
    gu_m::Float64 = 0.0
    gu_p::Float64 = 0.0

    N_points::Int64 = length(u)

    # Compute the first half 
    for i in 2 : div(N_points, 2)
        gu_m += (u[i] - u[i-1])^2
    end

    gu_m = 2 * gu_m / N_points

    # Compute the second half 
    for i in div(N_points, 2) + 1 : N_points
        gu_p += (u[i] - u[i-1])^2
    end

    gu_p = 2 * gu_p / N_points
    # println("gu_m = ", gu_m, ", gu_p = ", gu_p)
    return gu_m, gu_p
end
    
# Function that return a list of children of parent grid according to a tolerance 
function refine_grid(parent::Grid, tol::Float64, N_point_base::Int64)

    # Definition of the local variable 
    children::Vector{Grid}  = Vector{Grid}(undef, 0)
    boundaries::Vector{Tuple{Int64, Int64}} = []

    N_int::Int64 = length(parent.box.without_bz)

    N_patches::Int64 = div(N_int, N_point_base)
    
    # Loop on the patches 
    for i in 1:N_patches

        # Breakdown of the grid 
        index = parent.box.without_bz[ (i-1) * N_point_base + 1 : i * N_point_base  ]
        u = parent.u[index]

        lbnd, ubnd = refine_bnd(u, index, tol)

        # Update the boundary list 
        if (lbnd > 0) && (ubnd > 0)
            if length(boundaries) == 0
                push!(boundaries, (lbnd, ubnd))

            else
                if lbnd - 1 == boundaries[end][2] 
                    lbnd, _ = pop!(boundaries)

                end

                push!(boundaries, (lbnd, ubnd))
            end
        end
    end

    # Generate the new child grid list 
    if length(boundaries) == 0
        return nothing

    else
        for bnd in boundaries
            lbnd, ubnd = bnd 

            push!(children, Grid(parent, lbnd, ubnd)) 
        end  
    end

    return children 
end

# Function that return the boundaries between wich to refine 
function refine_bnd(u::Vector{Float64}, index::UnitRange{Int64}, tol::Float64)

    # Definition of the local variables
    lbnd::Int64 = -1
    ubnd::Int64 = -1

    N_points = length(u)
    # Compute the mean gradient
    gu_m, gu_p = compute_mean_grad(u)

    # Compare the gradient to the tolerance 
    test = [(gu_m >= tol), (gu_p >= tol)]


    if test == [true, true]
        lbnd = index[1]
        ubnd = index[end]

    elseif test == [true, false]
        lbnd = index[1]
        ubnd = index[div(N_points, 2)]

    elseif test == [false, true]
        lbnd = index[div(N_points, 2) + 1]
        ubnd = index[end]

    end
    
    return lbnd, ubnd
end

function prolongate_children!(parent::Grid, new_children::Vector{Grid})

    # Declaration of the local variables
    old_child::Union{Grid, Nothing} = nothing
    parent_ind::Int64 = 0
    slope::Float64 = 0.0


    # Loop on the new children 
    for new_child in new_children 

        parent_ind = new_child.box.lbnd

        # Loop on the inside control volumes 
        while parent_ind <= new_child.box.ubnd

            # Identifying the first cell index 
            new_ind = new_child.box.inside_ind[ (parent_ind - new_child.box.lbnd) * 2 + 1 ]

            # Check in wich child grid the cell was refined 
            old_child = is_refined(parent, parent_ind)

            if isnothing(old_child) # The parent control volume was not refined

                # Compute the slope 
                slope = compute_slope(parent, parent_ind)

                # Prolongate the new child
                new_child.u[new_ind]     = parent.u[parent_ind] - 0.25 * slope
                new_child.u[new_ind + 1] = parent.u[parent_ind] + 0.25 * slope

                # Update the parent index 
                parent_ind += 1

            else # The parent controle volume was refined

                # Compute the gap between the two grid
                gap::Int64 = 2 * (new_child.box.lbnd - old_child.box.lbnd)

                # Transfer the overlapping child grid of the old grid to the new one
                new_child_children::Union{Vector{Grid}, Nothing} = child_transfer!(new_child, old_child, gap)

                if !isnothing(new_child_children)
                    append!(new_child.children, new_child_children)
                end
                
                # Prolongate the child
                for old_ind in old_child.box.inside_ind
                    if old_ind - gap in new_child.box.inside_ind
                        new_child.u[old_ind - gap] = old_child.u[old_ind]
                    end
                end

                # Update the parent index 
                parent_ind = min(old_child.box.ubnd, new_child.box.ubnd) + 1

            end

        end

        # Prolongate the ghost zones
        set_ghost_zones!(parent, new_child)


    end
end

# Function that return the child in wich the control volume at index ind is refine or nothing if the controle is not refined 
function is_refined(grid::Grid, ind::Int64)

    result::Union{Grid, Nothing} = nothing 

    # Check in which child grid the refined controle volume is contained 
    for child in grid.children
        if child.box.lbnd <= ind <= child.box.ubnd 
            result = child
        end
    end

    return result
end

function compute_slope(grid::Grid, ind::Int64)

    # Declaration of the local variables 
    slope::Float64 = 0.0

    # Check for border index 
    if ind == grid.box.inside_ind[1]
        slope = grid.u[ind + 1] - grid.u[ind]

    elseif ind == grid.box.inside_ind[end]
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

function child_transfer!(new_grid::Grid, old_grid::Grid, gap::Int64)

    # Check if the two grid are overlapping 

    if max(old_grid.box.lbnd, new_grid.box.lbnd) <= min(old_grid.box.ubnd, new_grid.box.ubnd)

        # Identifying the overlapping range of the old grid on the new grid 

        parent_range::UnitRange{Int64} = max(old_grid.box.lbnd, new_grid.box.lbnd) : min(old_grid.box.ubnd, new_grid.box.ubnd)

        old_range::UnitRange{Int64} = old_grid.box.inside_ind[ 
            2 * (parent_range[1] - old_grid.box.lbnd) + 1 : 2 * (parent_range[end] - old_grid.box.lbnd) + 2 
        ]

        # Identifying the list of child of the old grid overlapping on the new grid 
        new_children::Vector{Grid} = []

        for old_child in old_grid.children

            # Check if one of the boun is in the overlapping area 
            if (old_child.box.lbnd in old_range) || (old_child.box.ubnd in old_range)
                
                # Updating the boundaries of the old_child
                new_lbnd = old_child.box.lbnd - gap
                new_ubnd = old_child.box.ubnd - gap

                old_child.box.lbnd = new_lbnd
                old_child.box.ubnd = new_ubnd

                # Update of the parent grid 
                old_child.parent = new_grid

                # Update the child list
                push!(new_children, old_child)

            end
        end
        
        return new_children

    else # The two grid are not overlapping
        return nothing
    end

end

function set_ghost_zones!(parent::Grid, child::Grid)
    
    # Declaration of the local variables 
    old_child::Union{Grid, Nothing} = nothing 
    slope::Float64 = 0.0
    gap::Int64 = 0

    i_start::Int64 = -1
    i_end::Int64   = -1

    new_ind::Int64 = -1

    # Left ghost zone

    # Check the boundary type
    if child.box.bbox[1] # Physical boundary
        
        # Apply the left boundary condition 
        for ind in child.box.lgzone_ind
            child.u[ind] = left_bound_cond(0.0)

        end

    else # Grid interface 

        # Declaration of the local variables
        i_start = child.box.lbnd - div(child.N_gzones + 1, 2)
        i_end   = child.box.lbnd - 1

        # Check if the number of ghost zone is odd
        if child.N_gzones % 2 == 1
            
            # Check if the parent control volume was refined 
            old_child = is_refined(parent, i_start)

            if isnothing(old_child) # The parent control volume was not refined 

                # Compute slope
                slope = compute_slope(parent, i_start)
                
                # Update the ghost zone value
                child.u[1] = parent.u[i_start] + 0.25 * slope

            else # The parent grid was refined 
                
                #  Compute the gap 
                gap = 2 * (child.box.lbnd - old_child.box.lbnd)

                # Update the ghost zone value 
                child.u[1] = old_child.u[1 - gap]
                
            end

            # Update the parent index 
            i_start += 1
            
        end

        # loop on the rest of the left ghost zone 
        for parent_ind in i_start : i_end

            #  Identifying the child first control volume index 
            new_ind = child.N_gzones + 1 - (child.box.lbnd - parent_ind) * 2

            #  Check if the pqrent control volume was refined 
            old_child = is_refined(parent, parent_ind)

            if isnothing(old_child) # The parent control volume was not refined

                # Compute the slope 
                slope = compute_slope(parent, parent_ind)

                # Update the ghost zone values
                child.u[new_ind]     = parent.u[parent_ind] - 0.25 * slope
                child.u[new_ind + 1] = parent.u[parent_ind] + 0.25 * slope

            else # The parent control volume was not refined 

                # Compute the gap 
                gap = 2 * (child.box.lbnd - old_child.box.lbnd)

                # Update the ghost zone values 
                child.u[new_ind]     = old_child.u[new_ind + gap]
                child.u[new_ind + 1] = old_child.u[new_ind + 1 + gap]
            end

        end 
        
        

    end

    # Right ghost zone 

    # Check the boundary type 
    if child.box.bbox[2] # Physical boundary 

        # Apply the right boundary condition 
        for ind in child.box.rgzone_ind
            child.u[ind] = left_bound_cond(0.0)
        end

    else # Grid interface

        # Declaration of the local variables 
        i_start = child.box.ubnd + 1
        i_end   = child.box.ubnd + div(child.N_gzones + 1, 2)

        #  Check if the number of ghost zone is odd 
        if child.N_gzones % 2 == 1 # Update the right most ghost zone 

            # Check if the parent control volume was refined 
            old_child = is_refined(parent, i_end)

            if isnothing(old_child) # The parent control volume was not refined 

                #  Compute the slope 
                slope = compute_slope(parent, i_end)

                # Update the ghost zone value 
                child.u[end] = parent.u[i_end] - 0.25 * slope

            else # The parent control volume was refined 

                # Identifying  the child control volume index 
                new_ind = child.box.rgzone_ind[end]

                # Compute the gap 
                gap = 2 * (child.box.lbnd - old_child.box.lbnd)

                # Update the ghost zone value 
                child.u[new_ind] = old_child.u[new_ind + gap]

            end

            # Update the parent index 
            i_end -= 1
        end

        # Loop on the rest of the right ghost zone 
        for parent_ind in i_start : i_end
            
            # Identifying the child first controle volume index 
            new_ind = child.box.rgzone_ind[1] + 2 * (parent_ind - child.box.ubnd - 1) 

            # Check if the parent control volume was refined 
            old_child = is_refined(parent, parent_ind)

            if isnothing(old_child) # The parent control volume was not refined 

                # Compute the slope 
                slope = compute_slope(parent, parent_ind)

                # Update the ghost zone values 
                child.u[new_ind]     = parent.u[parent_ind] - 0.25 * slope
                child.u[new_ind + 1] = parent.u[parent_ind] + 0.25 * slope

            else # The parent control volume was refined 

                # Compute the gap 
                gap = 2 * (child.box.lbnd - old_child.box.lbnd)

                # Update the ghost zone values
                child.u[new_ind] = old_child.u[new_ind + gap]
                child.u[new_ind + 1] = old_child.u[new_ind + 1 + gap]
            end
        end
    end



    
end