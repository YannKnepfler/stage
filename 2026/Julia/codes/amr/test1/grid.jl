include("box.jl")
include("functions.jl")

# Here will be defined the grid structur and the function linked to this structur


mutable struct Grid
    
    # Number of points and ghostzones
    N_points::Int64
    N_gzones::Int64
    N_buffer::Int64

    # iteration is the iteration this grid is on 
    # iteration_step is the number of point need to evolve for one time step

    iteration::Int64
    iteration_step::Int64

    # Location of the left boundarie of the grid
    x0::Float64

    # Time step and cell size
    dt::Float64
    dx::Float64

    # List of the solution, previous solution and local_error
    u::Vector{Float64}
    u_p::Vector{Float64}
    local_error::Vector{Float64}
    
    # List of the position of the centers of the cells
    x::Vector{Float64}

    # Box of the grid
    gbox::Box

    # List of the child grid
    children::Vector{Grid}

    # Constructor of the base grid 
    function Grid(N_points::Int64, N_gzones::Int64, x0::Float64, dt::Float64, dx::Float64) 
        # Defintion of the base grid
        iteration::Int64              = 0
        iteration_step::Int64         = 1

        # Definition of the vectors
        u::Vector{Float64}           = Vector{Float64}(undef, N_points + 2*N_gzones)
        u_p::Vector{Float64}         = Vector{Float64}(undef, N_points + 2*N_gzones)
        local_error::Vector{Float64} = Vector{Float64}(undef, N_points + 2*N_gzones)
        x::Vector{Float64}           = Vector{Float64}(undef, N_points + 2*N_gzones)

        for i in 1:(N_points + 2*N_gzones)
            ind = Float64(i - N_gzones)
            x[i] = x0 + (ind - 0.5) * dx
        end

        # Definition of the box
        # The base grid does not have parent grid 
        lbnd::Int64 = -1
        ubnd::Int64 = -1

        # The base grid has only physical boundaries
        bbox::Vector{Bool} = [true, true]

        gbox::Box = Box(lbnd, ubnd, N_points, N_gzones, bbox)

        # For now the base grid does not have child grid 
        children::Vector{Grid} = Vector{Grid}(undef, 0)

        # Inner constructor 
        new(N_points, N_gzones, iteration, iteration_step, x0, dt, dx, u, u_p, local_error, x, gbox, children)
    end

    # Constructor of a grid using a parent grid
    function Grid(parent::Grid, lbnd::Int64, ubnd::Int64)
        # Definition of the child grid
        N_points::Int64              = (ubnd - lbnd + 1) * 2
        N_gzones::Int64              = parent.N_gzones

        iteration::Int64             = parent.iteration
        iteration_step::Int64        = parent.iteration_step

        # Definition of the left boundary position of the child grid
        x0::Float64                  = parent.x[lbnd] - 0.5 * parent.dx

        dt::Float64                  = parent.dt
        dx::Float64                  = parent.dx * 0.5

        # Definition of the vectors
        u::Vector{Float64}           = Vector{Float64}(undef, N_points + 2*N_gzones)
        u_p::Vector{Float64}         = Vector{Float64}(undef, N_points + 2*N_gzones)
        local_error::Vector{Float64} = Vector{Float64}(undef, N_points + 2*N_gzones)

        # Definition of the position of the center of the cells 
        x::Vector{Float64}           = Vector{Float64}(undef, N_points + 2*N_gzones)

        for i in 1:(N_points + 2*N_gzones)
            ind = Float64(i - N_gzones)
            x[i] = x0 + (ind - 0.5) * dx
        end


        # Definition of the box
        # Definition of the boundary type
        bbox::Vector{Bool} = []

        # Left Boundary
        if (lbnd == (parent.N_gzones + 1))
            # if parent.gbox.bbox[1] # Physical boundary condition 
            #     push!(bbox, true) 

            # else # Grids interfaces -> does not ensure proper nesting 
            #     error("The proper nesting of a grid was not respected")
            #     exit()
            # end

            push!(bbox, parent.gbox.bbox[1])

        else
            push!(bbox, false)
        end

        # Right boundary
        if (ubnd == (parent.N_points + parent.N_gzones))
            # if parent.gbox.bbox[2] # Physical boundary condition 
            #     push!(bbox, true) 

            # else # Grids interfaces -> does not ensure proper nesting 
            #     println("The proper nesting of a grid was not respected")
            #     exit()
            # end

            push!(bbox, parent.gbox.bbox[2])

        else
            push!(bbox, false)
        end

        gbox::Box = Box(lbnd, ubnd, N_points, N_gzones, bbox)

        # For now the grid has no child grid 
        children::Vector{Grid} = Vector{Grid}(undef, 0)

        new(N_points, N_gzones, iteration, iteration_step, x0, dt, dx, u, u_p, local_error, x, gbox, children)

    end

end







# Function that initialize the values of the grid 
function initialize_grid!(base::Grid)

    for i in 1:(base.N_points + 2*base.N_gzones)
        base.u[i]           = init_cond(base.x[i])
        base.u_p[i]         = init_cond(base.x[i])
        base.local_error[i] = 0.0
    end
end



# Function that sets the value of a child grid using a linear prolongation of the parent grid values 
function prolongate_grid!(parent::Grid, child::Grid)

    # Linear prolongation using minmod
    i_start = 1
    i_end = child.N_points + 2 * child.N_gzones

    if child.N_gzones % 2 == 1
        # Left most ghost zone
        parent_i = child.gbox.lbnd - div(child.N_gzones + 1 , 2) 
        

        child.u[1] = parent.u[parent_i] + 0.25 * slope(parent, parent_i)

        # Right most ghost zone
        parent_i = child.gbox.ubnd + div(child.N_gzones + 1 , 2)

        child.u[end] = parent.u[parent_i] - 0.25 * slope(parent, parent_i)

        i_start += 1
        i_end   -= 1
    end

    for i in i_start:2:i_end

        parent_i = child.gbox.lbnd + div(i - child.N_gzones - 1, 2)

        sl = slope(parent, parent_i)

        child.u[i]   = parent.u[parent_i] - 0.25 * sl
        child.u[i+1] = parent.u[parent_i] + 0.25 * sl
        
    end
end

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
    
    return gu_m, gu_p
end

# Return a list of new child of the grid
function refine_grid(grid::Grid, tol::Float64, N_points_base::Int64)
    
    # Definition of the local variables 
    children::Vector{Grid} = Vector{Grid}(undef, 0)
    boundaries::Vector{Tuple{Int64, Int64}} = []

    N_patches = div(grid.N_points, N_points_base)

    # Loop on the patches of the grid 
    for i in 1:N_patches

        # Definition of the varaibles of the loop
        lbnd::Int64 = -1
        ubnd::Int64 = -1

        # Breakdown of the grid
        index::UnitRange{Int64} = grid.N_gzones + (i-1) * N_points_base + 1 : grid.N_gzones + i * N_points_base
        buff_u::Vector{Float64} = grid.u[index]

        # Find the boundarie of the refinement 
        lbnd, ubnd = refine_bnd(buff_u, index, tol)

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


    if length(boundaries) != 0
        for bnd in boundaries
            lbnd, ubnd = bnd 

            push!(children, Grid(grid, lbnd, ubnd)) 
        end
    else
        return nothing
    end

    return children
end

# Return lbnd and ubnd 
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


