include("box.jl")
include("functions.jl")

# Here will be defined the grid structur and the function linked to this structur


mutable struct grid
    
    # Number of points and ghostzones
    N_points::Int64
    N_gzones::Int64

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
    gbox::box

    # List of the child grid
    children::Vector{grid}
end





# Function that set up the base grid
function setup_base_grid(N_points::Int64, N_gzones::Int64, x0::Float64, dt::Float64, dx::Float64)

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

    gbox::box = box(lbnd, ubnd, N_points, N_gzones, bbox)

    # For now the base grid does not have child grid 
    children::Vector{grid} = Vector{grid}(undef, 0)


    base::grid = grid(N_points, N_gzones, iteration, iteration_step, x0, dt, dx, u, u_p, local_error, x, gbox, children)
    return base
end 

# Function that initialize the values of the grid 
function initialize_grid(base::grid)

    for i in 1:(base.N_points + 2*base.N_gzones)
        base.u[i]           = init_cond(base.x[i])
        base.u_p[i]         = init_cond(base.x[i])
        base.local_error[i] = 0.0
    end
end

# Function that set up a grid nested in a parent grid 
function setup_grid(parent::grid, lbnd::Int64, ubnd::Int64)

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
        if parent.gbox.bbox[1] # Physical boundary condition 
            push!(bbox, true) 

        else # Grids interfaces -> does not ensure proper nesting 
            println("The proper nesting of a grid was not respected")
            exit()
        end

    else
        push!(bbox, false)
    end

    # Right boundary
    if (ubnd == (parent.N_points + parent.N_gzones))
        if parent.gbox.bbox[2] # Physical boundary condition 
            push!(bbox, true) 

        else # Grids interfaces -> does not ensure proper nesting 
            println("The proper nesting of a grid was not respected")
            exit()
        end

    else
        push!(bbox, false)
    end

    gbox::box = box(lbnd, ubnd, N_points, N_gzones, bbox)

    # For now the grid has no child grid 
    children::Vector{grid} = Vector{grid}(undef, 0)

    child::grid = grid(N_points, N_gzones, iteration, iteration_step, x0, dt, dx, u, u_p, local_error, x, gbox, children)
end

# Function that sets the value ofa child grid using a linear prolongation of the parent grid values 
function prolongate_grid!(parent::grid, child::grid)

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

# Function that return the slope used to prolongate grid with minmod
function slope(parent::grid, ind::Int64)
    
    slope::Float64   = 0.0

    if ind == 1
        slope = parent.u[ind + 1] - parent.u[ind]

    elseif ind == parent.N_points + 2 * parent.N_gzones
        slope = parent.u[ind] - parent.u[ind - 1]

    else
        slope_l::Float64 = parent.u[ind] - parent.u[ind - 1]
        slope_r::Float64 = parent.u[ind + 1] - parent.u[ind]

        # minmod
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
