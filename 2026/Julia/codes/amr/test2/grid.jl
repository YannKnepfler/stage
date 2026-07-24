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

# Function that return the boundaries beetween wich to refine 
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