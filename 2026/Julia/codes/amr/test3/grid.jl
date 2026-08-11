# Here the structur Grid and the function related to this structur are defined 

include("box.jl")
include("functions.jl")

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
        local_error::Vector{Float64} = Vector{Float64}()

        # Definition of the position list 
        x::Vector{Float64} = Vector{Float64}(undef, n_points + 2*n_gzones)

        for i in 1 : n_points + 2*n_gzones
            ind = Float64(i - n_gzones)
            x[i] = x0 + (ind - 0.5) * dx
        end

        # Definition of the Box
        box::Box = Box(n_points, n_gzones)

        # The base grid does not have parent 
        parent::Union{Grid, Nothing} = nothing 

        # For nox the base grid does not have child grid
        children::Vector{Grid} = Vector{Grid}(undef, 0)

        new(n_points, n_gzones, x0, dt, dx, u, u_p, local_error, x, box, parent, children)
    end

    # Child grid Constructor
    function Grid(parent::Grid, lbnd::Int64, ubnd::Int64)

        # Definition of the child box 
        box::Box = Box(parent.box, lbnd, ubnd)

        # Number of points and ghost zones 
        n_points::Int64 = box.n_points 
        n_gzones::Int64 = box.n_gzones

        # Definition of the left border position 
        x0::Float64 = parent.x[box.lbnd] - 0.5 * parent.dx

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

        # For now the child grid does not have child grid itself
        children::Vector{Grid} = Vector{Grid}(undef, 0)

        new(n_points, n_gzones, x0, dt, dx, u, u_p, local_error, x, box, parent, children)
    end 
end

function initialize_grid!(base::Grid)

    # Initialization of the solution list and the previous solution list
    for ind in 1 : base.n_points + 2 * base.n_gzones
        base.u[ind]   = exact_solution(0.0, base.x[ind])
        base.u_p[ind] = exact_solution(-base.dt, base.x[ind])
    end
end