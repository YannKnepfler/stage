include("box.jl")
include("functions.jl")

# Ici on définit le type grid, servant à décrire un maiilage 


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
    childs::Vector{grid}
end





# Function that set up the base grid
function setup_base_grid(N_points::Int64, N_gzones::Int64, x0::Float64, dt::Float64, dx::Float64)

    # Defintion of the base_grid
    iteration::Int64              = 0
    iteration_step::Int64         = 1

    dt::Float64                  = 0.1
    dx::Float64                  = 0.1

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
    childs::Vector{grid} = Vector{grid}(undef, 0)


    base::grid = grid(N_points, N_gzones, iteration, iteration_step, x0, dt, dx, u, u_p, local_error, x, gbox, childs)
    return base
end 

function initialize_grid(base::grid)

    for i in 1:(base.N_points + 2*base.N_gzones)
        base.u[i]           = init_cond(base.x[i])
        base.u_p[i]         = init_cond(base.x[i])
        base.local_error[i] = 0.0
    end
end