include("box.jl")
include("functions.jl")

# Ici on définit le type grid, servant à décrire un maiilage 


mutable struct grid
    
    # Number of points and ghostzones
    N_points::Int8
    N_gzones::Int8

    # iteration is the iteration this grid is on 
    # iteration_step is the number of point need to evolve for one time step

    iteration::Int8
    iteration_step::Int8

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
function setup_base_grid(N_points::Int8, N_gzones::Int8, x0::Float64, dt::Float64, dx::Float64)

    # Défintion of the base_grid
    base::grid = grid(N_points, N_gzones, 0, 1, x0, dt, dx, Vector{Float64}(undef, N_points + 2*N_gzones), Vector{Float64}(undef, N_points + 2*N_gzones), Vector{Float64}(undef, N_gzones), Vector{Float64}(undef, N_points + 1), box(-1, -1, N_points, N_gzones, [true, true]), Vector{grid}(undef,1))

    # Definition of the list of the positions of the centers of the control volumes
    base.x = [x0 + 0.5*dx]

    for i in 1:(N_points-1)
        append!(base.x, base.x[i] + dx )
    end

    return base

end 

function initialize_grid(base::grid)
end