# Ici on définit le type grid, servant à décrire un maiilage 


@kwdef struct grid
    
    # Number of points and ghostzones
    N_points::Int8
    N_gzones::Int8

    # iteration is the iteration this f=grid is on 
    # iteration_step isthe number of point need to evolve for one time step

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
end