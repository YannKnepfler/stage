# Ici on définit le type box, servant a l'indexation des maillages


mutable struct box 
    # Position of the grid inside the grid parent
    lbnd::Int64
    ubnd::Int64

    # Number of points and ghost zones
    N_points::Int64
    N_gzones::Int64

    # Definition of the boundary condition of the mesh: true for a physical boundary and false for un grid interface boundary
    bbox::Vector{Bool}
end



