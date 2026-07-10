# Ici on définit le type box, servant a l'indexation des maillages


mutable struct box 
    # Position du maillage dans le maillage parent
    lbnd::Int64
    ubnd::Int64

    # Number of points and ghost zones
    N_points::Int64
    N_gzones::Int64

    # Définition des bords du maillages, true pour un bord physique, false pour une interface entre les maillages
    bbox::Vector{Bool}
end



