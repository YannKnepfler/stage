# Ici on définit le type box, servant a l'indexation des maillages


@kwdef struct box 
    # Position du maillage dans le maillage parent
    lbnd::Int8
    ubnd::Int8

    # Number of points and ghost zones
    N_points::Int8
    N_gzones::Int8

    # Définition des bords du maillages, true pour un bord physique, false pour une interface entre les maillages
    bbox::Vector{Bool}
end

