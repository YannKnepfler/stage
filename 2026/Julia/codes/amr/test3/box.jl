# Here the structure Box and the function related to this structur are defined 



# The structur Box contains information about the placement of the grid on the parent grid, and useful range of index 
mutable struct Box
    
    # Number of points and ghost zones of the grid
    n_points::Int64
    n_gzones::Int64

    # Position of the grid on the parent grid 
    lbnd::Int64 
    ubnd::Int64

    # Information on the boundary condition. true for a physical boundary and false for a grid interface
    bbox::Tuple{Bool, Bool}

    # Index of the solution points 
    inside_ind::UnitRange{Int64}
    lgzone_ind::UnitRange{Int64}
    rgzone_ind::UnitRange{Int64}


    # Constructor

    # Default constructor
    function Box()
        new()
    end

    # Base box constructor
    function Box(n_points::Int64, n_gzones::Int64)
        # The base grid does not have a parent grid 
        lbnd::Int64 = -1
        ubnd::Int64 = -1

        # The base grid has only physical boundary
        bbox::Tuple{Bool, Bool} = (true, true)

        # Definition of the list of index 
        inside_ind::UnitRange{Int64} = (n_gzones + 1) : (n_gzones + n_points)
        lgzone_ind::UnitRange{Int64} = 1 : n_gzones
        rgzone_ind::UnitRange{Int64} = (n_gzones + n_points) : (2 * n_gzones + n_points)

        new(n_points, n_gzones, lbnd, ubnd, bbox, inside_ind, lgzone_ind, rgzone_ind)
    end

    # Child box constructor 
    function Box(parent::Box, lbnd::Int64, ubnd::Int64)

        # Number of ghost zones 
        n_gzones::Int64 = parent.n_gzones

        # Defintion of the boundaries 
        bbox::Tuple{Bool, Bool} = (false, false)

        # Left boundary
        if lbnd == parent.inside_ind[1]
            if parent.bbox[1]
                bbox = (true, bbox[2])
            end
        end

        # Right boundary
        if ubnd == parent.inside_ind[end]
            if parent.bbox[2]
                bbox = (bbox[1], true)
            end
        end

        # Number of points 
        n_points::Int64 = (ubnd - lbnd + 1) * 2 
               

        # Definition of the index range 
        inside_ind::UnitRange{Int64} = (n_gzones + 1) : (n_gzones + n_points)
        lgzone_ind::UnitRange{Int64} = 1 : n_gzones
        rgzone_ind::UnitRange{Int64} = (n_gzones + n_points) : (2 * n_gzones + n_points)


        new(n_points, n_gzones, lbnd, ubnd, bbox, inside_ind, lgzone_ind, rgzone_ind)
    end
end

