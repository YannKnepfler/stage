# Here the structure Box and the function related to this structur are defined 



# The structur Box contains information about the placement of the grid on the parent grid, and useful range of index 
mutable struct Box
    
    # Number of points and ghost zones of the grid
    N_points::Int64
    N_gzones::Int64

    # Position of the grid on the parent grid 
    lbnd::Int64 
    ubnd::Int64

    # Information on the boundary condition. true for a physical boundary and false for a grid interface
    bbox::Tuple{Bool, Bool}

    # Index of the solution points 
    inside_ind::UnitRange{Int64}
    lgzone_ind::UnitRange{Int64}
    rgzone_ind::UnitRange{Int64}

    # The base box does not need buffer zones 
    without_bz::UnitRange{Int64}


    # Constructor

    # Default constructor
    function Box()
        new()
    end

    # Base box constructor
    function Box(N_points::Int64, N_gzones::Int64)
        # The base grid does not have a parent grid 
        lbnd::Int64 = -1
        ubnd::Int64 = -1

        # The base grid has only physical boundary
        bbox::Tuple{Bool, Bool} = (true, true)

        # Definition of the list of index 
        inside_ind::UnitRange{Int64} = (N_gzones + 1) : (N_gzones + N_points)
        lgzone_ind::UnitRange{Int64} = 1 : N_gzones
        rgzone_ind::UnitRange{Int64} = (N_gzones + N_points) : (2 * N_gzones + N_points)
        without_bz::UnitRange{Int64} = inside_ind

        new(N_points, N_gzones, lbnd, ubnd, bbox, inside_ind, lgzone_ind, rgzone_ind, without_bz)
    end

    # Child box constructor 
    function Box(parent::Box, lbnd::Int64, ubnd::Int64)

        # Number of ghost zones 
        N_gzones::Int64 = parent.N_gzones

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

        # Creation of buffer zone if needed and definition of the inside and non buffer cells index range
        new_lbnd::Int64 = lbnd
        new_ubnd::Int64 = ubnd

        # Left border
        i_start = N_gzones + 1
        if !(bbox[1])
            new_lbnd = lbnd - 1
            i_start += 2
        end

        # Right border
        i_end = i_start - 1 + (ubnd - lbnd + 1) * 2
        if !(bbox[2])
            new_ubnd = ubnd + 1
        end

        without_bz::UnitRange{Int64} = i_start : i_end

        # Number of points 
        N_points::Int64 = (new_ubnd - new_lbnd + 1) * 2
               

        # Definition of the index range 
        inside_ind::UnitRange{Int64} = (N_gzones + 1) : (N_gzones + N_points)
        lgzone_ind::UnitRange{Int64} = 1 : N_gzones
        rgzone_ind::UnitRange{Int64} = (N_gzones + N_points) : (2 * N_gzones + N_points)

 
        

        new(N_points, N_gzones, new_lbnd, new_ubnd, bbox, inside_ind, lgzone_ind, rgzone_ind, without_bz)
    end
end

