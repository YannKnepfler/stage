# here are defined the parameters of the problem 

# AMR parameters
const global refinement_tolerance::Float64 = 0.3
const global coarsening_tolerance::Float64 = 0.1
const global amr_buffer_width::Int64 = 4

global max_depth = 3

const global epsilon = 0.01
    
# 
const global cfl = 1.0