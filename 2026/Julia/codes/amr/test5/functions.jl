


# Function that return the exact solution of the problem
function exact_solution(t::Float64, x::Float64)
    return initial_condition(x - t)
end

function left_boundary_condition(t::Float64)
    return exact_solution(t, 0.0)
end

function right_boundary_condition(t::Float64)
    return exact_solution(t, 1.0)
end

function initial_condition(x::Float64)
    return 0.9 * exp(- ((x - 0.25)/0.05)^2) + 0.1
end

function compute_flux(u::Vector{Float64}, ind::Int64)
    return u[ind - 1]
end

