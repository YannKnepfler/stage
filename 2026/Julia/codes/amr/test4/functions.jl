


# Function that return the exact solution of the problem
function exact_solution(t::Float64, x::Float64)
    return initial_condition(x - t)
end

function initial_condition(x::Float64)
    return exp(-((x-0.25)/0.05)^2)
end

function compute_flux(u::Vector{Float64}, ind::Int64)
    return u[ind] - u[ind - 1]
end

