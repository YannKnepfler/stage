
function left_bound_cond(t::Float64)
    return t
end

function right_bound_cond(t::Float64)
    return t
end
function init_cond(x::Float64)

    return exp(-1000.0 * (x - 0.5)^2)

end

