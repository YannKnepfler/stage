
function bound_cond(t::Float64)
    return t
end
function init_cond(x::Float64)

    return 1000.0 * exp(  (x - 0.5)^2 )
    
end

