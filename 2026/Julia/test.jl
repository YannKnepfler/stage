# define the Coord type
@kwdef struct Coord 
    x::UInt16
    y::UInt16
end
# define the Plot keyword type
@kwdef struct Plot
    bottom_left::Coord
    top_right::Coord
end

function is_claim_staked(claim::Plot, register::Set{Plot})
        return claim in register
end

function stake_claim!(claim::Plot, register::Set{Plot})
    if is_claim_staked(claim, register)
        return false
    else
        push!(register, claim)
        return true
    end
end

function get_longest_side(claim::Plot)
    x_dist = claim.top_right.x - claim.bottom_left.x
    y_dist = claim.top_right.y - claim.bottom_left.y


    if x_dist > y_dist
        return x_dist
    else 
        return y_dist
    end    
end

function get_claim_with_longest_side(register::Set{Plot})
    res = Plot(Coord(0, 0), Coord(0,0))
    ref_side::UInt16 = 0


    for claim in register
        if ref_side < get_longest_side(claim)
            res = claim
            ref_side = get_longest_side(claim)
        end
    end
    return res
end

plot1 = Plot(bottom_left=Coord(5, 3), top_right=Coord(10, 5))
plot2 = Plot(bottom_left=Coord(0, 0), top_right=Coord(3, 3))

register = Set{Plot}([plot1])

print(get_claim_with_longest_side(register))