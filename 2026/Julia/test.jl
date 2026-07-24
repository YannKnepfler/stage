boundaries::Vector{Tuple{Int64, Int64}} = [(2, 11)]


println(boundaries)


a = 21

c, _ = pop!(boundaries)

push!(boundaries, (c, a))

println(boundaries)