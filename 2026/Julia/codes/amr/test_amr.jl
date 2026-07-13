using Plots

include("./grid.jl")
# # Programme qui test les fonctionnalité des modules de amr

base::grid = setup_base_grid(10, 0, 0.0, 0.1, 1.0)
initialize_grid(base)

push!(base.children, setup_grid(base, base.N_gzones + 1, base.N_gzones + base.N_points ))


prolongate_grid!(base, base.children[1])

child::grid = base.children[1]


println(base.x)
println(child.x)
# println(child.N_points)
fig = plot([child.x, base.x], [child.u, base.u], label = ["child" "base"])

savefig(fig, "./test.png")


