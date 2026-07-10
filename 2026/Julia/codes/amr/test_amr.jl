using Plots

include("./grid.jl")
# Programme qui test les fonctionnalité des modules de amr

base::grid = setup_base_grid(10, 1, 0.0, 0.1, 0.1)
initialize_grid(base)

# fig = plot(base.x, base.u, label = "test")

# savefig(fig, filename = "./test.png")

