using Plots

# ============================================================
# One-dimensional linear advection equation:
#
#     ∂u/∂t + a ∂u/∂x = 0
#
# Numerical method:
# - First-order upwind scheme in space
# - Explicit Euler method in time
# - Periodic boundary conditions
# ============================================================

# Spatial domain
x_min = 0.0
x_max = 1.0

# Advection velocity
a = 1.0

# Number of grid points
number_of_points = 64

# Spatial step
dx = (x_max - x_min) / number_of_points

# Grid coordinates
x = collect(range(x_min, x_max - dx, length = number_of_points))

# CFL number
CFL = 0.9

# Time step
dt = CFL * dx / a

# Final simulation time
final_time = 0.5

# Initial Gaussian wave with amplitude equal to 1
wave_center = 0.25
wave_width = 0.05

u = exp.(-((x .- wave_center) ./ wave_width).^2)

# Store the initial solution
u_initial = copy(u)

# Initialize time
time = 0.0

# Time integration loop
while time < final_time
    global time

    # Adjust the last time step
    current_dt = min(dt, final_time - time)

    # Store the previous solution
    u_old = copy(u)

    # Apply the upwind scheme to the internal grid points
    for i in 2:number_of_points
        u[i] = u_old[i] -
               a * current_dt / dx *
               (u_old[i] - u_old[i - 1])
    end

    # Periodic boundary condition at the left boundary
    u[1] = u_old[1] -
           a * current_dt / dx *
           (u_old[1] - u_old[end])

    # Advance time
    global time += current_dt
end

# Plot the initial and final solutions
plot(
    x,
    u_initial,
    label = "Initial wave",
    linewidth = 2,
    linestyle = :dash,
    xlabel = "x",
    ylabel = "u(x,t)",
    title = "Linear Advection with Explicit Upwind Scheme",
    xlims = (0.0, 1.0),
    ylims = (0.0, 1.1),
    legend = :topright
)

plot!(
    x,
    u,
    label = "Wave at t = $(final_time)",
    linewidth = 3
)

display(current())