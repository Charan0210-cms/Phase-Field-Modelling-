"""
8.1  2D Heat Equation – FDM and Dirichlet BC (Code 1)
Boundary conditions:
    u_top    = 100.0
    u_left   = 20.0
    u_bottom = 10.0
    u_right  = 0.0

Visualises the evolving temperature field as an animated heatmap
and saves it as a GIF.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Parameters 
plate_length  = 50       # grid points in each spatial direction
max_iter_time = 750      # number of time steps
alpha         = 2.0      # thermal diffusivity
delta_x       = 1        # spatial step

# Stability: gamma must be ≤ 0.25 for 2-D explicit scheme
delta_t = (delta_x ** 2) / (4 * alpha)
gamma   = alpha * delta_t / delta_x ** 2
print(f"gamma = {gamma:.4f}  ({'stable' if gamma <= 0.25 else 'UNSTABLE'})")

#  Boundary conditions 
u_top    = 100.0
u_left   = 20.0
u_bottom = 10.0
u_right  = 0.0

#  Initialise solution array u[time, y, x]
u         = np.empty((max_iter_time, plate_length, plate_length))
u_initial = 10              # interior initial condition

u.fill(u_initial)

# Apply boundary conditions (constant for all time steps)
u[:, (plate_length - 1):, :]  = u_top
u[:, :, :1]                   = u_left
u[:, :1, 1:]                  = u_bottom
u[:, :, (plate_length - 1):]  = u_right

# Time-stepping 
def calculate(u):
    for k in range(0, max_iter_time - 1):
        for i in range(1, plate_length - 1, delta_x):
            for j in range(1, plate_length - 1, delta_x):
                u[k + 1, i, j] = (gamma * (u[k][i + 1][j] + u[k][i - 1][j]
                                            + u[k][i][j + 1] + u[k][i][j - 1]
                                            - 4 * u[k][i][j])
                                   + u[k][i][j])
    return u

u = calculate(u)

#  Animation 
def plotheatmap(u_k, k):
    plt.clf()
    plt.title(f"Temperature at t = {k * delta_t:.3f} unit time")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.pcolormesh(u_k, cmap=plt.cm.jet, vmin=0, vmax=100)
    plt.colorbar()
    return plt

def animate(k):
    plotheatmap(u[k], k)

anim = animation.FuncAnimation(
    plt.figure(), animate,
    interval=1, frames=max_iter_time, repeat=False
)
anim.save("8_1_heat_equation_dirichlet.gif")
print("Saved: 8_1_heat_equation_dirichlet.gif")

# Static snapshot at final time step 
plotheatmap(u[-1], max_iter_time - 1)
plt.savefig("8_1_heat_equation_dirichlet_final.png", dpi=150)
plt.show()
print("Saved: 8_1_heat_equation_dirichlet_final.png")
