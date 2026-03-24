"""
8.3  2D Heat Equation – FDM and Mixed BC (Dirichlet and Neumann)
    u_top    (Dirichlet) = 0.0
    u_bottom (Dirichlet) = 0.0
    u_left   (Neumann)   = 50.0   (flux dT/dx at left)
    u_right  (Neumann)   = 50.0   (flux dT/dx at right)

Ghost-point technique applied on the Neumann (left/right) edges.
"""

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import os

# Interactive inputs 
w  = int(input("Enter length of the material (x-axis) [1]: ")   or 1)
h  = int(input("Enter width  of the material (y-axis) [1]: ")   or 1)
t  = int(input("Enter total time interval (ms) [500]: ")         or 500)

Nx = int(input("Enter number of grid points for length [50]: ")  or 50)
Ny = int(input("Enter number of grid points for width  [50]: ")  or 50)
Nt = int(input("Enter number of grid points for time   [500]: ") or 500)

D1 = int(input("Enter Thermal Diffusivity value (×10⁻⁵) [1]: ") or 1)
D  = D1 / 100000

ic     = int(input("Enter initial condition temperature [0]: ")        or 0)
top    = int(input("Enter Dirichlet temp at top    [0]: ")             or 0)
bottom = int(input("Enter Dirichlet temp at bottom [0]: ")             or 0)
left   = int(input("Enter Neumann dT/dx flux at left  [50]: ")        or 50)
right  = int(input("Enter Neumann dT/dx flux at right [50]: ")        or 50)

# Grid 
x_vec = np.linspace(0, w, Nx)
y_vec = np.linspace(0, h, Ny)

dx = x_vec[1] - x_vec[0]
dy = y_vec[1] - y_vec[0]
x2, y2 = dx**2, dy**2

dt = (x2 * y2) / (2 * D * (x2 + y2))
kx = D * dt / x2
ky = D * dt / y2

print(f"kx = {kx:.6f},  ky = {ky:.6f}")

T = np.zeros((Ny, Nx, Nt))
T[:, :, 0] = ic

# Dirichlet BCs (top & bottom, held for all time)
T[Nx - 1, :, :] = top
T[0,      :, :] = bottom

#  Time-stepping with Neumann on left/right via ghost points
for k in range(1, Nt):
    # Neumann left boundary (j = 0): ghost point u[-1] = u[1] - 2*dx*left
    for j in range(1, Ny - 1):
        T[j, 0, k] = (kx * (2 * T[j, 1, k - 1] - (2 * dx * left) - 2 * T[j, 0, k - 1])
                      + T[j, 0, k - 1]
                      + ky * (T[j + 1, 0, k - 1] + T[j - 1, 0, k - 1] - 2 * T[j, 0, k - 1]))

    # Interior nodes
    for i in range(1, Nx - 1):
        for j in range(1, Ny - 1):
            H = kx * (T[j, i + 1, k - 1] + T[j, i - 1, k - 1] - 2 * T[j, i, k - 1])
            V = ky * (T[j + 1, i, k - 1] + T[j - 1, i, k - 1] - 2 * T[j, i, k - 1])
            T[j, i, k] = T[j, i, k - 1] + H + V

    # Neumann right boundary (j = Ny-1): ghost point u[Ny] = u[Ny-2] + 2*dx*right
    for j in range(1, Ny - 1):
        T[j, Ny - 1, k] = (kx * (2 * T[j, Ny - 2, k - 1] + (2 * dx * right) - 2 * T[j, Ny - 1, k - 1])
                           + T[j, Ny - 1, k - 1]
                           + ky * (T[j + 1, Ny - 1, k - 1] + T[j - 1, Ny - 1, k - 1]
                                   - 2 * T[j, Ny - 1, k - 1]))

#  4-panel snapshot plot 
fig, axes = plt.subplots(2, 2, figsize=(10, 8))
steps = [0, Nt // 3, 2 * Nt // 3, Nt - 1]

for ax, s in zip(axes.flat, steps):
    cf = ax.contourf(x_vec, y_vec, T[:, :, s])
    fig.colorbar(cf, ax=ax)
    ax.set_title(f'Temperature at timestep {s + 1}')

plt.suptitle('2D Heat Equation – Mixed BC (Dirichlet top/bottom, Neumann left/right)')
plt.tight_layout()
plt.savefig("8_3_heat_equation_mixed_bc_panels.png", dpi=150)
plt.show()
print("Saved: 8_3_heat_equation_mixed_bc_panels.png")

# Save GIF 
os.makedirs("2D_mixed_BC_frames", exist_ok=True)
image_frames = []

for k in range(0, Nt):
    plt.contourf(x_vec, y_vec, T[:, :, k])
    plt.colorbar()
    plt.title(f'Temperature at timestep {k + 1}')
    fname = f"2D_mixed_BC_frames/frame_{k + 1:04d}.jpg"
    plt.savefig(fname)
    plt.clf()
    image_frames.append(Image.open(fname))

image_frames[0].save(
    "8_3_heat_equation_mixed_bc_evolution.gif",
    format='GIF',
    append_images=image_frames[1:],
    save_all=True,
    duration=t,
    loop=0
)
print("Saved: 8_3_heat_equation_mixed_bc_evolution.gif")
