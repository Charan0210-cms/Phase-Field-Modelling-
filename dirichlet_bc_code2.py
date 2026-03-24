"""
8.2  2D Heat Equation – FDM and Dirichlet BC (Code 2)
Interactive parameter entry.  Saves contour plots at selected time steps
and assembles them into a GIF.
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

ic     = int(input("Enter initial condition temperature [0]: ")   or 0)
top    = int(input("Enter temp value at top    [200]: ") or 200)
bottom = int(input("Enter temp value at bottom [0]: ")   or 0)
left   = int(input("Enter temp value at left   [0]: ")   or 0)
right  = int(input("Enter temp value at right  [0]: ")   or 0)

#  Grid 
x_vec = np.linspace(0, w, Nx)
y_vec = np.linspace(0, h, Ny)

dx = x_vec[1] - x_vec[0]
dy = y_vec[1] - y_vec[0]
x2, y2 = dx**2, dy**2

dt = (x2 * y2) / (2 * D * (x2 + y2))   # maximum stable time step
kx = D * dt / x2
ky = D * dt / y2

print(f"kx = {kx:.6f},  ky = {ky:.6f}  (both must be ≤ 0.5 for stability)")

T = np.zeros((Ny, Nx, Nt))

# Initial condition
T[:, :, 0] = ic

# Boundary conditions
T[Nx - 1, :, :] = top
T[0,      :, :] = bottom
T[:,  0,  :] = left
T[:, Ny-1, :] = right

#  Time-stepping 
for k in range(1, Nt):
    for i in range(1, Nx - 1):
        for j in range(1, Ny - 1):
            H = kx * (T[j, i + 1, k - 1] + T[j, i - 1, k - 1] - 2 * T[j, i, k - 1])
            V = ky * (T[j + 1, i, k - 1] + T[j - 1, i, k - 1] - 2 * T[j, i, k - 1])
            T[j, i, k] = T[j, i, k - 1] + H + V

#  Save contour plots at every time step & build GIF 
os.makedirs("2D_DN_BC_frames", exist_ok=True)
image_frames = []

for k in range(0, Nt):
    plt.contourf(x_vec, y_vec, T[:, :, k])
    plt.colorbar()
    plt.title(f'Temperature at timestep {k + 1}')
    fname = f"2D_DN_BC_frames/frame_{k + 1:04d}.jpg"
    plt.savefig(fname)
    plt.clf()
    image_frames.append(Image.open(fname))

image_frames[0].save(
    "8_2_heat_equation_dirichlet_evolution.gif",
    format='GIF',
    append_images=image_frames[1:],
    save_all=True,
    duration=t,
    loop=0
)
print("Saved: 8_2_heat_equation_dirichlet_evolution.gif")

# Quick 4-panel plot 
fig, axes = plt.subplots(2, 2, figsize=(10, 8))
steps = [0, Nt // 3, 2 * Nt // 3, Nt - 1]
titles = [f'timestep {s + 1}' for s in steps]

for ax, s, ttl in zip(axes.flat, steps, titles):
    cf = ax.contourf(x_vec, y_vec, T[:, :, s])
    fig.colorbar(cf, ax=ax)
    ax.set_title(f'Temperature at {ttl}')

plt.tight_layout()
plt.savefig("8_2_heat_equation_dirichlet_panels.png", dpi=150)
plt.show()
print("Saved: 8_2_heat_equation_dirichlet_panels.png")
