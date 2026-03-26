"""
10.1  Numerical Modelling of Pearlite Dissolution – Dirichlet Boundary Condition

Governing PDE: Fick's Law of Diffusion
    dC/dt = D * (d²C/dx² + d²C/dy²)

Domain represents two pearlite lamellae along the x-axis:
    - First  lamella : 0–44 % of Nx → Ferrite  (C = 0.02 wt%)
                       44–50% of Nx → Cementite (C = 6.67 wt%)
    - Second lamella : 50–94% of Nx → Ferrite  (C = 0.02 wt%)
                       94–100% of Nx→ Cementite (C = 6.67 wt%)

Dirichlet BCs are applied on all four edges.
"""

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import os


w  = int(input("Enter length of the material (x-axis) [1]: "))   
h  = int(input("Enter width  of the material (y-axis) [1]: "))  
t  = int(input("Enter total time interval (ms) [500]: "))         

Nx = int(input("Enter number of grid points for length [100]: ")) 
Ny = int(input("Enter number of grid points for width  [100]: ")) 
Nt = int(input("Enter number of grid points for time   [500]: ")) 

D1 = int(input("Enter Diffusion coefficient value (×10⁻⁵) [1]: ")) 
D  = D1 / 100000

top    = int(input("Enter Carbon conc. value at top    [0]: ")) 
bottom = int(input("Enter Carbon conc. value at bottom [0]: ")) 
left   = int(input("Enter Carbon conc. value at left   [0]: ")) 
right  = int(input("Enter Carbon conc. value at right  [0]: ")) 

#  Grid 
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

# Initial condition: two pearlite lamellae 
# Lamella 1: Ferrite  0 – 44% | Cementite 44 – 50%
# Lamella 2: Ferrite 50 – 94% | Cementite 94 – 100%
FERRITE_C   = 0.02   # wt% carbon in ferrite
CEMENTITE_C = 6.67   # wt% carbon in cementite

for c in range(0, Nx):
    if (c >= 0 and c < int(0.44 * Nx)) or (c >= int(0.50 * Nx) and c < int(0.94 * Nx)):
        T[:, c, 0] = FERRITE_C
    else:
        T[:, c, 0] = CEMENTITE_C

# Dirichlet BCs (all edges)
T[Nx - 1, :, :] = top
T[0,      :, :] = bottom
T[:, 0,   :] = left
T[:, Ny-1, :] = right

# Time-stepping 
for k in range(1, Nt):
    for i in range(1, Nx - 1):
        for j in range(1, Ny - 1):
            H = kx * (T[j, i + 1, k - 1] + T[j, i - 1, k - 1] - 2 * T[j, i, k - 1])
            V = ky * (T[j + 1, i, k - 1] + T[j - 1, i, k - 1] - 2 * T[j, i, k - 1])
            T[j, i, k] = T[j, i, k - 1] + H + V

#  4-panel snapshot 
fig, axes = plt.subplots(2, 2, figsize=(10, 8))
steps  = [0, Nt // 3, 2 * Nt // 3, Nt - 1]

for ax, s in zip(axes.flat, steps):
    cf = ax.contourf(x_vec, y_vec, T[:, :, s])
    fig.colorbar(cf, ax=ax)
    ax.set_title(f'Concentration profile at timestep {s + 1}')

plt.suptitle('Pearlite Dissolution – Dirichlet BC\nCarbon concentration evolution')
plt.tight_layout()
plt.savefig("10_1_pearlite_dirichlet_panels.png", dpi=150)
plt.show()
print("Saved: 10_1_pearlite_dirichlet_panels.png")

#  Save GIF 
os.makedirs("pearlite_dirichlet_frames", exist_ok=True)
image_frames = []

for k in range(0, Nt):
    plt.contourf(x_vec, y_vec, T[:, :, k])
    plt.colorbar()
    plt.title(f'Concentration profile at timestep {k + 1}')
    fname = f"pearlite_dirichlet_frames/frame_{k + 1:04d}.jpg"
    plt.savefig(fname)
    plt.clf()
    image_frames.append(Image.open(fname))

image_frames[0].save(
    "10_1_pearlite_dirichlet_evolution.gif",
    format='GIF',
    append_images=image_frames[1:],
    save_all=True,
    duration=t,
    loop=0
)
print("Saved: 10_1_pearlite_dirichlet_evolution.gif")
