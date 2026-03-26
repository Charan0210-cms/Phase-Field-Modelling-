"""
10.2  Numerical Modelling of Pearlite Dissolution – Neumann Boundary Condition

Same two-lamellae setup as section 10.1 but with Neumann (flux) BCs
applied on all four edges via the ghost-point technique.
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

top    = int(input("Enter dC/dx flux at top    [0]: ")) 
bottom = int(input("Enter dC/dx flux at bottom [0]: ")) 
left   = int(input("Enter dC/dx flux at left   [0]: ")) 
right  = int(input("Enter dC/dx flux at right  [0]: ")) 

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

# Initial condition: two pearlite lamellae 
FERRITE_C   = 0.02
CEMENTITE_C = 6.67

for c in range(0, Nx):
    if (c >= 0 and c < int(0.44 * Nx)) or (c >= int(0.50 * Nx) and c < int(0.94 * Nx)):
        T[:, c, 0] = FERRITE_C
    else:
        T[:, c, 0] = CEMENTITE_C

# Time-stepping with Neumann ghost-point BCs on all four edges 
for k in range(1, Nt):
    # Bottom edge  (j = 0):  ghost → T[-1] = T[1] - 2*dy*bottom
    for i in range(1, Nx - 1):
        T[0, i, k] = (ky * (2 * T[1, i, k - 1] - (2 * dx * bottom) - 2 * T[0, i, k - 1])
                      + T[0, i, k - 1]
                      + kx * (T[0, i + 1, k - 1] + T[0, i - 1, k - 1] - 2 * T[0, i, k - 1]))

    # Interior nodes
    for i in range(1, Nx - 1):
        for j in range(1, Ny - 1):
            H = kx * (T[j, i + 1, k - 1] + T[j, i - 1, k - 1] - 2 * T[j, i, k - 1])
            V = ky * (T[j + 1, i, k - 1] + T[j - 1, i, k - 1] - 2 * T[j, i, k - 1])
            T[j, i, k] = T[j, i, k - 1] + H + V

    # Right edge  (i = Ny-1): ghost → T[Ny] = T[Ny-2] + 2*dx*right
    for j in range(1, Ny - 1):
        T[j, Ny - 1, k] = (kx * (2 * T[j, Ny - 2, k - 1] + (2 * dx * right) - 2 * T[j, Ny - 1, k - 1])
                           + T[j, Ny - 1, k - 1]
                           + ky * (T[j + 1, Ny - 1, k - 1] + T[j - 1, Ny - 1, k - 1]
                                   - 2 * T[j, Ny - 1, k - 1]))

    # Top edge  (j = Nx-1): ghost → T[Nx] = T[Nx-2] + 2*dy*top
    for i in range(1, Nx - 1):
        T[Nx - 1, i, k] = (ky * (2 * T[Nx - 2, i, k - 1] + (2 * dx * top) - 2 * T[Nx - 1, i, k - 1])
                           + T[Nx - 1, i, k - 1]
                           + kx * (T[Nx - 1, i + 1, k - 1] + T[Nx - 1, i - 1, k - 1]
                                   - 2 * T[Nx - 1, i, k - 1]))

#  4-panel snapshot 
fig, axes = plt.subplots(2, 2, figsize=(10, 8))
steps = [0, Nt // 3, 2 * Nt // 3, Nt - 1]

for ax, s in zip(axes.flat, steps):
    cf = ax.contourf(x_vec, y_vec, T[:, :, s])
    fig.colorbar(cf, ax=ax)
    ax.set_title(f'Concentration profile at timestep {s + 1}')

plt.suptitle('Pearlite Dissolution – Neumann BC\nCarbon concentration evolution')
plt.tight_layout()
plt.savefig("10_2_pearlite_neumann_panels.png", dpi=150)
plt.show()
print("Saved: 10_2_pearlite_neumann_panels.png")

# Save GIF 
os.makedirs("pearlite_neumann_frames", exist_ok=True)
image_frames = []

for k in range(0, Nt):
    plt.contourf(x_vec, y_vec, T[:, :, k])
    plt.colorbar()
    plt.title(f'Concentration profile at timestep {k + 1}')
    fname = f"pearlite_neumann_frames/frame_{k + 1:04d}.jpg"
    plt.savefig(fname)
    plt.clf()
    image_frames.append(Image.open(fname))

image_frames[0].save(
    "10_2_pearlite_neumann_evolution.gif",
    format='GIF',
    append_images=image_frames[1:],
    save_all=True,
    duration=t,
    loop=0
)
print("Saved: 10_2_pearlite_neumann_evolution.gif")
