"""
5.1  1D Heat Equation – FDM and Dirichlet BC
Initial condition : u(x, 0) = 4x - 4x²
Boundary condition: u(0, t) = u(L, t) = 0
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Grid parameters 
M = 40      # number of spatial grid points
N = 200     # number of time grid points

x0, xL = 0, 1          # spatial domain [0, 1]
t0, tF = 0, 0.2        # time domain    [0, 0.2]

dx = (xL - x0) / (M - 1)
dt = (tF - t0) / (N - 1)

D = 0.3                 # thermal diffusivity
a = D * dt / dx**2      # stability parameter (must be ≤ 0.5)

print(f"Stability parameter a = {a:.4f}  ({'stable' if a <= 0.5 else 'UNSTABLE – reduce dt or increase dx'})")

#  Grid arrays 
tspan = np.linspace(t0, tF, N)
xspan = np.linspace(x0, xL, M)

#  Solution matrix  U[space, time] 
U = np.zeros((M, N))

# Initial condition: u(x, 0) = 4x - 4x²
U[:, 0] = 4 * xspan - 4 * xspan**2

# Dirichlet boundary conditions
U[0, :]  = 0   # u(0, t) = 0
U[-1, :] = 0   # u(L, t) = 0

#  Time-stepping (explicit Euler)
for k in range(0, N - 1):
    for i in range(1, M - 1):
        U[i, k + 1] = a * U[i - 1, k] + (1 - 2 * a) * U[i, k] + a * U[i + 1, k]

#  3-D surface plot 
T_grid, X_grid = np.meshgrid(tspan, xspan)

fig = plt.figure(figsize=(10, 6))
ax  = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(T_grid, X_grid, U, cmap='coolwarm')
fig.colorbar(surf, ax=ax, shrink=0.5, label='Temperature u')

ax.set_xlabel('Time')
ax.set_ylabel('Space')
ax.set_zlabel('U')
ax.set_title('1D Heat Equation – FDM & Dirichlet BC\nu(x,0)=4x−4x²,  u(0,t)=u(1,t)=0')
plt.tight_layout()
plt.savefig('5_1_dirichlet_bc_3D.png', dpi=150)
plt.show()
print("Saved: 5_1_dirichlet_bc_3D.png")
