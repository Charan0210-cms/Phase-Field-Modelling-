"""
5.2  1D Heat Equation – FDM and Dirichlet BC with Zero Boundary Condition
Constant initial condition (spike at midpoint) and zero BC at both ends.
Parameters are entered interactively (or set below as defaults).
"""

import numpy as np
import matplotlib.pyplot as plt


M  = int(input("Enter No. of grid points for space interval [3000]: ") 
N  = int(input("Enter No. of grid points for time  interval [3000]: ") 
x0 = float(input("Enter starting space grid point [0]: ")  
xL = float(input("Enter final   space grid point  [100]: ") 
t0 = float(input("Enter starting time  grid point [0]: ")   
tF = float(input("Enter final   time   grid point  [10]: ") 
D  = float(input("Enter Thermal Diffusivity value  [0.005]: ") 
bc1 = float(input("Enter boundary value at x0 for any t [0]: ") 
bc2 = float(input("Enter boundary value at xL for any t [0]: ") 

# Grid 
t_vec = np.linspace(t0, tF, N)
x_vec = np.linspace(x0, xL, M)
dx = x_vec[1] - x_vec[0]
dt = t_vec[1] - t_vec[0]
alpha = D * dt / dx**2

print(f"Stability parameter alpha = {alpha:.6f}  ({'stable' if alpha <= 0.5 else 'UNSTABLE'})")

T = np.zeros((M, N))

# Boundary conditions
T[0, :]    = bc1
T[M - 1, :] = bc2

# Initial condition: spike at midpoint (index 1500 assumes M≈3000)
mid = M // 2
T[mid, 0] = 100

#  Time-stepping 
for n in range(0, N - 1):
    for i in range(1, M - 1):
        T[i, n + 1] = alpha * (T[i - 1, n] + T[i + 1, n]) + (1 - 2 * alpha) * T[i, n]

#  Plot selected time snapshots 
snapshots = [0, N // 3, 2 * N // 3, N - 1]
labels    = [str(s + 1) for s in snapshots]

plt.figure(figsize=(9, 5))
for s, lbl in zip(snapshots, labels):
    plt.plot(x_vec, T[:, s], label=lbl)

plt.xlabel('Space')
plt.ylabel('Temperature')
plt.title('1D Heat Equation – FDM & Dirichlet BC (Zero BC)')
plt.legend(title='Time step')
plt.tight_layout()
plt.savefig('5_2_dirichlet_zero_bc_2D.png', dpi=150)
plt.show()
print("Saved: 5_2_dirichlet_zero_bc_2D.png")
