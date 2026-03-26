"""
5.3  1D Heat Equation – FDM and Dirichlet BC with Non-Zero Boundary Condition
Constant initial condition (spike at midpoint) and non-zero BC at both ends.
"""

import numpy as np
import matplotlib.pyplot as plt

M  = int(input("Enter No. of grid points for space interval [3000]: "))
N  = int(input("Enter No. of grid points for time  interval [3000]: ")) 
x0 = float(input("Enter starting space grid point [0]: ")) 
xL = float(input("Enter final   space grid point  [100]: ")) 
t0 = float(input("Enter starting time  grid point [0]: "))  
tF = float(input("Enter final   time   grid point  [10]: ")) 
D  = float(input("Enter Thermal Diffusivity value  [0.005]: ")) 
bc1 = float(input("Enter boundary value at x0 for any t [30]: ")) 
bc2 = float(input("Enter boundary value at xL for any t [50]: ")) 

#Grid
t_vec = np.linspace(t0, tF, N)
x_vec = np.linspace(x0, xL, M)
dx = x_vec[1] - x_vec[0]
dt = t_vec[1] - t_vec[0]
alpha = D * dt / dx**2

print(f"Stability parameter alpha = {alpha:.6f}  ({'stable' if alpha <= 0.5 else 'UNSTABLE'})")

T = np.zeros((M, N))

# Boundary conditions (non-zero)
T[0, :]     = bc1
T[M - 1, :] = bc2

# Initial condition: spike at midpoint
mid = M // 2
T[mid, 0] = 100

# Time-stepping 
for n in range(0, N - 1):
    for i in range(1, M - 1):
        T[i, n + 1] = alpha * (T[i - 1, n] + T[i + 1, n]) + (1 - 2 * alpha) * T[i, n]

#  Plot 
snapshots = [0, N // 3, 2 * N // 3, N - 1]

plt.figure(figsize=(9, 5))
for s in snapshots:
    plt.plot(x_vec, T[:, s], label=str(s + 1))

plt.xlabel('Space')
plt.ylabel('Temperature')
plt.title('1D Heat Equation – FDM & Dirichlet BC (Non-Zero BC)')
plt.legend(title='Time step')
plt.tight_layout()
plt.savefig('5_3_dirichlet_nonzero_bc_2D.png', dpi=150)
plt.show()
print("Saved: 5_3_dirichlet_nonzero_bc_2D.png")
