# Finite Difference Method for PDEs & Phase Field Modelling

This repo contains the work I did during my internship (June–July 2021) at NIT Tiruchirappalli, under Dr Prince Gideon Kubendran Amos from the Department of Metallurgical and Materials Engineering.

The internship was focused on learning how to numerically solve Partial Differential Equations using the Finite Difference Method, and eventually applying it to model the dissolution of Pearlite using Fick's Law of Diffusion.

### heat_equation_1D
Solving the 1D heat equation `u_t = c² u_xx` under different boundary conditions.

 File and its Description 

| `dirichlet_bc.py` | Parabolic IC `u(x,0) = 4x - 4x²`, zero BC at both ends. 3D surface plot output. |
| `dirichlet_bc_zero.py` | Spike IC, zero Dirichlet BC at both ends |
| `dirichlet_bc_nonzero.py` | Spike IC, non-zero Dirichlet BC at both ends |
| `neumann_bc_zero_flux_1.py` | Spike IC, zero flux Neumann BC |
| `neumann_bc_nonzero_flux_1.py` | Spike IC, non-zero flux Neumann BC |
| `neumann_bc_zero_flux_2.py` | Step function IC, zero flux Neumann BC |
| `neumann_bc_nonzero_flux_2.py` | Step function IC, non-zero flux Neumann BC |
| `1D_mixed_bc.py` | Dirichlet on left end, Neumann on right end |

### heat_equation_2D
Extending to the 2D heat equation `u_t = D(u_xx + u_yy)`.

File and its Description 

| `dirichlet_bc_code1.py` | Fixed boundary temperatures, animated heatmap output (saves as GIF) |
| `dirichlet_bc_code2.py` | Interactive inputs, contour plot output at each timestep |
| `mixed_bc.py` | Dirichlet on top/bottom, Neumann on left/right |

### pearlite_dissolution
Modelling how carbon diffuses during the dissolution of Pearlite into Austenite, using Fick's Law `dC/dt = D d²C/dx²`. The domain represents two pearlite lamellae — ferrite layers (0.02 wt% C) alternating with cementite layers (6.67 wt% C).

File and its Description 

| `dirichlet_bc.py` | Carbon concentration evolution with Dirichlet BC |
| `neumann_bc.py` | Same but with Neumann (flux) BC on all edges |


## Some context on the physics

Pearlite is a microstructure found in steel — it consists of alternating layers of ferrite and cementite. When you heat it above the eutectoid temperature, carbon starts diffusing out of the cementite into the ferrite, and the structure transforms into austenite. This project numerically models that carbon diffusion process.

The Finite Difference Method works by discretising the spatial and time domains into a grid, then approximating derivatives using Taylor series expansions at each grid point. The explicit Euler method was used for time-stepping, which requires the stability condition:

- 1D: `dt ≤ dx² / (2D)`
- 2D: `dt ≤ dx² / (4D)`


## References

Causon & Mingham — *Introductory Finite Difference Methods for PDEs*
Langtangen & Linge — *Finite Difference Computing with PDEs*
Mazumder — *Numerical Methods for PDEs: Finite Difference and Finite Volume Methods*
Rajan, Sharma & Sharma — *Heat Treatment: Principles and Techniques*

Internship submitted to: Dr Prince Gideon Kubendran Amos, NIT Tiruchirappalli  
Submitted by: Charan MS (2019110004), College of Engineering, Guindy, Anna University
