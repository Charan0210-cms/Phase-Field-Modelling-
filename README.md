# Numerical Modeling of PDEs using Finite Difference Method

## Overview
This repository documents my internship work on solving **partial differential equations (PDEs)** using the **Finite Difference Method (FDM)**, along with introductory **phase-field modeling** concepts and a materials-science application focused on the **dissolution of pearlite**.

The work covers:
- Classification of PDEs
- Initial and boundary conditions
- Finite Difference Method fundamentals
- 1D heat equation
- 2D heat equation
- Dirichlet, Neumann, and mixed boundary conditions
- Introductory phase-field modeling
- Numerical modeling of pearlite dissolution using **Fick's Law of Diffusion**

## Internship Context
This work was completed as part of an internship carried out from **1 June to 31 July 2021** under **Dr. Prince Gideon Kubendran Amos** at **National Institute of Technology, Tiruchirappalli (NITT)**.

## Project Objectives
The main goals of this work were:
1. To understand the basics of PDEs and their classifications.
2. To learn how initial and boundary conditions affect PDE solutions.
3. To implement numerical solutions using the Finite Difference Method.
4. To model temperature evolution in 1D and 2D heat conduction problems.
5. To explore diffusion-based microstructure evolution through pearlite dissolution modeling.

## Topics Covered

### 1. Partial Differential Equations (PDEs)
The report begins with an introduction to PDEs and their importance in physics, engineering, heat transfer, diffusion, fluid flow, and related applications.

### 2. Classification of PDEs
PDEs are categorized into:
- **Elliptic**
- **Parabolic**
- **Hyperbolic**

The report discusses the standard classification based on the discriminant \(B^2 - 4AC\).

### 3. Initial and Boundary Conditions
The project explains the role of:
- **Initial Conditions (ICs)**
- **Boundary Conditions (BCs)**

Boundary conditions studied include:
- **Dirichlet boundary condition**
- **Neumann boundary condition**
- **Mixed boundary condition**

### 4. Finite Difference Method (FDM)
The numerical approach is based on:
- Spatial and temporal discretization
- Taylor series expansion
- Forward, backward, and central difference approximations
- Explicit time-stepping schemes

## Numerical Problems Implemented

### 1D Heat Equation
The 1D heat equation was solved using FDM under multiple boundary-condition cases:

#### Dirichlet Boundary Condition
- Zero boundary condition
- Non-zero boundary condition

#### Neumann Boundary Condition
- Zero flux at both boundaries
- Non-zero flux at both boundaries
- Constant initial condition cases
- Step-function initial condition cases

#### Mixed Boundary Condition
- Dirichlet at one end
- Neumann at the other end

### 2D Heat Equation
The 2D heat equation was solved using an explicit finite-difference method.

Cases include:
- **2D heat equation with Dirichlet BC**
- **Alternative implementation (Code 2)**
- **2D heat equation with mixed Dirichlet-Neumann BC**

The results were visualized using:
- Heat maps
- Time-evolution plots
- Animated temperature evolution

## Phase-Field Modeling
The internship also included an introduction to **phase-field modeling**, particularly for representing microstructural evolution in materials systems.

The report discusses:
- Why phase-field methods are useful in materials science
- Their connection to PDE-based numerical modeling
- Example simulations showing microstructure evolution

## Application: Dissolution of Pearlite
A key materials-science application in this work is the **numerical modeling of pearlite dissolution**.

### Background
Pearlite is a eutectoid mixture of:
- **Ferrite**
- **Cementite**

It has a lamellar microstructure with alternating layers. The report explains how carbon diffusion and phase transformation occur as pearlite transforms toward austenite at elevated temperatures.

### Governing Equation
The model uses **Fick's Law of Diffusion**:

\[
\frac{dC}{dt} = D \frac{d^2 C}{dx^2}
\]

where:
- \(C\) = concentration
- \(D\) = diffusion coefficient

### Numerical Modeling Cases
The pearlite dissolution problem was implemented for:
- **Dirichlet boundary condition**
- **Neumann boundary condition**

The simulations track the **evolution of carbon concentration** over time.

## Key Learning Outcomes
Through this work, I gained experience in:
- Understanding PDE fundamentals
- Applying numerical discretization techniques
- Implementing explicit FDM schemes
- Working with different types of boundary conditions
- Visualizing numerical solutions in 2D and 3D
- Connecting diffusion equations to materials-science problems
- Modeling pearlite dissolution conceptually through concentration evolution

## Suggested Repository Structure
Since the original Python files are currently unavailable, this repository can be organized as follows for documentation purposes:

```text
.
├── README.md
├── report/
│   └── PF-Internship report.pdf
├── images/
│   ├── 1d_heat_equation_results.png
│   ├── 2d_heatmap_results.png
│   └── pearlite_dissolution_results.png
└── notes/
    └── methodology.md
```

## Software / Tools
Based on the report, the numerical work appears to have involved:
- **Python**
- Numerical arrays / matrix-based computation
- Plotting for 2D and 3D visualization

If the code is rebuilt later, likely dependencies may include:
- `numpy`
- `matplotlib`

## Future Improvements
Possible future additions to this repository:
- Rebuild the original Python scripts
- Add cleaned and commented code for each numerical case
- Include stability-condition explanations for explicit schemes
- Add animations and output plots
- Extend the work to implicit methods and more advanced phase-field simulations

## Reference
This README is reconstructed from the internship report and is intended to document the scope of the work even though the original code files are currently missing.

## Author
**Charan MS**  
College of Engineering, Guindy, Anna University
