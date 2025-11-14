# MATLAB Spectral-Galerkin thermal model for (cylindrical) battery cells

## Overview

This project implements a 2D (r-z coordinates) thermal model for a cylindrical battery cell. The model uses the Spectral-Galerkin method with Chebyshev polynomials to solve the heat equation PDE.

The core of the model is a state-space representation:

* **Continuous-time:**  $G \dot{x}(t) = A   x(t) + B   u(t) + F   w(t)$
*   **Discrete-time:**        $x(k+1) = A_d x(k) + B_d u(k) + F_d w(k)$

It uses a "lifting" technique by splitting the temperature $T$ into a homogeneous solution $T_h$ (the states) and a particular solution $T_p$ (driven by boundary conditions).



![Schematic of the cylindrical battery model](./Images/Battery_diag.jpg)


---

## 📋 Requirements

* MATLAB R2023b or newer
* Control System Toolbox (for `expm`)

---

## 🚀 How to Run

1.  Open MATLAB.
2.  Add the functions to your path:
    ```matlab
    addpath(genpath('./Functions_ThermalMdl'));
    ```
3.  Run the main simulation script:
    ```matlab
    Main_ThermMdl_Simulation
    ```
4.  The script will generate a plot of key battery temperatures over time.

---

## 📁 File Descriptions

* `Functions_ThermalMdl/`: Folder containing all helper functions.
    * `clenCurt_Quad.m`: Provides weights for Clenshaw-Curtis quadrature.
    * `chebdif.m`: Provides Chebyshev differentiation matrices.
    * `Cheby_BasisFxn.m`: Generates the Chebyshev basis functions. 
    * `Parameters_ThermalMdl.m`: Defines all physical constants (density $\rho$, conductivity      $k_r$, geometry $H$, $R$, etc.).
    * `Inputs_ThermalMdl.m`: Defines simulation inputs (time vector `p.t`, heat generation `p.Q`).
    * `SpectralGalerkin_ThermalMdl.m`: The core function. Builds the continuous-time matrices (`A`, `B`, `G`, `C`).
    * `MdlSetup_ThermalMdl.m`: Constructs the discrete-time state-space model (`Ad`, `Bd`, `C`, `Du`).
    * `Main_ThermMdl_Simulation.m`: The main executable script. It sets up parameters, runs the time-step simulation, and plots the results.
    

---

## 📚 Citation 

**If you use this thermal model (code) in your research, please cite the following paper:**

> **"Thermal Modelling of Battery Cells for Optimal Tab and Surface Cooling Control"**
> Godwin K. Peprah, Yicun Huang, Torsten Wik, Faisal Altaf, Changfu Zou
> * arXiv:2409.08974v2*
> https://doi.org/10.48550/arXiv.2409.08974

## Acknowledgments

The methodology in this model was adapted from:
R. R. Richardson's "On-board monitoring of 2-D spatially-resolved temperatures in cylindrical lithium-ion batteries: Part I. Low-order thermal modelling" paper.

---

## 📄 License

This project is licensed under the MIT License - see the `LICENSE` file for details.