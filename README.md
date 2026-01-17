# Mixedstate_topo Repository

This repository contains **raw data and scripts** generated from DQMC simulations and numerical calculations used in the paper:  

**_Diagnosis of mixed-state topological phases in strongly correlated systems via disorder parameters_**.  

The repository includes the following:

---

## 1. freecase_code

The `freecase_code` folder provides tools to generate **1D and 2D disorder parameters** and  **second derivatives of DP at twisted angle π  **, as well as to plot related results. Specifically:  

- **1D model**: SSH model  
- **2D models**: QWZ model and Haldane model  
- Computes and visualizes the disorder parameter (DP) and its second derivatives  
- Generates plots of relevant **band structures**  

---

## 2. Strongly correlated models

The repository also contains **DQMC simulation results** for interacting systems in `Kane_Mele` and `HaldaneHubbard` folders:  

- **Models included**:  
  - Kane-Mele Hubbard model  
  - Haldane Hubbard model  

- **Parameters considered**:  
  - Lattice size (`Lx`, `Ly`)  
  - Spin-orbit coupling (`λ`)  
  - Inverse temperature (`β`)  
  - Interaction strength (`U`)  

- **Data contents**:  
  - Second derivatives of DP at twisted angle π  
  - Raw data under various parameter sets  

- **Data format**:  
Each data file is a text file with three columns:  

| Column | Description |
|--------|-------------|
| 1      | `Lx` (lattice size along x) |
| 2      | DP second derivative at twisted angle π |
| 3      | Error of the DP second derivative |

These raw data can be used to generate the **topological indicator** (`\mathcal{F}`) and ΔLx plots in the paper.  












