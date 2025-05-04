# ThinFilm-Drying

This repository contains numerical code for simulating the drying of a thin aqueous film with a dynamic evaporation front. The model couples a nonlinear diffusion PDE with a time-dependent ODE governing the moving boundary.

## 🔍 Description

- Implicit finite difference scheme
- Newton-Raphson iteration for nonlinearity
- Time-dependent domain evolution
- Reproducible examples

## 📦 Dependencies

Create the environment using:

```bash
conda env create -f environment.yml
conda activate thinfilm-drying
