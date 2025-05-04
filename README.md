# Simulation of Thin Film Drying with Moving Boundaries

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
```

Or install manually:

```bash
pip install -r requirements.txt
```

## ▶️ Usage

```bash
python run_simulation.py
```

## 📁 Folder Structure

- `src/`: Source code (PDE solver, boundary evolution, utilities)
- `examples/`: Example scripts to reproduce figures and tables
- `data/`: Input data (if needed)
- `results/`: Output plots and simulations
- `docs/`: Extended documentation (optional)

## 🧪 Reproducibility

To reproduce the main results of the publication, run:

```bash
python examples/reproduce_main_figure.py
```

## 📜 License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.

## 📄 Citation

Please cite this code as follows:

```
@misc{Heyd2025,
  author    = {Rodolphe Heyd},
  title     = {Simulation Code for Drying of Thin Films with Moving Boundaries},
  year      = {2025},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.1234567},
  url       = {https://doi.org/10.5281/zenodo.1234567}
}
```

Or via the DOI: [https://doi.org/10.5281/zenodo.1234567](https://doi.org/10.5281/zenodo.1234567)
