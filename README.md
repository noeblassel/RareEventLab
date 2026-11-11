# Navigating the design parameter sea for rare-event sampling algorithms

Practical session for the "Learning the Dynamics of Biological Macromolecules" research school, Cargèse, November 2026 (author: Noé Blassel).

A Julia notebook to experiment with the hyperparameters of dynamical rare-event sampling methods (Adaptive Multilevel Splitting, Parallel Replica) on a 2D overdamped Langevin toy model, to explore how they affect correctness and efficiency of numerical results.

## Setup

Linux or macOS (on Windows, use WSL):

```bash
./setup.sh [number_of_threads]   # default: number of available logical cores
```

This downloads Julia and the packages into local folders (`julia/`, `.julia-depot/`) and installs a Jupyter kernel named **Julia RareEventLab**. Then open `notebook.ipynb` with that kernel.

No Jupyter installed? `setup.sh` prints the command to launch one from Julia.

To remove everything: `./uninstall.sh`.

## Contents

- `notebook.ipynb` — the practical session
- `RareEventLab/` — the supporting Julia package (dynamics, AMS, ParRep, plotting)
- `figures/` — figures used in the notebook