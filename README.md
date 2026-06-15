[![arXiv](https://img.shields.io/badge/arXiv-2601.22459-b31b1b.svg)](https://arxiv.org/abs/2601.22459)

# Multilplex social contagion contact based simulation in Fortran, orchestrated from Python.

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Requirements & Installation](#requirements--installation)
- [Quick Start](#quick-start)
- [sim_manager.py – All Commands](#sim_managerpy--all-commands)
  - [run](#run--single-simulation)
  - [beta-loop](#beta-loop--sweep-inpr)
  - [eta-loop](#eta-loop--sweep-eta)
  - [density](#density--density-sweep)
  - [corr-loop](#corr-loop--correlation-run)
- [All CLI Flags (Shared Parameters)](#all-cli-flags-shared-parameters)
- [Supporting Scripts](#supporting-scripts)
  - [matrixGeneration/adjacencymat.py](#matrixgenerationadjacencymatpy)
  - [correlationCalc/corrtrue.py](#correlationcalccorrtrruepy)
  - [simulation/ (Fortran core)](#simulation-fortran-core)
  - [mstereq/ – Master Equation](#mstereq--master-equation-numerical-integration)
  - [experiment/ – Sudden Beta-Change](#experiment--sudden-beta-change-experiment)
- [Output Files & Data](#output-files--data)
- [Temporary Files](#temporary-files)
- [Folder Reference](#folder-reference)

---

## Overview

This project implements a **multiplex social contagion contact-based model. The simulation core is written in **Fortran** for performance and is compiled and invoked automatically by a Python orchestration layer (`sim_manager.py`). All workflows — single runs, parameter sweeps, and correlation analyses — are accessible through a single entry point.

The model studies how intra-layer contagion ratio (`INPR`) and interlayer contagion ratio (`ETA`) in order to study the regime transition signatiures in correlation. Results are stored as binary data and post-processed into HDF5 correlation files and PDF/PNG plots.

For more information, see the [References](#references).

---

## Repository Structure

```
TFM/
├── sim_manager.py              ← single entry point for ALL workflows
│
├── simulation/                 ← Fortran source + compiled binary
│   ├── sync.f90                ← main simulation logic
│   ├── r1279.f90               ← RNG (LFSR generator, period 2^1279−1)
│   ├── ran2.f                  ← alternative RNG (Numerical Recipes ran2)
│   └── simulation(.exe)        ← compiled binary (auto-built on first run)
│
├── matrixGeneration/
│   └── adjacencymat.py         ← generates two-layer adjacency matrix (ER or RRN)
│
├── correlationCalc/
│   └── corrtrue.py             ← reads output.bin, computes lag correlations → HDF5
│
├── experiment/                 ← experiment sudden beta change in simulation
├── images generator/           ← scripts for generating result plots
├── images/                     ← output plots (PDF/PNG)
├── Data/                       ← result data (TXT/NPY)
├── mstereq/                    ← mean-field / master-equation utilities
│
├── .gitattributes
├── .gitignore
└── README.md
```

> **Note:** `INITIAL.txt`, `adjacency_matrix.txt`, and `output.bin` are temporary files created during a run. They are listed in `.gitignore` and should not be committed.

---

## Requirements & Installation

### Python dependencies

```bash
pip install numpy networkx h5py joblib
```

### Fortran compiler

`gfortran` must be on your `PATH`. The simulation binary is compiled automatically on the first run. To compile manually:

```bash
cd simulation/
gfortran -O2 -o simulation sync.f90 r1279.f90 ran2.f
```

### Verifying setup

```bash
python sim_manager.py --help
```

---

## Quick Start

```bash
# Single simulation
python sim_manager.py run --ETA 25 --INPR 0.05

# Sweep inter-layer coupling (beta-loop)
python sim_manager.py beta-loop --ETA 25 --repeats 3

# Sweep noise parameter (eta-loop)
python sim_manager.py eta-loop --INPR 0.034 --repeats 3

# Density sweep over ETA
python sim_manager.py density --INPR 0.034 --N 5000 --steps 50 --sweep E

# Single-INPR correlation run (10 repeats)
python sim_manager.py corr-loop --ETA 25 --INPR 0.05
```

---

## sim_manager.py – All Commands

`sim_manager.py` is the **single entry point** for all simulation workflows. It handles compilation of the Fortran binary if needed, matrix generation, invoking the simulation, and post-processing results.

---

### `run` – Single Simulation

Runs a single simulation with the given parameters and computes the output correlations.

```bash
python sim_manager.py run [options]
```

**What it does:**
1. Generates (or reuses) the two-layer adjacency matrix via `adjacencymat.py`.
2. Writes `INITIAL.txt` with all simulation parameters.
3. Invokes the compiled Fortran binary, producing `output.bin`.
4. Calls `corrtrue.py` to compute lag-0 and lag-k correlations, saving results to `DATA/`.

**Example:**
```bash
python sim_manager.py run --ETA 25 --INPR 0.05 --N 1000 --STIME 40000
```

---

### `beta-loop` – Sweep INPR

Sweeps the inter-layer coupling parameter `INPR` (beta) over a logarithmic range, running multiple repeats at each value.

```bash
python sim_manager.py beta-loop [options]
```

**What it does:**
- Iterates over `x = np.logspace(-2, -1, 20)` (20 values from 0.01 to 0.1 by default).
- For each value of `INPR`, repeats the simulation `--repeats` times.
- Appends results to `DATA/eta{eta_str}_mu{mu_str}_N{N}/final_corr{1..6}.txt`.

**To change the sweep range**, edit the `logspace` call inside the `beta-loop` handler in `sim_manager.py`:
```python
x = np.logspace(-2, -1, 20)  # modify start, stop, or num points here
```

**Example:**
```bash
python sim_manager.py beta-loop --ETA 25 --repeats 3
python sim_manager.py beta-loop --ETA 1 --repeats 5 --N 2000
```

**Additional flag:**

| Flag | Default | Description |
|------|---------|-------------|
| `--repeats` | 3 | Number of independent repeats per INPR value |

---

### `eta-loop` – Sweep ETA

Sweeps the noise parameter `ETA` over a logarithmic range, for a fixed `INPR`, running multiple repeats at each value.

```bash
python sim_manager.py eta-loop [options]
```

**What it does:**
- Iterates over `x = np.logspace(0, 1.778, 20)` (20 values, roughly 1 to 60) by default.
- For each value of `ETA`, repeats the simulation `--repeats` times.
- Appends results to `DATA/eta{eta_str}_mu{mu_str}_N{N}/final_corr{1..6}.txt`.

**To change the sweep range**, edit the `logspace` call inside the `eta-loop` handler in `sim_manager.py`:
```python
x = np.logspace(0, 1.778, 20)  # modify to change the ETA range
```

**Example:**
```bash
python sim_manager.py eta-loop --INPR 0.034 --repeats 3
python sim_manager.py eta-loop --INPR 0.05 --repeats 10 --N 500
```

**Additional flag:**

| Flag | Default | Description |
|------|---------|-------------|
| `--repeats` | 3 | Number of independent repeats per ETA value |

---

### `density` – Density Sweep

Sweeps either `ETA` (`E`) or `INPR`/beta (`B`) and records the mean per-layer node density (fraction of active nodes) at each step.

```bash
python sim_manager.py density [options]
```

**What it does:**
- Iterates over a range of the chosen parameter (`--sweep E` for ETA, `--sweep B` for INPR/beta).
- At each step, runs the simulation and records the mean density of active nodes per layer.
- Saves a `.npy` array with per-layer mean densities to `DATA/`.

**Example:**
```bash
# Sweep ETA, N=5000 nodes, 50 parameter steps
python sim_manager.py density --INPR 0.034 --N 5000 --steps 50 --sweep E

# Sweep beta (INPR), default N
python sim_manager.py density --ETA 25 --steps 30 --sweep B
```

**Additional flags:**

| Flag | Default | Description |
|------|---------|-------------|
| `--steps` | 50 | Number of points in the parameter sweep |
| `--sweep` | — | `E` to sweep ETA, `B` to sweep INPR/beta (required) |

---

### `corr-loop` – Correlation Run

Runs 10 independent simulations at a fixed `(ETA, INPR)` pair and computes lag correlations for each, useful for averaging over stochastic realisations.

```bash
python sim_manager.py corr-loop [options]
```

**What it does:**
- Runs the simulation 10 times with identical parameters.
- Computes lag-0 and lag-k correlations via `corrtrue.py` for each run.
- Results are stored in `DATA/` for subsequent averaging or plotting.

**Example:**
```bash
python sim_manager.py corr-loop --ETA 25 --INPR 0.05
python sim_manager.py corr-loop --ETA 1 --INPR 0.034 --N 2000
```

---

## All CLI Flags (Shared Parameters)

All `sim_manager.py` commands accept the following flags:

| Flag | Default | Description |
|------|---------|-------------|
| `--N` | 1000 | Network size (nodes per layer) |
| `--DIM` | 2 | Number of network layers |
| `--STIME` | 40000 | Total simulation time steps |
| `--ETA` | 25.0 | Noise parameter (η); controls stochasticity of spin flips |
| `--INPR` | 0.05 | Inter-layer coupling strength (β); links the two layers |
| `--M` | 1000 | Record state every M steps |
| `--GAMMA1` | 100000 | Degree parameter for layer 1 (network connectivity) |
| `--GAMMA2` | 100000 | Degree parameter for layer 2 (network connectivity) |
| `--DISCTIME` | 20000 | Warm-up (discard) steps before recording begins |
| `--MU` | 0.05 | Death/rewiring rate μ |

---

## Supporting Scripts

These scripts are called automatically by `sim_manager.py` but can also be run directly.

---

### `matrixGeneration/adjacencymat.py`

Generates the two-layer network adjacency matrix used as input to the simulation.

**What it does:**
- Generates an Erdős–Rényi (ER) random graph adjacency matrix for each layer by default.
- Can be configured (by editing comments in the script) to produce a **Regular Random Network (RRN)** instead.
- Writes the result to `adjacency_matrix.txt` in the repo root.

**Usage (direct):**
```bash
python matrixGeneration/adjacencymat.py --N 1000 --GAMMA1 100000 --GAMMA2 100000
```

**To switch to RRN:** open `adjacencymat.py` and follow the commented instructions to change the graph generation mode.

---

### `correlationCalc/corrtrue.py`

Reads the raw binary simulation output and computes lag-k cross-correlations between node activity time series.

**What it does:**
- Reads `output.bin` produced by the Fortran simulation.
- Computes pairwise (or mean-field) lag correlations across both layers.

**Usage (direct):**
```bash
python correlationCalc/corrtrue.py
```

Adjust file paths and parameters at the top of the script as needed.

---

### `simulation/` (Fortran core)

The computational heart of the project. The Fortran source files are:

| File | Role |
|------|------|
| `sync.f90` | Main simulation: Ising-like spin dynamics on a two-layer network with death/rewiring (rate μ). Reads `INITIAL.txt` and `adjacency_matrix.txt`, writes `output.bin`. |
| `r1279.f90` | High-quality LFSR pseudorandom number generator (period 2^1279−1), used for large-scale runs. |
| `ran2.f` | Alternative RNG from Numerical Recipes (ran2 algorithm), available as a fallback. |
| `simulation` / `simulation.exe` | Compiled binary. Created automatically by `sim_manager.py` on first run. |

**Manual compilation:**
```bash
cd simulation/
gfortran -O2 -o simulation sync.f90 r1279.f90 ran2.f
```

---

### `mstereq/` – Master Equation Numerical Integration

Provides a **mean-field analytical reference** for the simulation by numerically integrating the master equation of the model.

**What it does:**
- Evolves the probability distribution over system states forward in time by evolving the master equation numerically.
- Because this is a deterministic numerical integration rather than a noisy simulation, it produces smooth, noise-free trajectories that serve as a theoretical baseline.
- Results can be overlaid on top of simulation data (from `sim_manager.py`) to assess how well the mean-field approximation captures the true dynamics at finite network sizes.

**Typical use:**
- Run the master equation integrator for a given `(ETA, INPR, MU)` parameter set.
- Compare the resulting steady-state densities or time-evolution curves against the Fortran simulation output.

**Usage:**
```bash
python mstereq/<script_name>.py
```

Edit parameters (ETA, INPR, MU, integration step size, total time) directly at the top of the script.

---

### `experiment/` – Sudden Beta-Change Experiment

Simulates a **sudden, instantaneous shift in the intra-layer coupling** (β / `INPR`) at a randomly chosen time during the run. This models real-world scenarios where the coupling between two network layers changes abruptly.

**What it does:**
- Runs the Fortran simulation normally up to a random (or fixed) switching time `t*`.
- At `t*`, the intra-layer coupling is changed from an initial value β₀ to a new value β₁ without interrupting the spin dynamics.
- Records the system's response — how node activity, correlations, and synchronisation evolve before and after the switch.
- Useful for studying transient dynamics and recovery times after a structural perturbation.

**Typical use:**
```bash
python experiment/<script_name>.py
```

Edit `beta_before`, `beta_after`, `switch_time` (or set `switch_time = None` for a random switch), and other parameters directly inside the script.

**Output:** time-series of per-layer densities and correlations split into pre- and post-switch windows, saved to `DATA/` or `images/` as configured in the script.

---

## References

- Tey, J. H., & Cozzo, E. (2026). *Correlation-Based Diagnostics of Social Contagion Dynamics in Multiplex Networks.* arXiv preprint. [arXiv:2601.22459](https://arxiv.org/abs/2601.22459)

- Tey, J. H., & Cozzo, E. (2024). *Eigenvector Localization and Universal Regime Transitions in Multiplex Networks: A Perturbative Approach.* arXiv preprint. [arXiv:2408.04784](https://arxiv.org/abs/2408.04784)

---

## Output Files & Data

| Path | Content |
|------|---------|
| `DATA/eta{η}_mu{μ}_N{N}/final_corr{1..6}.txt` | Appended correlation results from beta-loop and eta-loop sweeps |
| `DATA/*.npy` | Per-layer density arrays from density sweeps |
| `images/*.pdf` / `images/*.png` | Result plots (correlation matrices, density curves, etc.) |

---

## Temporary Files

The following files are generated during a run and are excluded from version control via `.gitignore`:

| File | Created by | Description |
|------|-----------|-------------|
| `INITIAL.txt` | `sim_manager.py` | Simulation parameters passed to Fortran binary |
| `adjacency_matrix.txt` | `adjacencymat.py` | Two-layer network connectivity |
| `output.bin` | Fortran `simulation` binary | Raw time-series state output |

Do **not** commit these files.

---

## Folder Reference

| Folder | Contents |
|--------|----------|
| `simulation/` | Fortran source code and compiled binary |
| `matrixGeneration/` | Network (adjacency matrix) generation scripts |
| `correlationCalc/` | Post-processing: lag-correlation computation |
| `experiment/` | Sudden-beta-change experiment: simulates an abrupt shift in inter-layer coupling at a random time, mimicking real-network perturbations |
| `images generator/` | Scripts for producing plots from DATA/ |
| `images/` | Output plots (PDF/PNG) |
| `Data/` | Output data (TXT/NPY) |
| `mstereq/` | Numerical integration of the master equation for the model; provides mean-field reference trajectories to compare against simulation results |


