# TFM – Simulation Manager

Multilayer Ising-like synchronisation simulation in Fortran, orchestrated from Python.

---

## Repository Structure

```
TFM/
├── sim_manager.py              ← single entry point for ALL workflows
│
├── simulation/                 ← Fortran source + compiled binary
│   ├── sync.f90
│   ├── r1279.f90
│   ├── ran2.f
│   └── simulation(.exe)        (compiled automatically on first run)
│
├── matrixGeneration/
│   └── adjacencymat.py         ← generates two-layer Erdős–Rényi (or RRN if you change some comments) adjacency matrix
│
├── correlationCalc/
│   └── corrtrue.py             ← reads output.bin, computes lag correlations → HDF5
│
├── images/                     ← result plots (PDF/PNG)
├── Data/                     ← results data (TXT/NPY)
│
├── experiment/
│   └── 
│
└── README.md
```

> **Note:** `INITIAL.txt`, `adjacency_matrix.txt`, and `output.bin` are
> temporary files created and moved around during a run. They are listed in
> `.gitignore` and should not be committed.

---

## Quick Start

### Requirements

```bash
pip install numpy networkx h5py joblib
# gfortran must be on your PATH
```

### Run a single simulation

```bash
python sim_manager.py run --ETA 25 --INPR 0.05
```

### Sweep INPR (beta-loop study)

```bash
python sim_manager.py beta-loop --ETA 25 --repeats 3
```

Results are appended to `DATA/eta{eta_str}_mu{mu_str}_N{N}/final_corr{1..6}.txt` in the repo root for a beta sweep.

To change the sweeping range, modify the x = np.logspace(-2, -1, 20) in the sim_manager.py beta-loop function.


### Sweep ETA (eta-loop study)

```bash
python sim_manager.py eta-loop --INPR 0.034  --repeats 3
```

Results are appended to `DATA/eta{eta_str}_mu{mu_str}_N{N}/final_corr{1..6}.txt` in the repo root for a eta sweep.

To change the sweeping range, modify the x = np.logspace(0,1.778,20) in the sim_manager.py eta-loop function.

### Sweep ETA 'E' (or BETA 'B') (density study)

```bash
python sim_manager.py density --INPR 0.034 --N 5000 --steps 50 --sweep E 
```

Saves a `.npy` array with per-layer mean densities.

### Single-INPR correlation run (10 repeats)

```bash
python sim_manager.py corr-loop --ETA 25 --INPR 0.05
```

---

## All Options

Every command shares these parameters (shown with defaults):

| Flag | Default | Description |
|------|---------|-------------|
| `--N` | 1000 | Network size (nodes per layer) |
| `--DIM` | 2 | Number of layers |
| `--STIME` | 40000 | Total simulation steps |
| `--ETA` | 25.0 | Noise parameter |
| `--INPR` | 0.05 | Inter-layer coupling |
| `--M` | 1000 | Save every M steps |
| `--GAMMA1` | 100000 | Layer 1 gamma |
| `--GAMMA2` | 100000 | Layer 2 gamma |
| `--DISCTIME` | 20000 | Discard (warm-up) steps |
| `--MU` | 0.05 | Death rate |

Additional per-command flags:

- `beta-loop`: `--repeats N` (default 3)
- `eta-loop`: `--repeats N` (default 3)
- `density`: `--steps N` (default 50)
- `density`: `--sweep STR` (Option E or B) 


