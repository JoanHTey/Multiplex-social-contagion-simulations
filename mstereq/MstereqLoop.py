"""
theory/MstereqLoop.py
=====================
Runs the mean-field steady-state equation solver (mstereq.f90 / meq.exe)
over a sweep of INPR (beta) values and saves per-layer correlation
approximations and densities.

Usage
-----
  python MstereqLoop.py [options]

Options
-------
  --N        int    Network size per layer        (default 5000)
  --ETA      float  Noise parameter               (default 0.01)
  --MU       float  Recovery probability          (default 0.5)
  --GAMMA1   float  Layer-1 gamma                 (default 100000)
  --GAMMA2   float  Layer-2 gamma                 (default 100000)
  --STIME    int    Simulation steps              (default 20000)
  --M        int    Memory size                   (default 1000)
  --DIM      int    Number of layers              (default 2)

Output
------
  Data/eta{ETA}_mu{MU}_N{N}/
      final_corr1msteq.txt
      final_corr2msteq.txt
      den1msteq.txt
      den2msteq.txt
"""

import argparse
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
HERE     = Path(__file__).parent
ROOT     = HERE.parent
MATRIX_DIR = ROOT / "matrixGeneration"
SIM_DIR    = ROOT / "simulation"

EXE_NAME = "meq.exe" if platform.system() == "Windows" else "meq"
MEQ_EXE  = HERE / EXE_NAME

INITIAL_TXT     = ROOT / "INITIAL.txt"
INITIAL_IN_SIM  = SIM_DIR / "INITIAL.txt"
ADJ_TXT         = ROOT / "adjacency_matrix.txt"
OUTPUT2_TXT     = HERE / "output2INPR.txt"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def write_initial(params: list, path: Path = INITIAL_TXT):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write(" ".join(map(str, params)) + "\n")


def compile_mstereq(force: bool = False):
    """Compile mstereq.f90 if the executable is missing."""
    if MEQ_EXE.exists() and not force:
        print("  [compile] meq executable already exists, skipping.")
        return
    print("  [compile] Compiling mstereq.f90…")
    result = subprocess.run(
        ["gfortran", "-O2", "-o", str(MEQ_EXE), "mstereq.f90"],
        cwd=HERE, capture_output=True, text=True,
    )
    if result.returncode != 0:
        print("  [compile] ERROR:\n", result.stderr)
        sys.exit(1)
    print("  [compile] Done.")


def generate_adjacency_matrix():
    """Generate adjacency_matrix.txt via matrixGeneration/adjacencymat.py."""
    print("  [matrix] Generating adjacency matrix…")
    result = subprocess.run(
        [sys.executable, str(MATRIX_DIR / "adjacencymat.py")],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print("  [matrix] ERROR:\n", result.stderr)
        sys.exit(1)


def run_mstereq():
    """Move inputs into theory dir and run meq executable."""
    # mstereq.f90 reads INITIAL.txt and adjacency_matrix.txt from its own dir
    shutil.copy(str(INITIAL_TXT),  str(HERE / "INITIAL.txt"))
    # also keep a copy in simulation/ so adjacencymat.py can read N from there
    shutil.copy(str(INITIAL_TXT),  str(INITIAL_IN_SIM))
    shutil.copy(str(ADJ_TXT),      str(HERE / "adjacency_matrix.txt"))

    result = subprocess.run(
        [str(MEQ_EXE)], cwd=HERE, capture_output=True, text=True,
    )
    if result.returncode != 0:
        print("  [meq] ERROR:\n", result.stderr)
        sys.exit(1)


def personcorr_approx(filename: Path, h: np.ndarray, layer: int,
                      mu: float, N: int) -> np.ndarray:
    """
    Approximate Pearson correlation from mstereq output.

    Parameters
    ----------
    filename : path to output2INPR.txt
    h        : array of lag values (integers)
    layer    : 1 or 2
    mu       : recovery probability
    N        : nodes per layer
    """
    mat = np.loadtxt(str(filename), dtype=float)
    b = layer
    x_mean  = np.mean(mat[0, N * (b - 1):N * b])   # mean steady-state density
    q_mean  = np.mean(mat[1, N * (b - 1):N * b])   # mean Q factor

    term1 = (1 - q_mean) + mu * q_mean * (q_mean * (1 - mu)) ** h
    term2 = 1 - q_mean * (1 - mu)
    corr  = ((term1 / term2) * x_mean - x_mean ** 2) / (x_mean - x_mean ** 2)
    return corr


# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------

def run_beta_sweep(args):
    N      = args.N
    mu     = args.MU
    h      = np.linspace(0, 14, 15, dtype=int)

    eta_str  = str(args.ETA).replace(".", "p")
    mu_str   = str(args.MU).replace(".", "p")
    out_dir  = ROOT / "Data" / f"eta{eta_str}_mu{mu_str}_N{N}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"  [output] Results will be saved to {out_dir}")

    compile_mstereq()

    inpr_values = np.logspace(-2, -0.3, 60)
    t0 = time.time()
    print(f"\n=== Mstereq beta sweep | {len(inpr_values)} INPR values ===")

    for idx, INPR in enumerate(inpr_values, 1):
        print(f"  INPR={INPR:.6f}  ({idx}/{len(inpr_values)})")
        params = [N, args.DIM, args.STIME, args.ETA, INPR,
                  args.M, args.GAMMA1, args.GAMMA2,args.MU]
        write_initial(params)
        generate_adjacency_matrix()
        run_mstereq()

        mat = np.loadtxt(str(OUTPUT2_TXT), dtype=float)

        # Correlation approximations
        c1 = personcorr_approx(OUTPUT2_TXT, h, 1, mu, N)
        c2 = personcorr_approx(OUTPUT2_TXT, h, 2, mu, N)

        with open(out_dir / "final_corr1msteq.txt", "a") as f:
            np.savetxt(f, [np.concatenate([[INPR], c1])], fmt="%.6f")
        with open(out_dir / "final_corr2msteq.txt", "a") as f:
            np.savetxt(f, [np.concatenate([[INPR], c2])], fmt="%.6f")

        # Densities (first row of output = X values)
        den1 = np.mean(mat[0, N * 0:N * 1])
        den2 = np.mean(mat[0, N * 1:N * 2])

        with open(out_dir / "den1msteq.txt", "a") as f:
            np.savetxt(f, [[INPR, den1]], fmt="%.6f")
        with open(out_dir / "den2msteq.txt", "a") as f:
            np.savetxt(f, [[INPR, den2]], fmt="%.6f")

    print(f"\nMstereq sweep finished in {time.time() - t0:.1f}s")
    print(f"Results in {out_dir}")

def run_eta_sweep(args):
    N      = args.N
    mu     = args.MU
    h      = np.linspace(0, 14, 15, dtype=int)

    beta_str  = str(args.INPR).replace(".", "p")
    mu_str   = str(args.MU).replace(".", "p")
    out_dir  = ROOT / "Data" / f"beta{beta_str}_mu{mu_str}_N{N}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"  [output] Results will be saved to {out_dir}")

    compile_mstereq()

    eta_values = np.logspace(0, 2, 60)
    t0 = time.time()
    print(f"\n=== Mstereq eta sweep | {len(eta_values)} eta values ===")

    for idx, ETA in enumerate(eta_values, 1):
        print(f"  ETA={ETA:.6f}  ({idx}/{len(eta_values)})")
        params = [N, args.DIM, args.STIME, ETA, args.INPR,
                  args.M, args.GAMMA1, args.GAMMA2, args.MU]
        write_initial(params)
        generate_adjacency_matrix()
        run_mstereq()

        mat = np.loadtxt(str(OUTPUT2_TXT), dtype=float)

        # Correlation approximations
        c1 = personcorr_approx(OUTPUT2_TXT, h, 1, mu, N)
        c2 = personcorr_approx(OUTPUT2_TXT, h, 2, mu, N)

        with open(out_dir / "final_corr1msteq.txt", "a") as f:
            np.savetxt(f, [np.concatenate([[ETA], c1])], fmt="%.6f")
        with open(out_dir / "final_corr2msteq.txt", "a") as f:
            np.savetxt(f, [np.concatenate([[ETA], c2])], fmt="%.6f")

        # Densities (first row of output = X values)
        den1 = np.mean(mat[0, N * 0:N * 1])
        den2 = np.mean(mat[0, N * 1:N * 2])

        with open(out_dir / "den1msteq.txt", "a") as f:
            np.savetxt(f, [[ETA, den1]], fmt="%.6f")
        with open(out_dir / "den2msteq.txt", "a") as f:
            np.savetxt(f, [[ETA, den2]], fmt="%.6f")

    print(f"\nMstereq sweep finished in {time.time() - t0:.1f}s")
    print(f"Results in {out_dir}")



# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser():
    parser = argparse.ArgumentParser(
        description="Mstereq mean-field beta sweep",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    sub = parser.add_subparsers(dest="command", required=True)

    # Shared parameter defaults
    shared = argparse.ArgumentParser(add_help=False)
    shared.add_argument("--N",      type=int,   default=5000)
    shared.add_argument("--DIM",    type=int,   default=2)
    shared.add_argument("--STIME",  type=int,   default=20000)
    shared.add_argument("--ETA",    type=float, default=0.01)
    shared.add_argument("--INPR",    type=float, default=0.034)
    shared.add_argument("--MU",     type=float, default=0.5)
    shared.add_argument("--M",      type=int,   default=1000)
    shared.add_argument("--GAMMA1", type=float, default=100000.0)
    shared.add_argument("--GAMMA2", type=float, default=100000.0)
    
    p_run_eta = sub.add_parser("eta", parents=[shared], help="ETA sweep")
    p_run_eta.set_defaults(func=run_eta_sweep)

    p_run_beta = sub.add_parser("run", parents=[shared], help="Beta sweep")
    p_run_beta.set_defaults(func=run_beta_sweep)
    
    return parser


if __name__ == "__main__":
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)
