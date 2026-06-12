"""
experiment/run_experiment.py
============================
Replaces: Initializer.py, corr.py, image.py

Runs the multi-beta experiment (experiment.f90), computes lag correlations
and densities across beta segments, and plots the results.

Usage
-----
  # Run full pipeline (simulate + analyse + plot)
  python run_experiment.py run

  # Only analyse an existing output.bin
  python run_experiment.py analyse

  # Only plot existing .npy files
  python run_experiment.py plot

Options (shared)
----------------
  --N         int     Nodes per layer              (default 10000)
  --ETA       float   Noise                        (default 1.0)
  --STIME     int     Steps per beta segment       (default 5000)
  --DISCTIME  int     Discard steps per segment    (default 50)
  --LAGS      int     Lags to record per segment   (default 50)
  --RLAGS     int     Lags to compute correlation  (default 3)
  --M         int     Memory reservoir size        (default 1000)
  --GAMMA1    float                                (default 100000)
  --GAMMA2    float                                (default 100000)
  --betas     floats  Space-separated beta values  (default 0.01 0.04 0.07 0.11 0.15)
  --MU        float   Recovery probability         (default 0.5)
  --out-dir   path    Where to save results        (default Data/experiment/)
"""

import argparse
import platform
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
HERE      = Path(__file__).parent
ROOT      = HERE.parent
MATRIX_DIR = ROOT / "matrixGeneration"

EXE_NAME  = "expe.exe" if platform.system() == "Windows" else "expe"
EXPE_EXE  = HERE / EXE_NAME

INITIAL_TXT = HERE / "INITIAL.txt"
BETAS_TXT   = HERE / "betas.txt"
ADJ_TXT     = HERE / "adjacency_matrix.txt"
OUTPUT_BIN  = HERE / "output.bin"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def write_initial(params: list):
    with open(INITIAL_TXT, "w") as f:
        f.write(" ".join(map(str, params)) + "\n")


def write_betas(betas: np.ndarray):
    np.savetxt(str(BETAS_TXT), betas, fmt="%.6f")


def compile_experiment(force: bool = False):
    if EXPE_EXE.exists() and not force:
        print("  [compile] expe executable already exists, skipping.")
        return
    print("  [compile] Compiling experiment.f90…")
    result = subprocess.run(
        ["gfortran", "-O2", "-o", str(EXPE_EXE),
         "experiment.f90", "r1279.f90", "ran2.f"],
        cwd=HERE, capture_output=True, text=True,
    )
    if result.returncode != 0:
        print("  [compile] ERROR:\n", result.stderr)
        sys.exit(1)
    print("  [compile] Done.")


def generate_adjacency_matrix():
    """Use the experiment-local adjacencymat.py (regular random graphs)."""
    print("  [matrix] Generating adjacency matrix…")
    result = subprocess.run(
        [sys.executable, str(HERE / "experiment_adjacencymat.py")],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print("  [matrix] ERROR:\n", result.stderr)
        sys.exit(1)


def run_simulation():
    print("  [sim] Running experiment…")
    result = subprocess.run(
        [str(EXPE_EXE)], cwd=HERE, capture_output=True, text=True,
    )
    if result.returncode != 0:
        print("  [sim] ERROR:\n", result.stderr)
        sys.exit(1)
    print("  [sim] Done.")


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

def analyse(args, out_dir: Path):
    """Read output.bin and compute correlations + densities per beta segment."""
    init   = np.loadtxt(str(INITIAL_TXT))
    betas  = np.loadtxt(str(BETAS_TXT))
    n_beta = len(betas)

    NTOT   = int(init[0] * 2)
    N      = NTOT // 2
    lags   = int(init[9])
    rlags  = args.rlags

    adja   = np.loadtxt(str(ADJ_TXT), dtype=np.int32)

    matrix = np.fromfile(str(OUTPUT_BIN), dtype=np.int32)
    wrap   = lags * 100
    Ts     = 10 * n_beta        # time-series segments
    T      = wrap * n_beta

    matrix = matrix.reshape((T, NTOT + 3), order="C")[:, 1:NTOT + 2]

    # Fix time index wrapping across beta segments
    max_steps = 100 * int(init[8])
    for i in range(n_beta - 1):
        matrix[wrap * i:wrap * (i + 1), 0] += max_steps * i

    lags_long = lags * 10

    shape = (Ts, rlags)
    corr  = [np.zeros(shape) for _ in range(6)]
    ecorr = [np.zeros(shape) for _ in range(6)]
    den1  = np.zeros(Ts)
    den2  = np.zeros(Ts)

    for i in range(Ts):
        sl  = slice(i * lags_long, (i + 1) * lags_long)
        for j in range(rlags):
            sl0 = slice(i * lags_long, (i + 1) * lags_long - j)
            slj = slice(i * lags_long + j, (i + 1) * lags_long)

            def _corr(a, b):
                num = np.mean(a * b) - np.mean(a) * np.mean(b)
                den = np.std(a) * np.std(b)
                return num / den if den > 0 else 0.0

            def _ecorr(a, b):
                num = np.mean(a * b, axis=1) - np.mean(a, axis=1) * np.mean(b, axis=1)
                den = np.std(a, axis=1) * np.std(b, axis=1) + 1e-7
                return np.std(num / den)

            L1   = matrix[sl0, 1:N + 1]
            L1j  = matrix[slj, 1:N + 1]
            L2   = matrix[sl0, N + 1:NTOT + 1]
            L2j  = matrix[slj, N + 1:NTOT + 1]
            nb1  = matrix[slj, adja[:N,  1]]
            nb2  = matrix[slj, adja[N:NTOT, 1]]

            corr[0][i, j]  = _corr(L1,  L1j)
            corr[1][i, j]  = _corr(L2,  L2j)
            corr[2][i, j]  = _corr(L1,  L2j)
            corr[3][i, j]  = _corr(L1,  nb1)
            corr[4][i, j]  = _corr(L2,  nb2)
            corr[5][i, j]  = _corr(matrix[sl0, 1:N + 1], nb2)

            ecorr[0][i, j] = _ecorr(matrix[sl0, 1:N+1],      matrix[slj, 1:N+1])
            ecorr[1][i, j] = _ecorr(matrix[sl0, N+1:NTOT+1], matrix[slj, N+1:NTOT+1])
            ecorr[2][i, j] = _ecorr(matrix[sl0, 1:N+1],      matrix[slj, N+1:NTOT+1])
            ecorr[3][i, j] = _ecorr(matrix[sl0, 1:N+1],      matrix[slj, adja[:N, 1]])
            ecorr[4][i, j] = _ecorr(matrix[sl0, N+1:NTOT+1], matrix[slj, adja[N:NTOT, 1]])
            ecorr[5][i, j] = _ecorr(matrix[sl0, 1:N+1],      matrix[slj, adja[N:NTOT, 1]])

        den1[i] = np.mean(matrix[sl, 1:N + 1])
        den2[i] = np.mean(matrix[sl, N + 1:NTOT + 1])

    out_dir.mkdir(parents=True, exist_ok=True)
    for k in range(6):
        np.save(str(out_dir / f"corr{k+1}.npy"),  corr[k])
        np.save(str(out_dir / f"ecorr{k+1}.npy"), ecorr[k])
    np.save(str(out_dir / "den1.npy"), den1)
    np.save(str(out_dir / "den2.npy"), den2)
    print(f"  [analyse] Results saved to {out_dir}")


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot(args, out_dir: Path):
    betas = np.loadtxt(str(BETAS_TXT))
    n_beta = len(betas)
    lags   = int(np.loadtxt(str(INITIAL_TXT))[9])
    Ts     = 10 * n_beta
    t      = np.linspace(0, (Ts - 1) * lags * 100 / 10, Ts)

    c  = [np.load(str(out_dir / f"corr{k+1}.npy"))  for k in range(6)]
    ec = [np.load(str(out_dir / f"ecorr{k+1}.npy")) for k in range(6)]
    d1 = np.load(str(out_dir / "den1.npy"))
    d2 = np.load(str(out_dir / "den2.npy"))

    labels = [
        r"$\rho_{uu}^{11}(1)$", r"$\rho_{uu}^{22}(1)$",
        r"$\rho_{uu}^{12}(1)$", r"$\rho_{uv}^{11}(1)$",
        r"$\rho_{uv}^{22}(1)$", r"$\rho_{uv}^{12}(1)$",
    ]
    colors = ["tab:blue", "tab:orange", "tab:green",
              "tab:red",  "tab:purple", "tab:brown"]
    vline_colors = ["red", "green", "purple", "blue", "gold",
                    "cyan", "magenta", "black"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    for k in range(6):
        ax1.plot(t, c[k][:, 1], label=labels[k],
                 color=colors[k], marker="o", markersize=4, linewidth=1.5)
        ax1.fill_between(
            t,
            c[k][:, 1] - ec[k][:, 1],
            c[k][:, 1] + ec[k][:, 1],
            color=colors[k], alpha=0.15,
        )

    # Vertical lines at beta transitions
    segment_len = Ts // n_beta
    for bi, (beta, vc) in enumerate(zip(betas, vline_colors)):
        tx = t[bi * segment_len]
        ax1.axvline(tx, linestyle="--", color=vc, linewidth=1,
                    alpha=0.7, label=fr"$\beta={beta:.2f}$")
        ax2.axvline(tx, linestyle="--", color=vc, linewidth=1, alpha=0.7,
                    label=fr"$\beta={beta:.2f}$")

    ax1.set_xlabel("Time steps", fontsize=14)
    ax1.set_ylabel(r"$\rho$", fontsize=14)
    ax1.legend(fontsize=9, loc="best")
    ax1.set_ylim(0, 1.05)

    ax2.plot(t, d1, label=r"$\langle X \rangle_1$",
             color="tab:blue", marker="o", markersize=4, linewidth=1.5)
    ax2.plot(t, d2, label=r"$\langle X \rangle_2$",
             color="tab:orange", marker="o", markersize=4, linewidth=1.5)
    ax2.set_xlabel("Time steps", fontsize=14)
    ax2.set_ylabel(r"$\langle X \rangle$", fontsize=14)
    ax2.legend(fontsize=10, loc="best")

    plt.tight_layout()
    fig_path = out_dir / "expe.pdf"
    plt.savefig(str(fig_path))
    print(f"  [plot] Figure saved to {fig_path}")
    plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser():
    p = argparse.ArgumentParser(
        description="Experiment runner (multi-beta simulation)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    sub = p.add_subparsers(dest="command", required=True)

    shared = argparse.ArgumentParser(add_help=False)
    shared.add_argument("--N",        type=int,   default=10000)
    shared.add_argument("--DIM",      type=int,   default=2)
    shared.add_argument("--ETA",      type=float, default=1.0)
    shared.add_argument("--MU",       type=float, default=0.5)
    shared.add_argument("--STIME",    type=int,   default=5000)
    shared.add_argument("--DISCTIME", type=int,   default=50)
    shared.add_argument("--LAGS",     type=int,   default=50)
    shared.add_argument("--M",        type=int,   default=1000)
    shared.add_argument("--GAMMA1",   type=float, default=100000.0)
    shared.add_argument("--GAMMA2",   type=float, default=100000.0)
    shared.add_argument("--betas",    type=float, nargs="+",
                        default=[0.01, 0.04, 0.07, 0.11, 0.15])
    shared.add_argument("--rlags",    type=int,   default=3,
                        help="Number of lags to compute correlations for")
    shared.add_argument("--out-dir",  type=Path,
                        default=None,
                        help="Output directory (default: Data/experiment/eta{ETA}_mu{MU}_N{N})")

    for cmd in ("run", "analyse", "plot"):
        sub.add_parser(cmd, parents=[shared],
                       help={"run": "Simulate + analyse + plot",
                             "analyse": "Analyse existing output.bin",
                             "plot": "Plot existing .npy files"}[cmd])

    return p


def resolve_out_dir(args) -> Path:
    if args.out_dir:
        return Path(args.out_dir)
    eta_str = str(args.ETA).replace(".", "p")
    mu_str  = str(args.MU).replace(".", "p")
    return ROOT / "Data" / "experiment" / f"eta{eta_str}_mu{mu_str}_N{args.N}"


if __name__ == "__main__":
    parser = build_parser()
    args   = parser.parse_args()
    out_dir = resolve_out_dir(args)

    if args.command == "run":
        params = [args.N, args.DIM, args.STIME, args.ETA, 0,
                  args.M, args.GAMMA1, args.GAMMA2, args.DISCTIME, args.LAGS]
        write_initial(params)
        write_betas(np.array(args.betas))
        compile_experiment()
        generate_adjacency_matrix()
        run_simulation()
        analyse(args, out_dir)
        plot(args, out_dir)

    elif args.command == "analyse":
        analyse(args, out_dir)

    elif args.command == "plot":
        plot(args, out_dir)
