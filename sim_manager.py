"""
sim_manager.py  –  Unified simulation manager for the TFM project
===================================================================
Replaces: loop.py, corrLoop.py, FinalBetaLoop.py, density.py,
          scripts/sync_script_corr.bat, scripts/betaLoop.bat,
          scripts/den_script.bat

Usage
-----
  python sim_manager.py <command> [options]

Commands
--------
  run          Run one simulation with given parameters
  beta-loop    Sweep INPR values (correlation vs beta study)
  density      Sweep ETA values (density study)
  corr-loop    Single INPR correlation run (10 repeats, collects h5 data)

Examples
--------
  python sim_manager.py run --N 1000 --ETA 25 --INPR 0.05
  python sim_manager.py beta-loop --ETA 25 --repeats 3
  python sim_manager.py density --INPR 0.034 --N 5000
  python sim_manager.py corr-loop --ETA 25 --INPR 0.05
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
# Paths  (all relative to the repo root – adjust if needed)
# ---------------------------------------------------------------------------
ROOT          = Path(__file__).parent
SIM_DIR       = ROOT / "simulation"
MATRIX_DIR    = ROOT / "matrixGeneration"
CORR_DIR      = ROOT / "correlationCalc"
SCRIPTS_DIR   = ROOT / "scripts"

INITIAL_TXT   = ROOT / "INITIAL.txt"
ADJ_TXT       = ROOT / "adjacency_matrix.txt"
OUTPUT_BIN    = SIM_DIR / "output.bin"
ADJ_IN_SIM    = SIM_DIR / "adjacency_matrix.txt"
INITIAL_IN_SIM= SIM_DIR / "INITIAL.txt"
SIM_EXE       = SIM_DIR / ("simulation.exe" if platform.system() == "Windows" else "simulation")

CORR_FILES = [
    "samelayer1_t0.h5",   "samelayer1_th.h5",
    "samelayer2_t0.h5",   "samelayer2_th.h5",
    "dif_samenode_layer1_t0.h5", "dif_samenode_layer2_th.h5",
    "samelayer1_difnode_t0.h5",  "samelayer1_difnode_th.h5",
    "difflayer_difnode_t0.h5",   "difflayer_difnode_th.h5",
    "samelayer2_difnode_t0.h5",  "samelayer2_difnode_th.h5",
]
FINAL_CORR_FILES = [f"final_corr{i}.txt" for i in range(1, 7)]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def write_initial(params: list, path: Path = INITIAL_TXT):
    """Write simulation parameters to an INITIAL.txt file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write(" ".join(map(str, params)) + "\n")
    print(f"  [params] Written to {path}")


def compile_fortran(force: bool = False):
    """Compile the Fortran simulation if the executable doesn't exist."""
    if SIM_EXE.exists() and not force:
        print("  [compile] Executable already exists, skipping compilation.")
        return
    print("  [compile] Compiling Fortran simulation…")
    result = subprocess.run(
        ["gfortran", "-O2", "-o", str(SIM_EXE), "sync.f90", "r1279.f90", "ran2.f"],
        cwd=SIM_DIR,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print("  [compile] ERROR:\n", result.stderr)
        sys.exit(1)
    print("  [compile] Compilation successful.")


def generate_adjacency_matrix():
    """Run adjacencymat.py to produce adjacency_matrix.txt in ROOT."""
    print("  [matrix] Generating adjacency matrix…")
    result = subprocess.run(
        [sys.executable, str(MATRIX_DIR / "adjacencymat.py")],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print("  [matrix] ERROR:\n", result.stderr)
        sys.exit(1)
    print("  [matrix] Done.")


def run_fortran_simulation():
    """Move input files, run the executable, move output back."""
    # Move inputs into sim dir
    shutil.move(str(INITIAL_TXT),  str(INITIAL_IN_SIM))
    shutil.move(str(ADJ_TXT),      str(ADJ_IN_SIM))

    print("  [sim] Running Fortran simulation…")
    result = subprocess.run(
        [str(SIM_EXE)], cwd=SIM_DIR, capture_output=True, text=True
    )
    if result.returncode != 0:
        print("  [sim] ERROR:\n", result.stderr)
        sys.exit(1)
    print("  [sim] Simulation finished.")


def move_output_to_corr():
    """Move output.bin and adjacency_matrix from sim dir to correlationCalc."""
    CORR_DIR.mkdir(parents=True, exist_ok=True)
    shutil.move(str(OUTPUT_BIN),  str(CORR_DIR / "output.bin"))
    shutil.move(str(ADJ_IN_SIM),  str(CORR_DIR / "adjacency_matrix.txt"))
    print("  [output] Moved output to correlationCalc/")


def run_corrtrue():
    """Run correlationCalc/corrtrue.py."""
    print("  [corr] Running corrtrue.py…")
    result = subprocess.run(
        [sys.executable, str(CORR_DIR / "corrtrue.py")],
        cwd=CORR_DIR, capture_output=True, text=True,
    )
    if result.returncode != 0:
        print("  [corr] ERROR:\n", result.stderr)
        sys.exit(1)
    print("  [corr] Correlation calculation done.")


def init_h5_files(N: int):
    """Create (or reset) the 12 HDF5 accumulator files in correlationCalc."""
    try:
        import h5py
    except ImportError:
        print("  [h5] h5py not installed – run: pip install h5py")
        sys.exit(1)

    CORR_DIR.mkdir(parents=True, exist_ok=True)
    for fname in CORR_FILES:
        fpath = CORR_DIR / fname
        if fpath.exists():
            fpath.unlink()
        with h5py.File(fpath, "w") as f:
            f.create_dataset(
                "X", shape=(0, 15), maxshape=(None, 15),
                dtype="int32", chunks=(N * 500, 15),
            )
    print(f"  [h5] Initialised {len(CORR_FILES)} HDF5 files.")


def corrcalc(X1, X2) -> np.ndarray:
    """Pearson correlation between two arrays."""
    mu1, mu2 = np.mean(X1, axis=0), np.mean(X2, axis=0)
    cov  = np.mean(X1 * X2, axis=0) - mu1 * mu2
    denom = np.std(X1, axis=0) * np.std(X2, axis=0)
    return np.where(denom == 0, 0.0, cov / denom)


def full_pipeline(params: list, with_corr: bool = True):
    """Generate matrix → compile → simulate → (optionally) calc correlation."""
    write_initial(params)
    generate_adjacency_matrix()
    compile_fortran()
    run_fortran_simulation()
    if with_corr:
        move_output_to_corr()
        run_corrtrue()


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def cmd_run(args):
    """Run a single simulation (no correlation post-processing)."""
    params = [
        args.N, args.DIM, args.STIME, args.ETA,
        args.INPR, args.M, args.GAMMA1, args.GAMMA2, args.DISCTIME,args.MU
    ]
    print(f"\n=== Single simulation | ETA={args.ETA} INPR={args.INPR} ===")
    full_pipeline(params, with_corr=False)
    print("Done.\n")


def cmd_beta_loop(args):
    """
    Sweep INPR values: for each value run corrLoop pipeline and accumulate
    results into final_corr{1..6}.txt.
    """
    try:
        import h5py
    except ImportError:
        print("h5py required. Install with: pip install h5py")
        sys.exit(1)

    x = np.logspace(-2, -1, 20)
    log_mids = 0.5 * (np.log10(x[1:]) + np.log10(x[:-1]))
    inpr_values = np.sort(np.append(x, 10 ** log_mids))

    t0 = time.time()
    print(f"\n=== Beta loop | {len(inpr_values)} INPR values × {args.repeats} repeats ===")

    compile_fortran()
    N = args.N

    # Build output directory name from ETA, MU, N
    eta_str = str(args.ETA).replace(".", "p")
    mu_str  = str(args.MU).replace(".", "p")
    out_dir = ROOT / "Data" / f"eta{eta_str}_mu{mu_str}_N{N}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"  [output] Results will be saved to {out_dir}")

    for rep in range(args.repeats):
        print(f"\n--- Repeat {rep + 1}/{args.repeats} ---")
        for idx, INPR in enumerate(inpr_values, 1):
            print(f"  INPR={INPR:.6f}  ({idx}/{len(inpr_values)})")
            params = [
                N, args.DIM, args.STIME, args.ETA,
                INPR, args.M, args.GAMMA1, args.GAMMA2, args.DISCTIME, args.MU
            ]


            # Init h5 accumulators
            init_h5_files(N)

            # Run simulation 10 times accumulating data
            for run_i in range(10):
                write_initial(params)
                generate_adjacency_matrix()
                run_fortran_simulation()
                move_output_to_corr()
                run_corrtrue()

            # Read accumulated h5 data and compute correlations
            data = {}
            for i, fname in enumerate(CORR_FILES, 1):
                with h5py.File(CORR_DIR / fname, "r") as f:
                    data[f"X{i}"] = f["X"][:]

            pairs = [(1,2), (3,4), (5,6), (7,8), (9,10), (11,12)]
            size = data["X1"].shape
            for ci, (a, b) in enumerate(pairs, 1):
                corr = np.array([
                    corrcalc(data[f"X{a}"][:, i], data[f"X{b}"][:, i])
                    for i in range(size[1])
                ])
                with open(out_dir / f"final_corr{ci}.txt", "a") as f:
                    np.savetxt(f, [np.concatenate([[INPR], corr])], fmt="%.6f")

    print(f"\nEta loop finished in {time.time() - t0:.1f}s")

def cmd_eta_loop(args):
    """
    Sweep ETA values: for each value run corrLoop pipeline and accumulate
    results into final_corr{1..6}.txt.
    """
    try:
        import h5py
    except ImportError:
        print("h5py required. Install with: pip install h5py")
        sys.exit(1)

    x = np.logspace(0,1.778,20) #BETA0034 The maximum eta needs to satisfy 1 GE ETA*INPR*MU
    log_mids = 0.5 * (np.log10(x[1:]) + np.log10(x[:-1]))
    eta_values = np.sort(np.append(x, 10 ** log_mids))

    t0 = time.time()
    print(f"\n=== Eta loop | {len(eta_values)} ETA values × {args.repeats} repeats ===")

    compile_fortran()
    N = args.N

    # Build output directory name from BETA, MU, N
    beta_str = str(args.INPR).replace(".", "p")
    mu_str  = str(args.MU).replace(".", "p")
    out_dir = ROOT / "Data" / f"beta{beta_str}_mu{mu_str}_N{N}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"  [output] Results will be saved to {out_dir}")


    for rep in range(args.repeats):
        print(f"\n--- Repeat {rep + 1}/{args.repeats} ---")
        for idx, ETA in enumerate(eta_values, 1):
            print(f"  ETA={ETA:.6f}  ({idx}/{len(eta_values)})")
            params = [
                N, args.DIM, args.STIME, ETA,
                args.INPR, args.M, args.GAMMA1, args.GAMMA2, args.DISCTIME,args.MU
            ]


            # Init h5 accumulators
            init_h5_files(N)

            # Run simulation 10 times accumulating data
            for run_i in range(10):
                write_initial(params)
                generate_adjacency_matrix()
                run_fortran_simulation()
                move_output_to_corr()
                run_corrtrue()

            # Read accumulated h5 data and compute correlations
            data = {}
            for i, fname in enumerate(CORR_FILES, 1):
                with h5py.File(CORR_DIR / fname, "r") as f:
                    data[f"X{i}"] = f["X"][:]

            pairs = [(1,2), (3,4), (5,6), (7,8), (9,10), (11,12)]
            size = data["X1"].shape
            for ci, (a, b) in enumerate(pairs, 1):
                corr = np.array([
                    corrcalc(data[f"X{a}"][:, i], data[f"X{b}"][:, i])
                    for i in range(size[1])
                ])
                with open(out_dir / f"final_corr{ci}.txt", "a") as f:
                    np.savetxt(f, [np.concatenate([[ETA], corr])], fmt="%.6f")

    print(f"\nEta loop finished in {time.time() - t0:.1f}s")


def cmd_density(args):

    """Sweep ETA or BETA values and compute mean layer densities."""
    compile_fortran()

    if args.sweep == 'E':

        eta_values = np.logspace(0, 1.778, args.steps)
        N, T = args.N, args.STIME - args.DISCTIME
        densities = np.zeros((len(eta_values), 2))

        print(f"\n=== Density sweep | {len(eta_values)} ETA values ===")

        # Build output directory name from BETA, MU, N
        beta_str = str(args.INPR).replace(".", "p")
        mu_str  = str(args.MU).replace(".", "p")
        out_dir = ROOT / "Data" / f"beta{beta_str}_mu{mu_str}_N{N}"
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"  [output] Results will be saved to {out_dir}")

        for i, ETA in enumerate(eta_values):
            print(f"  ETA={ETA:.4f}  ({i+1}/{len(eta_values)})")
            params = [
                N, args.DIM, args.STIME, ETA,
                args.INPR, args.M, args.GAMMA1, args.GAMMA2, args.DISCTIME,args.MU
            ]
            write_initial(params)
            generate_adjacency_matrix()
            run_fortran_simulation()

            matrix = np.fromfile(str(OUTPUT_BIN), dtype=np.int32)
            matrix = matrix.reshape((T, N * 2 + 2), order="C")[:, 1:N * 2 + 1]
            densities[i, 0] = np.mean(matrix[:, :N])
            densities[i, 1] = np.mean(matrix[:, N:])

        out_path = out_dir / f"density_beta{str(args.INPR).replace('.', '')}_N{N}.npy"
        np.save(str(out_path), densities)
        print(f"\nDensity results saved to {out_path}")
    
    elif args.sweep == "B":

        inpr_values = np.logspace(-2, -1, args.steps)
        N, T = args.N, args.STIME - args.DISCTIME
        densities = np.zeros((len(inpr_values), 2))

        print(f"\n=== Density sweep | {len(inpr_values)} BETA values ===")

        # Build output directory name from ETA, MU, N
        eta_str = str(args.ETA).replace(".", "p")
        mu_str  = str(args.MU).replace(".", "p")
        out_dir = ROOT / "Data" / f"eta{eta_str}_mu{mu_str}_N{N}"
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"  [output] Results will be saved to {out_dir}")

        for i, INPR in enumerate(inpr_values):
            print(f"  INPR={INPR:.4f}  ({i+1}/{len(inpr_values)})")
            params = [
                N, args.DIM, args.STIME, args.ETA,
                INPR, args.M, args.GAMMA1, args.GAMMA2, args.DISCTIME,args.MU
            ]
            write_initial(params)
            generate_adjacency_matrix()
            run_fortran_simulation()

            matrix = np.fromfile(str(OUTPUT_BIN), dtype=np.int32)
            matrix = matrix.reshape((T, N * 2 + 2), order="C")[:, 1:N * 2 + 1]
            densities[i, 0] = np.mean(matrix[:, :N])
            densities[i, 1] = np.mean(matrix[:, N:])

        out_path = out_dir / f"density_eta{str(args.ETA).replace('.', '')}_N{N}.npy"
        np.save(str(out_path), densities)
        print(f"\nDensity results saved to {out_path}")
    else:
        print("--sweep only has two possible values, B for Beta sweeps and E for Eta sweeps.")


def cmd_corr_loop(args):
    """Single-INPR correlation run (10 repeats, matches old corrLoop.py)."""
    try:
        import h5py
    except ImportError:
        print("h5py required. Install with: pip install h5py")
        sys.exit(1)

    params = [
        args.N, args.DIM, args.STIME, args.ETA,
        args.INPR, args.M, args.GAMMA1, args.GAMMA2, args.DISCTIME,args.MU
    ]
    print(f"\n=== Corr loop | ETA={args.ETA} INPR={args.INPR} ===")
    compile_fortran()
    init_h5_files(N)
    write_initial(params)

    for run_i in range(10):
        print(f"  Run {run_i + 1}/10")
        generate_adjacency_matrix()
        run_fortran_simulation()
        move_output_to_corr()
        run_corrtrue()

    print("Corr loop done.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser():
    parser = argparse.ArgumentParser(
        description="TFM simulation manager",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # Shared parameter defaults
    shared = argparse.ArgumentParser(add_help=False)
    shared.add_argument("--N",        type=int,   default=1000)
    shared.add_argument("--DIM",      type=int,   default=2)
    shared.add_argument("--STIME",    type=int,   default=40000)
    shared.add_argument("--ETA",      type=float, default=25.0)
    shared.add_argument("--INPR",     type=float, default=0.034)
    shared.add_argument("--M",        type=int,   default=1000)
    shared.add_argument("--GAMMA1",   type=float, default=100000.0)
    shared.add_argument("--GAMMA2",   type=float, default=100000.0)
    shared.add_argument("--DISCTIME", type=int,   default=20000)
    shared.add_argument("--MU", type=float, default=0.5)

    # run
    p_run = sub.add_parser("run", parents=[shared], help="Single simulation run")
    p_run.set_defaults(func=cmd_run)

    # beta-loop
    p_beta = sub.add_parser("beta-loop", parents=[shared], help="Sweep INPR values")
    p_beta.add_argument("--repeats", type=int, default=3)
    p_beta.set_defaults(func=cmd_beta_loop)

    # eta-loop
    p_eta = sub.add_parser("eta-loop", parents=[shared], help="Sweep ETA values")
    p_eta.add_argument("--repeats", type=int, default=3)
    p_eta.set_defaults(func=cmd_eta_loop)

    # density
    p_den = sub.add_parser("density", parents=[shared], help="Sweep ETA or BETA for density")
    p_den.add_argument("--steps", type=int, default=50, help="Number of steps")
    p_den.add_argument("--sweep", type=str, default="E", help="Decide if its a beta (B) or an eta (E) sweep")
    p_den.set_defaults(func=cmd_density)

    # corr-loop
    p_corr = sub.add_parser("corr-loop", parents=[shared], help="10-repeat correlation run")
    p_corr.set_defaults(func=cmd_corr_loop)

    return parser


if __name__ == "__main__":
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)
