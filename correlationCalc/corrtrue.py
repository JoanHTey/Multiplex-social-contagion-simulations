"""
correlationCalc/corrtrue.py
Reads output.bin and adjacency_matrix.txt (already placed here by
sim_manager.py), computes pair correlations over multiple lags, and
appends the results to the 12 HDF5 accumulators.

Called automatically by sim_manager.py – can also be run standalone from
inside the correlationCalc/ directory.
"""

from pathlib import Path
import numpy as np
import h5py
from joblib import Parallel, delayed

HERE = Path(__file__).parent
ROOT = HERE.parent

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
output_bin = HERE / "output.bin"
adj_txt    = HERE / "adjacency_matrix.txt"
ini_txt    = ROOT / "simulation" / "INITIAL.txt"

if not ini_txt.exists():
    ini_txt = ROOT / "INITIAL.txt"

matrix    = np.fromfile(str(output_bin), dtype=np.int32)
init      = np.loadtxt(str(ini_txt))
NTOT      = int(init[0] * 2)
T         = int(init[2] - init[8])

matrix    = matrix.reshape((T, NTOT + 2), order="C")[:, 1:NTOT + 1]
adj       = np.loadtxt(str(adj_txt), dtype=int)

N          = NTOT // 2
num_lags   = 15
h          = np.linspace(0, 14, num_lags, dtype=int)
num_rand_i = 500
gap        = T // num_rand_i

# ---------------------------------------------------------------------------
# Per-node correlation extraction (parallelised)
# ---------------------------------------------------------------------------

def process_j(j):
    lk = [np.zeros((num_rand_i, num_lags), dtype=int) for _ in range(12)]

    for i in range(num_lags):
        t0 = np.arange(num_rand_i) * gap

        lk[0][:, i]  = matrix[t0,     j]
        lk[1][:, i]  = matrix[t0 + i, j]

        lk[2][:, i]  = matrix[t0,     N + j]
        lk[3][:, i]  = matrix[t0 + i, N + j]

        lk[4][:, i]  = matrix[t0,     j]
        lk[5][:, i]  = matrix[t0 + i, N + j]

        nb1 = adj[j,     1] - 1
        nb2 = adj[N + j, 1] - 1

        lk[6][:, i]  = matrix[t0,     j]
        lk[7][:, i]  = matrix[t0 + i, nb1]

        lk[8][:, i]  = matrix[t0,     j]
        lk[9][:, i]  = matrix[t0 + i, nb2]

        lk[10][:, i] = matrix[t0,     N + j]
        lk[11][:, i] = matrix[t0 + i, nb2]

    return tuple(lk)


if __name__ == "__main__":
    results = Parallel(n_jobs=-1)(delayed(process_j)(j) for j in range(N))
    arrays  = [np.vstack(arr) for arr in zip(*results)]

    filenames = [
        "samelayer1_t0.h5",        "samelayer1_th.h5",
        "samelayer2_t0.h5",        "samelayer2_th.h5",
        "dif_samenode_layer1_t0.h5","dif_samenode_layer2_th.h5",
        "samelayer1_difnode_t0.h5", "samelayer1_difnode_th.h5",
        "difflayer_difnode_t0.h5",  "difflayer_difnode_th.h5",
        "samelayer2_difnode_t0.h5", "samelayer2_difnode_th.h5",
    ]

    for fname, data in zip(filenames, arrays):
        with h5py.File(HERE / fname, "a") as f:
            dset     = f["X"]
            old_size = dset.shape[0]
            dset.resize(old_size + data.shape[0], axis=0)
            dset[old_size:] = data

    print("corrtrue.py: results appended to HDF5 files.")
