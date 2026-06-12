"""
experiment/adjacencymat.py
Generates two-layer REGULAR RANDOM graph adjacency matrix.
(Unlike matrixGeneration/adjacencymat.py which uses Erdős–Rényi.)

Reads INITIAL.txt from the same directory (experiment/).
"""

from pathlib import Path
import sys

import networkx as nx
import numpy as np

HERE     = Path(__file__).parent
INI_PATH = HERE / "INITIAL.txt"

if not INI_PATH.exists():
    print(f"ERROR: Cannot find INITIAL.txt at {INI_PATH}", file=sys.stderr)
    sys.exit(1)

info = np.loadtxt(str(INI_PATH))
N    = int(info[0])
k1   = 30
k2   = 10

print(f"Generating regular random graphs  N={N}  k1={k1}  k2={k2}…")

G1 = nx.random_regular_graph(k1, N)
while not nx.is_connected(G1):
    G1 = nx.random_regular_graph(k1, N)

G2 = nx.random_regular_graph(k2, N)
while not nx.is_connected(G2):
    G2 = nx.random_regular_graph(k2, N)

A1 = nx.to_numpy_array(G1)
A2 = nx.to_numpy_array(G2)

Adja = np.zeros((2 * N, N), dtype=int)

for i in range(N):
    nbrs = np.where(A1[i] == 1)[0]
    Adja[i, 0] = len(nbrs)
    Adja[i, 1:len(nbrs) + 1] = nbrs + 1

for i in range(N):
    nbrs = np.where(A2[i] == 1)[0]
    Adja[N + i, 0] = len(nbrs)
    Adja[N + i, 1:len(nbrs) + 1] = N + nbrs + 1

out_path = HERE / "adjacency_matrix.txt"
np.savetxt(str(out_path), Adja, fmt="%d")
print(f"Adjacency matrix saved to {out_path}")
