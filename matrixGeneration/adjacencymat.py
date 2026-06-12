"""
matrixGeneration/adjacencymat.py
Generates the two-layer Erdős–Rényi adjacency matrix and writes it to
adjacency_matrix.txt in the repo root.

Reads INITIAL.txt from simulation/ (placed there by sim_manager.py before
this script is called).
"""

import os
import sys
from pathlib import Path

import networkx as nx
import numpy as np

# ---------------------------------------------------------------------------
# Locate INITIAL.txt robustly (works whether called directly or via manager)
# ---------------------------------------------------------------------------
THIS_DIR = Path(__file__).parent
ROOT     = THIS_DIR.parent
INI_PATH = ROOT / "simulation" / "INITIAL.txt"

if not INI_PATH.exists():
    # Fall back: look in ROOT
    INI_PATH = ROOT / "INITIAL.txt"

if not INI_PATH.exists():
    print(f"ERROR: Cannot find INITIAL.txt (looked in {INI_PATH})", file=sys.stderr)
    sys.exit(1)

info = np.loadtxt(str(INI_PATH))
N    = int(info[0])
k1   = 30   # mean degree layer 1
k2   = 10   # mean degree layer 2

#print(f"Generating Random-Regular graphs  N={N}  k1={k1}  k2={k2}…")

#G1_2=nx.random_regular_graph(k1,N)
#G2_2=nx.random_regular_graph(k2,N)
#while nx.is_connected(G1_2)!=True:
#  G1_2=nx.random_regular_graph(k1,N)
#while nx.is_connected(G2_2)!=True:
#  G2_2=nx.random_regular_graph(k2,N)
#print("Graphs generated")




print(f"Generating Erdos-Renyi graphs  N={N}  k1={k1}  k2={k2}…")

# Generate connected graphs
G1 = nx.erdos_renyi_graph(N, k1 / N)
while not nx.is_connected(G1):
    G1 = nx.erdos_renyi_graph(N, k1 / N)

G2 = nx.erdos_renyi_graph(N, k2 / N)
while not nx.is_connected(G2):
    G2 = nx.erdos_renyi_graph(N, k2 / N)

A1 = nx.to_numpy_array(G1)
A2 = nx.to_numpy_array(G2)

# Build compact adjacency list (Fortran 1-based indices)
Adja = np.zeros((2 * N, N), dtype=int)

for i in range(N):
    nbrs = np.where(A1[i] == 1)[0]
    Adja[i, 0] = len(nbrs)
    Adja[i, 1:len(nbrs) + 1] = nbrs + 1          # 1-based

for i in range(N):
    nbrs = np.where(A2[i] == 1)[0]
    Adja[N + i, 0] = len(nbrs)
    Adja[N + i, 1:len(nbrs) + 1] = N + nbrs + 1  # 1-based, offset by N

out_path = ROOT / "adjacency_matrix.txt"
np.savetxt(str(out_path), Adja, fmt="%d")
print(f"Adjacency matrix saved to {out_path}")
