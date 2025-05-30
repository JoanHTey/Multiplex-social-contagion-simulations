
print("Generating graphs")
import networkx as nx
import numpy as np    

info = np.loadtxt('C:/Users/Usuario/Desktop/Master/TFM/Code/simulation/INITIAL.txt')
k1=30
k2=10
N=int(info[0])
print("Generating graphs")
G1_2=nx.random_regular_graph(k1,N)
G2_2=nx.random_regular_graph(k2,N)
while nx.is_connected(G1_2)!=True:
  G1_2=nx.random_regular_graph(k1,N)
while nx.is_connected(G2_2)!=True:
  G2_2=nx.random_regular_graph(k2,N)
print("Graphs generated")

# G1_2=nx.erdos_renyi_graph(N,k1/N)
# G2_2=nx.erdos_renyi_graph(N,k2/N)
# while nx.is_connected(G1_2)!=True:
#         G1_2=nx.erdos_renyi_graph(N,k1/N)
# while nx.is_connected(G2_2)!=True:
#         G2_2=nx.erdos_renyi_graph(N,k2/N)

    # Defining adjacency matrix
A_1_2 = nx.to_numpy_array(G1_2)
A_2_2 = nx.to_numpy_array(G2_2)



positions1 = np.argwhere(A_1_2 == 1)
positions2 = np.argwhere(A_2_2 == 1)

Adja = np.zeros((2*N,N))

# Fill Adja matrix for G1_2
for i in range(N):
    neighbors = np.where(A_1_2[i] == 1)[0]
    Adja[i, 0] = len(neighbors)
    Adja[i, 1:len(neighbors)+1] = neighbors + 1  # +1 to match Fortran 1-based indexing

# Fill Adja matrix for G2_2
for i in range(N):
    neighbors = np.where(A_2_2[i] == 1)[0]
    Adja[N+i, 0] = len(neighbors)
    Adja[N+i, 1:len(neighbors)+1] = N+neighbors + 1  # +1 to match Fortran 1-based indexing
np.savetxt('adjacency_matrix.txt', Adja, fmt='%d')
print("Adjacency matrix saved")
