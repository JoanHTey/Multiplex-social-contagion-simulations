#%%
import h5py
import numpy as np
from joblib import Parallel, delayed




if __name__ == "__main__":
    filename = 'output.bin'
    ini = r'C:\Users\Usuario\Desktop\Master\TFM\Code\simulation\INITIAL.txt'
    matrix =  np.fromfile(filename, dtype=np.int32)  # or np.int8 if you know the type
    init = np.loadtxt(ini)
    NTOT = int(init[0]*2)
    T = int(init[2]- init[8])  # Total time steps minus initial time
    matrix = matrix.reshape((T, NTOT+2), order='C')  # Fortran order
    matrix = matrix[:, 1:NTOT+1] 
    adjacency_matrix = np.loadtxt('adjacency_matrix.txt', dtype=int)

    N = NTOT // 2
    mu = 0.5
    h = np.linspace(0, 14, 15)
    num_rands_i = 500
    num_rands = num_rands_i*N
    gap = T//500

    # Compute ACF for multiple lags and average for the first 1000 particles
    k1 = np.zeros((num_rands, 15), dtype=int)
    k2 = np.zeros((num_rands, 15), dtype=int)
    k3 = np.zeros((num_rands, 15), dtype=int)
    k4 = np.zeros((num_rands, 15), dtype=int)
    k5 = np.zeros((num_rands, 15), dtype=int)
    k6 = np.zeros((num_rands, 15), dtype=int)
    k7 = np.zeros((num_rands, 15), dtype=int)
    k8 = np.zeros((num_rands, 15), dtype=int)
    k9 = np.zeros((num_rands, 15), dtype=int)
    k10 = np.zeros((num_rands, 15), dtype=int)
    k11 = np.zeros((num_rands, 15), dtype=int)
    k12 = np.zeros((num_rands, 15), dtype=int)

    def process_j(j):
        # Local arrays for a single j
        local_k1 = np.zeros((num_rands_i, 15), dtype=int)
        local_k2 = np.zeros((num_rands_i, 15), dtype=int)
        local_k3 = np.zeros((num_rands_i, 15), dtype=int)
        local_k4 = np.zeros((num_rands_i, 15), dtype=int)
        local_k5 = np.zeros((num_rands_i, 15), dtype=int)
        local_k6 = np.zeros((num_rands_i, 15), dtype=int)
        local_k7 = np.zeros((num_rands_i, 15), dtype=int)
        local_k8 = np.zeros((num_rands_i, 15), dtype=int)
        local_k9 = np.zeros((num_rands_i, 15), dtype=int)
        local_k10 = np.zeros((num_rands_i, 15), dtype=int)
        local_k11 = np.zeros((num_rands_i, 15), dtype=int)
        local_k12 = np.zeros((num_rands_i, 15), dtype=int)

        for i in range(15):

            RANDOM1 = np.linspace(0,num_rands_i-1,num_rands_i, dtype=int)*gap
            RANDOM2 = np.linspace(0,num_rands_i-1,num_rands_i, dtype=int)*gap
            RANDOM3 = np.linspace(0,num_rands_i-1,num_rands_i, dtype=int)*gap
            RANDOM4 = np.linspace(0,num_rands_i-1,num_rands_i, dtype=int)*gap
            RANDOM5 = np.linspace(0,num_rands_i-1,num_rands_i, dtype=int)*gap
            RANDOM6 = np.linspace(0,num_rands_i-1,num_rands_i, dtype=int)*gap

            # First layer correlation
            local_k1[:, i] = matrix[RANDOM1, j]
            local_k2[:, i] = matrix[RANDOM1 + i, j]

            # Second layer correlation
            local_k3[:, i] = matrix[RANDOM2, N + j]
            local_k4[:, i] = matrix[RANDOM2 + i, N + j]

            # Correlation between layers same node
            local_k5[:, i] = matrix[RANDOM3, j]
            local_k6[:, i] = matrix[RANDOM3 + i, N + j]

            # Correlation same layer diff node
            local_k7[:, i] = matrix[RANDOM4, j]
            local_k8[:, i] = matrix[RANDOM4 + i, adjacency_matrix[j, 1] - 1]

            # Correlation diff layer diff node
            local_k9[:, i] = matrix[RANDOM5, j]
            local_k10[:, i] = matrix[RANDOM5 + i, adjacency_matrix[N + j, 1] - 1]

            # Correlation same layer diff node (second layer)
            local_k11[:, i] = matrix[RANDOM6, N + j]
            local_k12[:, i] = matrix[RANDOM6 + i, adjacency_matrix[N + j, 1] - 1]

        return (
            local_k1, local_k2, local_k3, local_k4, local_k5, local_k6,
            local_k7, local_k8, local_k9, local_k10, local_k11, local_k12
        )

    # Run in parallel
    results = Parallel(n_jobs=10)(delayed(process_j)(j) for j in range(N))

    # Stack each group of results
    k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12 = [
        np.vstack(arr) for arr in zip(*results)
]
    # Write results to files
    def save_results(filename, data):
        with h5py.File(filename, "a") as f:
            dset = f["X"]
            old_size = dset.shape[0]
            dset.resize(old_size + data.shape[0], axis=0)
            dset[old_size:] = data

    filenames = [
        'samelayer1_t0.h5', 'samelayer1_th.h5', 'samelayer2_t0.h5', 'samelayer2_th.h5',
        'dif_samenode_layer1_t0.h5', 'dif_samenode_layer2_th.h5', 'samelayer1_difnode_t0.h5',
        'samelayer1_difnode_th.h5', 'difflayer_difnode_t0.h5', 'difflayer_difnode_th.h5',
        'samelayer2_difnode_t0.h5', 'samelayer2_difnode_th.h5'
    ]

    data_arrays = [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12]

    for i in range(12):
        save_results(filenames[i], data_arrays[i])


