#%%
import h5py
import numpy as np
from joblib import Parallel, delayed




if __name__ == "__main__":
    filename = 'output.txt'
    matrix = np.loadtxt(filename, dtype=int)
    adjacency_matrix = np.loadtxt('adjacency_matrix.txt', dtype=int)

    T,NTOT = matrix.shape

    SK = 5000
    mu = 0.5
    h = np.linspace(0, 14, 15)
    num_rands = 500*1000
    
    # Compute ACF for multiple lags and average for the first 1000 particles
    k1 = np.zeros((500 * 1000, 15), dtype=int)
    k2 = np.zeros((500 * 1000, 15), dtype=int)
    k3 = np.zeros((500 * 1000, 15), dtype=int)
    k4 = np.zeros((500 * 1000, 15), dtype=int)
    k5 = np.zeros((500 * 1000, 15), dtype=int)
    k6 = np.zeros((500 * 1000, 15), dtype=int)
    k7 = np.zeros((500 * 1000, 15), dtype=int)
    k8 = np.zeros((500 * 1000, 15), dtype=int)
    k9 = np.zeros((500 * 1000, 15), dtype=int)
    k10 = np.zeros((500 * 1000, 15), dtype=int)
    k11 = np.zeros((500 * 1000, 15), dtype=int)
    k12 = np.zeros((500 * 1000, 15), dtype=int)

    def process_j(j):
        # Local arrays for a single j
        local_k1 = np.zeros((500, 15), dtype=int)
        local_k2 = np.zeros((500, 15), dtype=int)
        local_k3 = np.zeros((500, 15), dtype=int)
        local_k4 = np.zeros((500, 15), dtype=int)
        local_k5 = np.zeros((500, 15), dtype=int)
        local_k6 = np.zeros((500, 15), dtype=int)
        local_k7 = np.zeros((500, 15), dtype=int)
        local_k8 = np.zeros((500, 15), dtype=int)
        local_k9 = np.zeros((500, 15), dtype=int)
        local_k10 = np.zeros((500, 15), dtype=int)
        local_k11 = np.zeros((500, 15), dtype=int)
        local_k12 = np.zeros((500, 15), dtype=int)

        for i in range(15):
            RANDOM1 = np.random.choice(np.arange(SK, T - i), size=500, replace=False)
            RANDOM2 = np.random.choice(np.arange(SK, T - i), size=500, replace=False)
            RANDOM3 = np.random.choice(np.arange(SK, T - i), size=500, replace=False)
            RANDOM4 = np.random.choice(np.arange(SK, T - i), size=500, replace=False)
            RANDOM5 = np.random.choice(np.arange(SK, T - i), size=500, replace=False)
            RANDOM6 = np.random.choice(np.arange(SK, T - i), size=500, replace=False)

            # First layer correlation
            local_k1[:, i] = matrix[RANDOM1, j]
            local_k2[:, i] = matrix[RANDOM1 + i, j]

            # Second layer correlation
            local_k3[:, i] = matrix[RANDOM2, 1000 + j]
            local_k4[:, i] = matrix[RANDOM2 + i, 1000 + j]

            # Correlation between layers same node
            local_k5[:, i] = matrix[RANDOM3, j]
            local_k6[:, i] = matrix[RANDOM3 + i, 1000 + j]

            # Correlation same layer diff node
            local_k7[:, i] = matrix[RANDOM4, j]
            local_k8[:, i] = matrix[RANDOM4 + i, adjacency_matrix[j, 1] - 1]

            # Correlation diff layer diff node
            local_k9[:, i] = matrix[RANDOM5, j]
            local_k10[:, i] = matrix[RANDOM5 + i, adjacency_matrix[1000 + j, 1] - 1]

            # Correlation same layer diff node (second layer)
            local_k11[:, i] = matrix[RANDOM6, 1000 + j]
            local_k12[:, i] = matrix[RANDOM6 + i, adjacency_matrix[1000 + j, 1] - 1]

        return (
            local_k1, local_k2, local_k3, local_k4, local_k5, local_k6,
            local_k7, local_k8, local_k9, local_k10, local_k11, local_k12
        )

    # Run in parallel
    results = Parallel(n_jobs=10)(delayed(process_j)(j) for j in range(1000))

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


# %%
