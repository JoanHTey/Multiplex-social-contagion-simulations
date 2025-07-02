import numpy as np


def write_initial_file(params, filename='INITIAL.txt'):
    with open(filename, 'w') as file:
        file.write(' '.join(map(str, params)) + '\n')


N = 10000
DIM = 2
STIME= 5000
ETA = 1
M = 1000
GAMMA1 = 100000
GAMMA2 = 100000
INPR = 0
DISCTIME = 50
LAGS = 50


params = [N, DIM, STIME, ETA, INPR, M, GAMMA1, GAMMA2, DISCTIME,LAGS]
write_initial_file(params)

betas = np.linspace(0.01,0.15, 5)
np.savetxt('betas.txt', betas, fmt='%.6f')