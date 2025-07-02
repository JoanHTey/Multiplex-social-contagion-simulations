import numpy as np
import subprocess
import os
import time
import shutil

def personcorr_aprox3(filename,h, b, mu,N):
    mat = np.loadtxt(filename, dtype=float)
    
    mean_col1 = np.mean(mat[1, N * (b - 1):N * b])
    mean_col0 = np.mean(mat[0, N * (b - 1):N * b])
    
    term1 = (1 - mean_col1) + mu * mean_col1 * (mean_col1 * (1 - mu)) ** h
    term2 = 1 - mean_col1 * (1 - mu)
    term3 = mean_col0 - mean_col0 ** 2
    
    x = ((term1 / term2) * mean_col0 - mean_col0 ** 2) / term3
    #x = (term1 / term2) * mean_col0
    return x

def write_initial_file(params, filename='INITIAL.txt'):
    with open(filename, 'w') as file:
        file.write(' '.join(map(str, params)) + '\n')


N = 5000
DIM = 2
STIME= 20000
ETA = 0.01
M = 1000
GAMMA1 = 100000
GAMMA2 = 100000
INPR = 0.034
h = np.linspace(0, 14, 15,dtype=int)
mu=0.5

start_time = time.time()
for i in np.logspace(-2,-0.3,60):
    INPR = i
    params = [N, DIM, STIME, ETA, INPR, M, GAMMA1, GAMMA2]
    write_initial_file(params)
    shutil.copy('INITIAL.txt', os.path.join('C:/Users/Usuario/Desktop/Master/TFM/Code/simulation', 'INITIAL.txt'))
    subprocess.run(['python', r'C:/Users/Usuario/Desktop/Master/TFM/Code/matrixGeneration/adjacencymat.py'], check=True)
    result = subprocess.run([r'meq.exe'], capture_output=True, text=True)
    if result.returncode != 0:
        print("error MEQLoop:", result.stderr)
    else:
        print("Succesfull MEQLoop.")
    mat = np.loadtxt('output2INPR.txt', dtype=float)
    with open('final_corr1msteq.txt', 'a') as f:
        np.savetxt(f, [np.concatenate([[INPR], personcorr_aprox3('output2INPR.txt',h, 1, mu,N)])], fmt='%.6f')
    with open('final_corr2msteq.txt', 'a') as f:
        np.savetxt(f, [np.concatenate([[INPR], personcorr_aprox3('output2INPR.txt',h, 2, mu,N)])], fmt='%.6f')
    with open('den1msteq.txt', 'a') as f:
        np.savetxt(f, [np.concatenate([[INPR], [np.mean(mat[0, N * (1 - 1):N * 1])]])], fmt='%.6f')
    with open('den2msteq.txt', 'a') as f:
        np.savetxt(f, [np.concatenate([[INPR], [np.mean(mat[0, N * (2 - 1):N * 2])]])], fmt='%.6f')
print("Time taken for the loop: ", time.time() - start_time)