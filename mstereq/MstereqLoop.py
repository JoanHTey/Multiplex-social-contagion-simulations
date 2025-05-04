import numpy as np
import subprocess
import os
import time

def write_initial_file(params, filename='INITIAL.txt'):
    with open(filename, 'w') as file:
        file.write(' '.join(map(str, params)) + '\n')

N = 1000
DIM = 2
STIME= 20000
ETA = 0.01
M = 1000
GAMMA1 = 100000
GAMMA2 = 100000


start_time = time.time()
for i in np.logspace(-2,-0.5,25):
    INPR = i
    params = [N, DIM, STIME, ETA, INPR, M, GAMMA1, GAMMA2]
    write_initial_file(params)
    subprocess.run(['python', r'C:\Users\Usuario\Desktop\Master\TFM\Code\matrixGeneration\adjacencymat.py'], check=True)
    result = subprocess.run([r'meq.exe'], capture_output=True, text=True)
    if result.returncode != 0:
        print("error MEQLoop:", result.stderr)
    else:
        print("Succesfull MEQLoop.")

print("Time taken for the loop: ", time.time() - start_time)