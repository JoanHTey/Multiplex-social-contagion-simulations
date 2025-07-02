import numpy as np
import subprocess
import h5py
import os

def write_initial_file(params, filename='INITIAL.txt'):
    with open(filename, 'w') as file:
        file.write(' '.join(map(str, params)) + '\n')

def run_simulation():
    # Run the Fortran executable
    result = subprocess.run([r'scripts\den_script.bat'], capture_output=True, text=True)
    if result.returncode != 0:
        print("Error running simulation:", result.stderr)
    else:
        print("Simulation completed successfully.")


# Define the range of parameters you want to test
N = 5000
DIM = 2
STIME= 40000
ETA = 1
M = 1000
GAMMA1 = 100000
GAMMA2 = 100000
INPR = 0.034
DISCTIME = 20000
T = STIME-DISCTIME
den = np.zeros([50,2])
i=0

filename='simulation/output.bin'
for ETA in np.logspace(0,1.778,50):
    params = [N, DIM, STIME, ETA, INPR, M, GAMMA1, GAMMA2, DISCTIME]
    write_initial_file(params)
    run_simulation()
    matrix = np.fromfile(filename, dtype=np.int32)
    matrix = matrix.reshape((T, N*2+2), order='C')  # Fortran order
    matrix = matrix[:, 1:N*2+1] 
    den[i,0] = np.mean(matrix[:,0:N])
    den[i,1] = np.mean(matrix[:,N:N*2])
    i=i+1
np.save('Data/density_beta0034_N5000',den)
