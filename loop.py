import numpy as np
import subprocess
import os

def write_initial_file(params, filename='INITIAL.txt'):
    with open(filename, 'w') as file:
        file.write(' '.join(map(str, params)) + '\n')

def run_simulation():
    # Run the Fortran executable
    result = subprocess.run(['scripts\sync_script.bat'], capture_output=True, text=True)
    if result.returncode != 0:
        print("Error running simulation:", result.stderr)
    else:
        print("Simulation completed successfully.")

def main():
    # Define the range of parameters you want to test
    N = 1000
    DIM = 2
    STIME= 10000
    ETA = 50
    INPR_values = np.logspace(-2,-0.5,100)
    M = 1000
    GAMMA1 = 100000
    GAMMA2 = 100000

    # Run the simulation for each combination of parameters
    t=1
    for INPR in INPR_values:
        params = [N, DIM, STIME, ETA, INPR, M, GAMMA1, GAMMA2]
        write_initial_file(params)
        run_simulation()
        print(t)
        t+=1
if __name__ == "__main__":
    main()