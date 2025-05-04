import numpy as np
import subprocess
import h5py
import os

def write_initial_file(params, filename='INITIAL.txt'):
    with open(filename, 'w') as file:
        file.write(' '.join(map(str, params)) + '\n')

def run_simulation():
    # Run the Fortran executable
    result = subprocess.run([r'scripts\sync_script_corr.bat'], capture_output=True, text=True)
    if result.returncode != 0:
        print("Error running simulation:", result.stderr)
    else:
        print("Simulation completed successfully.")
def corrcalc(X1,X2):
    MX1 = np.mean(X1)
    MX2 = np.mean(X2)
    MX1X2 = np.mean(X1*X2)
    sigX1 = np.std(X1)
    sigX2 = np.std(X2)
    corr = (MX1X2 - MX1*MX2)/(sigX1*sigX2)
    return corr


def main():
    # Define the range of parameters you want to test
    N = 1000
    DIM = 2
    STIME= 20000
    ETA = 0.01
    M = 1000
    GAMMA1 = 100000
    GAMMA2 = 100000
    INPR = np.loadtxt('INPR.txt', dtype=float)
    
    # Check if the files exist and delete them if they do
    files_to_check = [
        r'correlationCalc\samelayer1_t0.h5', r'correlationCalc\samelayer1_th.h5', r'correlationCalc\samelayer2_t0.h5', r'correlationCalc\samelayer2_th.h5',
        r'correlationCalc\dif_samenode_layer1_t0.h5', r'correlationCalc\dif_samenode_layer2_th.h5', r'correlationCalc\samelayer1_difnode_t0.h5',
        r'correlationCalc\samelayer1_difnode_th.h5', r'correlationCalc\difflayer_difnode_t0.h5', r'correlationCalc\difflayer_difnode_th.h5',
        r'correlationCalc\samelayer2_difnode_t0.h5', r'correlationCalc\samelayer2_difnode_th.h5'
    ]
    
    for file in files_to_check:
        if os.path.exists(file):
            os.remove(file)

    for file in files_to_check:
        with h5py.File(file, "w") as f:
            f.create_dataset(
            "X", shape=(0, 15), maxshape=(None, 15),
            dtype='int', chunks=True
        )
            
    params = [N, DIM, STIME, ETA, INPR, M, GAMMA1, GAMMA2]
    write_initial_file(params)

    for i in range(0,10):
        run_simulation()

    # Read the output files and calculate the correlation
    data_dict = {}

    for i, file in enumerate(files_to_check, start=1):
        with h5py.File(file, "r") as f:
             data_dict[f"X{i}"] = f["X"][:]

    # Example usage:
    X1 = data_dict["X1"]
    X2 = data_dict["X2"]
    X3 = data_dict["X3"]
    X4 = data_dict["X4"]
    X5 = data_dict["X5"]
    X6 = data_dict["X6"]
    X7 = data_dict["X7"]
    X8 = data_dict["X8"]
    X9 = data_dict["X9"]
    X10 = data_dict["X10"]
    X11 = data_dict["X11"]
    X12 = data_dict["X12"]



    size = X1.shape
    corr1 = np.zeros(size[1])
    corr2 = np.zeros(size[1])
    corr3 = np.zeros(size[1])
    corr4 = np.zeros(size[1])
    corr5 = np.zeros(size[1])
    corr6 = np.zeros(size[1])

    for i in range(size[1]):
        corr1[i] = corrcalc(X1[:,i],X2[:,i])
        corr2[i] = corrcalc(X3[:,i],X4[:,i])
        corr3[i] = corrcalc(X5[:,i],X6[:,i])
        corr4[i] = corrcalc(X7[:,i],X8[:,i])
        corr5[i] = corrcalc(X9[:,i],X10[:,i])
        corr6[i] = corrcalc(X11[:,i],X12[:,i])


    

    with open('final_corr1.txt', 'a') as f:
        np.savetxt(f, [np.concatenate([[INPR], corr1.flatten()])], fmt='%.6f')
    with open('final_corr2.txt', 'a') as f:
        np.savetxt(f, [np.concatenate([[INPR], corr2.flatten()])], fmt='%.6f')
    with open('final_corr3.txt', 'a') as f:
        np.savetxt(f, [np.concatenate([[INPR], corr3.flatten()])], fmt='%.6f')
    with open('final_corr4.txt', 'a') as f:
        np.savetxt(f, [np.concatenate([[INPR], corr4.flatten()])], fmt='%.6f')
    with open('final_corr5.txt', 'a') as f:
        np.savetxt(f, [np.concatenate([[INPR], corr5.flatten()])], fmt='%.6f')
    with open('final_corr6.txt', 'a') as f:
        np.savetxt(f, [np.concatenate([[INPR], corr6.flatten()])], fmt='%.6f')
if __name__ == "__main__":
    main()