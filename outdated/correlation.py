import numpy as np
import matplotlib.pyplot as plt

def read_output_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Determine the number of time steps (T) and the number of elements (NTOT)
    T = len(lines)
    NTOT = len(lines[0].split())

    # Initialize the matrix
    matrix = np.zeros((T, NTOT), dtype=int)

    # Fill the matrix with data from the file
    for i, line in enumerate(lines):
        matrix[i, :] = np.array(line.split(), dtype=int)

    return matrix,T,NTOT

if __name__ == "__main__":
    filename = 'output.txt'
    matrix,T,NTOT = read_output_file(filename)
    SK = 1000
    TOTden = np.mean(matrix[:,:NTOT],axis=1) 
    L1den = np.mean(matrix[:,:NTOT//2],axis=1)
    L2den = np.mean(matrix[:,NTOT//2:NTOT],axis=1)

    TOTdenF = np.mean(TOTden[SK:])
    L1denF = np.mean(L1den[SK:])
    L2denF = np.mean(L2den[SK:])
    
    TOTsus = (np.mean(TOTden[SK:]**2)-TOTdenF**2)/TOTdenF
    L1sus = (np.mean(L1den[SK:]**2)-L1denF**2)/L1denF 
    L2sus = (np.mean(L2den[SK:]**2)-L2denF**2)/L2denF
              
    print(TOTdenF, L1denF, L2denF, TOTsus, L1sus, L2sus)    
    with open('results.txt', 'a') as file:
        file.write(f"{TOTdenF} {L1denF} {L2denF} {TOTsus} {L1sus} {L2sus}\n")