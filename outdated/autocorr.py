#%%
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
        matrix[i, :] = np.array(line.split(), dtype=float)

    return matrix, T, NTOT

def autocorrelation(x, lag=1):
    x = np.array(x)
    return np.corrcoef(x[:-lag], x[lag:])[0, 1]

def personcorr_aprox1(h, psi, mu):
    x = ((psi + (1 - mu * (1 + psi)) ** h) / (1 + psi) * psi - psi ** 2) / (psi - psi ** 2)
    return x

def personcorr_aprox2(h, ps, mu,eps,beta,N):
    
    eigenvalue1 = 30
    eigenvalue2 = 10
    
    v1 = 1/N**(1/2)
    v2 = eps/(eigenvalue1-eigenvalue2)*v1
    
    A = (1-mu)**2*v1*v2
    B = (1-mu)*(v1+v2)*(mu/(beta*eigenvalue1)-N*(1-mu)*v1*v2)
    C = mu**2/(beta**2*eigenvalue1**2)-N*mu/(beta*eigenvalue1)*(v1**2+v2**2)
    
    x1 = (-B + (B**2-4*A*C)**(1/2))/(2*A)
    x2 = (-B - (B**2-4*A*C)**(1/2))/(2*A)
    
    psi = x1*ps
    
    x = (((psi + (1 - mu * psi) * ((1 - mu) + mu ** 2 * psi - mu * psi) ** h) / (1 + psi - mu * psi)) * psi - psi ** 2) / (psi - psi ** 2)
    #x = ((psi + (1 - mu * psi) * ((1 - mu) + mu ** 2 * psi - mu * psi) ** h) / (1 + psi - mu * psi)) * psi 
    return x

def personcorr_aprox3(h, b, mu):
    mat = np.loadtxt('output2.txt')
    
    mean_col1 = np.mean(mat[1, 1000 * (b - 1):1000 * b])
    mean_col0 = np.mean(mat[0, 1000 * (b - 1):1000 * b])
    
    term1 = (1 - mean_col1) + mu * mean_col1 * (mean_col1 * (1 - mu)) ** h
    term2 = 1 - mean_col1 * (1 - mu)
    term3 = mean_col0 - mean_col0 ** 2
    
    x = ((term1 / term2) * mean_col0 - mean_col0 ** 2) / term3
    #x = (term1 / term2) * mean_col0
    return x

if __name__ == "__main__":
    filename = 'output.txt'
    matrix, T, NTOT = read_output_file(filename)
    
    SK = 1000
    mu = 0.5
    h = np.linspace(0, 14, 15)
    
    # Compute ACF for multiple lags and average for the first 1000 particles
    corr_val = np.zeros(15)
    psi = 1 / 1000
    
    for i in range(15):
        for j in range(1000):
            #k = np.corrcoef(matrix[SK:T - i, j], matrix[SK + i:T, j])
            #corr_val[i] += k[1, 0]
            k = np.mean(matrix[SK:T - i, j]*matrix[SK + i:T, j])
            corr_val[i] += k
            
        corr_val[i] = corr_val[i] / 1000

    print(corr_val)
    plt.plot(corr_val, label='sim. Corr.')
    plt.plot(h, personcorr_aprox3(h, 1, mu), color='black', label='Theo. Corr.')
    plt.ylabel('Autocorr')
    plt.xlabel('lag')
    plt.legend(loc='upper right')
    plt.show()
    
    # Compute ACF for multiple lags and average for the second 1000 particles
    corr_val = np.zeros(15)
    psi = 0.000001
    
    for i in range(15):
        g = 0
        for j in range(1000, 2000):
            
            #k = np.corrcoef(matrix[SK:T - i, j], matrix[SK + i:T, j])
            #if not np.isnan(k[1, 0]):
            #    corr_val[i] += k[1, 0]
            #    g += 1
            
            k = np.mean(matrix[SK:T - i, j]*matrix[SK + i:T, j])
            if not np.isnan(k):
                corr_val[i] += k
                g += 1
                    
                    
        corr_val[i] = corr_val[i] / g

    print(corr_val)
    plt.plot(corr_val, label='sim. Corr.')
    plt.plot(h, personcorr_aprox3(h, 2, mu), color='black', label='Theo. Corr.')
    plt.ylabel('Autocorr')
    plt.xlabel('lag')
    plt.legend(loc='upper right')
    plt.show()
    
#%%
mu = 0.5
for i in np.logspace(-4, 0, num=10):
    plt.plot(h, personcorr_aprox2(h, i, mu), label=rf'$\Psi$ {i}')
plt.legend(loc='upper right')
plt.show()


#%%

corr_val = np.mean(np.loadtxt('corr_val_output1.txt'),axis=0)

plt.plot(corr_val, label='sim. Corr.')
plt.plot(h, personcorr_aprox3(h, 1, mu), color='black', label='Theo. Corr.')
plt.ylabel('Autocorr')
plt.xlabel('lag')
plt.legend(loc='upper right')
plt.show()


corr_val = np.mean(np.loadtxt('corr_val_output2.txt'),axis=0)

plt.plot(corr_val, label='sim. Corr.')
plt.plot(h, personcorr_aprox3(h, 2, mu), color='black', label='Theo. Corr.')
plt.ylabel('Autocorr')
plt.xlabel('lag')
plt.legend(loc='upper right')
plt.show()

#%%
corr_val = np.mean(np.loadtxt('corr_val_output1.txt'),axis=0)

plt.plot(corr_val, label='sim. Corr.')
plt.plot(h, personcorr_aprox2(h, 1/1000**(1/2),mu,0.01,0.05,1000), color='black', label='Theo. Corr.')
plt.ylabel('Autocorr')
plt.xlabel('lag')
plt.legend(loc='upper right')
plt.show()


corr_val = np.mean(np.loadtxt('corr_val_output2.txt'),axis=0)

plt.plot(corr_val, label='sim. Corr.')
plt.plot(h, personcorr_aprox2(h, 0.01/(20*1000**(1/2)),mu,0.01,0.05,1000), color='black', label='Theo. Corr.')
plt.ylabel('Autocorr')
plt.xlabel('lag')
plt.legend(loc='upper right')
plt.show()
#%%
N=1000
ps=0.01/(20*N**(1/2))
ps = 1/N**(1/2)
beta=0.05
eps =0.01
eigenvalue1 = 30
eigenvalue2 = 10

v1 = 1/N**(1/2)
v2 = eps/(eigenvalue1-eigenvalue2)*v1

A = (1-mu)**2*v1*v2
B = (1-mu)*(v1+v2)*(mu/(beta*eigenvalue1)-N*(1-mu)*v1*v2)
C = mu**2/(beta**2*eigenvalue1**2)-N*mu/(beta*eigenvalue1)*(v1**2+v2**2)

x1 = (-B + (B**2-4*A*C)**(1/2))/(2*A)
x2 = (-B - (B**2-4*A*C)**(1/2))/(2*A)

psi = x1*ps

x = (((psi + (1 - mu * psi) * ((1 - mu) + mu ** 2 * psi - mu * psi) ** h) / (1 + psi - mu * psi)) * psi - psi ** 2) / (psi - psi ** 2)

 