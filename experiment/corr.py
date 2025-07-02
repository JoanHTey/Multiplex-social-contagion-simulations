import numpy as np


filename = 'output.bin'
matrix =  np.fromfile(filename, dtype=np.int32)
ini = 'INITIAL.txt'
init = np.loadtxt(ini)

adja = np.loadtxt('adjacency_matrix.txt',dtype=np.int32)
NTOT = int(init[0]*2)
N = NTOT//2
dist = int(init[8])
lags = int(init[9])
beta_value = 5
T=100*lags*beta_value
Ts = 10*beta_value
matrix = matrix.reshape((T, NTOT+3), order='C')  # Fortran order
matrix = matrix[:, 1:NTOT+2] 
# Sum along the time axis, wrapping every 15000 time steps, up to 40000 steps
wrap = lags*100
max_steps = 100*dist
rlags = 3
corr1 = np.zeros([Ts,rlags])
corr2 = np.zeros([Ts,rlags])
corr3 = np.zeros([Ts,rlags])
corr4 = np.zeros([Ts,rlags])
corr5 = np.zeros([Ts,rlags])
corr6 = np.zeros([Ts,rlags])
ecorr1 = np.zeros([Ts,rlags])
ecorr2 = np.zeros([Ts,rlags])
ecorr3 = np.zeros([Ts,rlags])
ecorr4 = np.zeros([Ts,rlags])
ecorr5 = np.zeros([Ts,rlags])
ecorr6 = np.zeros([Ts,rlags])
den1 = np.zeros(Ts)
den2 = np.zeros(Ts)

for i in range(0,beta_value-1):
    matrix[wrap*i:wrap*(i+1),0] = matrix[wrap*i:wrap*(i+1),0] + max_steps*i
lags= lags*10
for i in range(0,Ts):
    for j in range(0,rlags):
        corr1[i,j] = (np.mean(matrix[i*lags:(i+1)*lags-j,1:(N+1)] * matrix[(i*lags+j):(i+1)*lags,1:(N+1)])-np.mean(matrix[i*lags:(i+1)*lags-j,1:(N+1)])*np.mean(matrix[(i*lags+j):(i+1)*lags,1:(N+1)]))/(np.std(matrix[i*lags:(i+1)*lags-j,1:(N+1)])*np.std(matrix[(i*lags+j):(i+1)*lags,1:(N+1)]))
        corr2[i,j] = (np.mean(matrix[i*lags:(i+1)*lags-j,N+1:NTOT+1]*matrix[(i*lags+j):(i+1)*lags,N+1:NTOT+1])-np.mean(matrix[i*lags:(i+1)*lags-j,N+1:NTOT+1])*np.mean(matrix[(i*lags+j):(i+1)*lags,N+1:NTOT+1]))/(np.std(matrix[i*lags:(i+1)*lags-j,N+1:NTOT+1])*np.std(matrix[(i*lags+j):(i+1)*lags,N+1:NTOT+1]))
        corr3[i,j] = (np.mean(matrix[i*lags:(i+1)*lags-j,1:N+1]*matrix[(i*lags+j):(i+1)*lags,N+1:NTOT+1])-np.mean(matrix[i*lags:(i+1)*lags-j,1:N+1])*np.mean(matrix[(i*lags+j):(i+1)*lags,N+1:NTOT+1]))/(np.std(matrix[i*lags:(i+1)*lags-j,1:N+1])*np.std(matrix[(i*lags+j):(i+1)*lags,N+1:NTOT+1]))
        corr4[i,j] = (np.mean( matrix[i*lags:(i+1)*lags-j,1:N+1] * matrix[(i*lags+j):(i+1)*lags,adja[:N,1]]) - np.mean( matrix[i*lags:(i+1)*lags-j,1:N+1])*np.mean(matrix[(i*lags+j):(i+1)*lags,adja[:N,1]]))/(np.std( matrix[i*lags:(i+1)*lags-j,1:N+1])*np.std(matrix[(i*lags+j):(i+1)*lags,adja[:N,1]]))
        corr5[i,j] = (np.mean(matrix[i*lags:(i+1)*lags-j,N+1:NTOT+1] * matrix[(i*lags+j):(i+1)*lags,adja[N:NTOT,1]])-np.mean(matrix[i*lags:(i+1)*lags-j,N+1:NTOT+1])*np.mean(matrix[(i*lags+j):(i+1)*lags,adja[N:NTOT,1]]))/(np.std(matrix[i*lags:(i+1)*lags-j,N+1:NTOT+1])*np.std(matrix[(i*lags+j):(i+1)*lags,adja[N:NTOT,1]]))
        corr6[i,j] = (np.mean(matrix[i*lags:(i+1)*lags-j,1:N+1]*matrix[(i*lags+j):(i+1)*lags,adja[N:NTOT,1]])-np.mean(matrix[i*lags:(i+1)*lags-j,1:N+1])*np.mean(matrix[(i*lags+j):(i+1)*lags,adja[N:NTOT,1]]))/(np.std(matrix[i*lags:(i+1)*lags-j,1:N+1])*np.std(matrix[(i*lags+j):(i+1)*lags,adja[N:NTOT,1]]))
        
        ecorr1[i,j] = np.std((np.mean(matrix[i*lags:(i+1)*lags-j,1:(N+1)] * matrix[(i*lags+j):(i+1)*lags,1:(N+1)],axis=1)-np.mean(matrix[i*lags:(i+1)*lags-j,1:(N+1)],axis=1)*np.mean(matrix[(i*lags+j):(i+1)*lags,1:(N+1)],axis=1))/(np.std(matrix[i*lags:(i+1)*lags-j,1:(N+1)],axis=1)*np.std(matrix[(i*lags+j):(i+1)*lags,1:(N+1)],axis=1)+0.0000001))
        ecorr2[i,j] = np.std((np.mean(matrix[i*lags:(i+1)*lags-j,N+1:NTOT+1]*matrix[(i*lags+j):(i+1)*lags,N+1:NTOT+1],axis=1)-np.mean(matrix[i*lags:(i+1)*lags-j,N+1:NTOT+1],axis=1)*np.mean(matrix[(i*lags+j):(i+1)*lags,N+1:NTOT+1],axis=1))/(np.std(matrix[i*lags:(i+1)*lags-j,N+1:NTOT+1],axis=1)*np.std(matrix[(i*lags+j):(i+1)*lags,N+1:NTOT+1],axis=1)+0.0000001))
        ecorr3[i,j] = np.std((np.mean(matrix[i*lags:(i+1)*lags-j,1:N+1]*matrix[(i*lags+j):(i+1)*lags,N+1:NTOT+1],axis=1)-np.mean(matrix[i*lags:(i+1)*lags-j,1:N+1],axis=1)*np.mean(matrix[(i*lags+j):(i+1)*lags,N+1:NTOT+1],axis=1))/(np.std(matrix[i*lags:(i+1)*lags-j,1:N+1],axis=1)*np.std(matrix[(i*lags+j):(i+1)*lags,N+1:NTOT+1],axis=1)+0.0000001))
        ecorr4[i,j] = np.std((np.mean( matrix[i*lags:(i+1)*lags-j,1:N+1] * matrix[(i*lags+j):(i+1)*lags,adja[:N,1]],axis=1) - np.mean( matrix[i*lags:(i+1)*lags-j,1:N+1],axis=1)*np.mean(matrix[(i*lags+j):(i+1)*lags,adja[:N,1]],axis=1))/(np.std( matrix[i*lags:(i+1)*lags-j,1:N+1],axis=1)*np.std(matrix[(i*lags+j):(i+1)*lags,adja[:N,1]],axis=1)+0.0000001))
        ecorr5[i,j] = np.std((np.mean(matrix[i*lags:(i+1)*lags-j,N+1:NTOT+1] * matrix[(i*lags+j):(i+1)*lags,adja[N:NTOT,1]],axis=1)-np.mean(matrix[i*lags:(i+1)*lags-j,N+1:NTOT+1],axis=1)*np.mean(matrix[(i*lags+j):(i+1)*lags,adja[N:NTOT,1]],axis=1))/(np.std(matrix[i*lags:(i+1)*lags-j,N+1:NTOT+1],axis=1)*np.std(matrix[(i*lags+j):(i+1)*lags,adja[N:NTOT,1]],axis=1)+0.0000001))
        ecorr6[i,j] = np.std((np.mean(matrix[i*lags:(i+1)*lags-j,1:N+1]*matrix[(i*lags+j):(i+1)*lags,adja[N:NTOT,1]],axis=1)-np.mean(matrix[i*lags:(i+1)*lags-j,1:N+1],axis=1)*np.mean(matrix[(i*lags+j):(i+1)*lags,adja[N:NTOT,1]],axis=1))/(np.std(matrix[i*lags:(i+1)*lags-j,1:N+1],axis=1)*np.std(matrix[(i*lags+j):(i+1)*lags,adja[N:NTOT,1]],axis=1)+0.0000001))   
    den1[i] = np.mean(matrix[i*lags:(i+1)*lags,1:(N+1)])
    den2[i] = np.mean(matrix[i*lags:(i+1)*lags,N+1:NTOT+1])

np.save('corr1.npy',corr1)
np.save('corr2.npy',corr2)
np.save('corr3.npy',corr3)
np.save('corr4.npy',corr4)
np.save('corr5.npy',corr5)
np.save('corr6.npy',corr6)

np.save('ecorr1.npy',ecorr1)
np.save('ecorr2.npy',ecorr2)
np.save('ecorr3.npy',ecorr3)
np.save('ecorr4.npy',ecorr4)
np.save('ecorr5.npy',ecorr5)
np.save('ecorr6.npy',ecorr6)

np.save('den1.npy',den1)
np.save('den2.npy',den2)