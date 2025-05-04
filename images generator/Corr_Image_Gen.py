import numpy as np
import matplotlib.pyplot as plt

def personcorr_aprox3(filename,h, b, mu):
    mat = np.loadtxt(filename, dtype=float)
    
    mean_col1 = np.mean(mat[1, 1000 * (b - 1):1000 * b])
    mean_col0 = np.mean(mat[0, 1000 * (b - 1):1000 * b])
    
    term1 = (1 - mean_col1) + mu * mean_col1 * (mean_col1 * (1 - mu)) ** h
    term2 = 1 - mean_col1 * (1 - mu)
    term3 = mean_col0 - mean_col0 ** 2
    
    x = ((term1 / term2) * mean_col0 - mean_col0 ** 2) / term3
    #x = (term1 / term2) * mean_col0
    return x

SamelayerSameNode1 = np.loadtxt('final_corr1.txt')
SamelayerSameNode2 = np.loadtxt('final_corr2.txt')
DiffLayerSameNode = np.loadtxt('final_corr3.txt')
SameLayerDiffNode1 = np.loadtxt('final_corr4.txt')
SamelayerDiffNode2 = np.loadtxt('final_corr6.txt')
DiffLayerDiffNode = np.loadtxt('final_corr5.txt')
h = np.linspace(0, 14, 15,dtype=int)

#%%
filename = [r'C:\Users\Usuario\Desktop\Master\TFM\Code\mstereq\output2INPR0.005.txt',
            r'C:\Users\Usuario\Desktop\Master\TFM\Code\mstereq\output2INPR0.006.txt',
            r'C:\Users\Usuario\Desktop\Master\TFM\Code\mstereq\output2INPR0.007.txt',
            r'C:\Users\Usuario\Desktop\Master\TFM\Code\mstereq\output2INPR0.008.txt',
            r'C:\Users\Usuario\Desktop\Master\TFM\Code\mstereq\output2INPR0.009.txt',
            r'C:\Users\Usuario\Desktop\Master\TFM\Code\mstereq\output2INPR0.010.txt',
            r'C:\Users\Usuario\Desktop\Master\TFM\Code\mstereq\output2INPR0.012.txt',
            r'C:\Users\Usuario\Desktop\Master\TFM\Code\mstereq\output2INPR0.014.txt',
            r'C:\Users\Usuario\Desktop\Master\TFM\Code\mstereq\output2INPR0.016.txt',
            r'C:\Users\Usuario\Desktop\Master\TFM\Code\mstereq\output2INPR0.018.txt',
            r'C:\Users\Usuario\Desktop\Master\TFM\Code\mstereq\output2INPR0.021.txt',
            r'C:\Users\Usuario\Desktop\Master\TFM\Code\mstereq\output2INPR0.024.txt',
            r'C:\Users\Usuario\Desktop\Master\TFM\Code\mstereq\output2INPR0.028.txt',
            r'C:\Users\Usuario\Desktop\Master\TFM\Code\mstereq\output2INPR0.032.txt',
            r'C:\Users\Usuario\Desktop\Master\TFM\Code\mstereq\output2INPR0.037.txt',
            r'C:\Users\Usuario\Desktop\Master\TFM\Code\mstereq\output2INPR0.043.txt',
            r'C:\Users\Usuario\Desktop\Master\TFM\Code\mstereq\output2INPR0.050.txt',
            r'C:\Users\Usuario\Desktop\Master\TFM\Code\mstereq\output2INPR0.058.txt',
            r'C:\Users\Usuario\Desktop\Master\TFM\Code\mstereq\output2INPR0.067.txt',
            r'C:\Users\Usuario\Desktop\Master\TFM\Code\mstereq\output2INPR0.077.txt',
            r'C:\Users\Usuario\Desktop\Master\TFM\Code\mstereq\output2INPR0.089.txt',
            r'C:\Users\Usuario\Desktop\Master\TFM\Code\mstereq\output2INPR0.103.txt',
            r'C:\Users\Usuario\Desktop\Master\TFM\Code\mstereq\output2INPR0.119.txt',
            r'C:\Users\Usuario\Desktop\Master\TFM\Code\mstereq\output2INPR0.137.txt',
            r'C:\Users\Usuario\Desktop\Master\TFM\Code\mstereq\output2INPR0.158.txt']


for i in range(25):
    plt.plot(h, SamelayerSameNode1[i,1:], label='SamelayerSameNode1', color='blue')
    plt.plot(h, personcorr_aprox3(filename[i],h,1,0.5), label='SamelayerSameNode1 theoric', color='black')
    plt.show()


#%%
for i in h:
    plt.title(f'Indent space {i}')
    plt.plot(SamelayerSameNode1[:,0],SamelayerSameNode1[:,i])
    plt.show()