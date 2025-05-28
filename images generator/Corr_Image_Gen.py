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

SamelayerSameNode1_ETA1 = np.loadtxt('Data/final_corr1_eta1.txt')
SamelayerSameNode2_ETA1 = np.loadtxt('Data/final_corr2_eta1.txt')
DiffLayerSameNode_ETA1 = np.loadtxt('Data/final_corr3_eta1.txt')
SameLayerDiffNode1_ETA1 = np.loadtxt('Data/final_corr4_eta1.txt')
SamelayerDiffNode2_ETA1 = np.loadtxt('Data/final_corr6_eta1.txt')
DiffLayerDiffNode_ETA1 = np.loadtxt('Data/final_corr5_eta1.txt')

# # Load the data from the text files
SamelayerSameNode1_ETA25 = np.loadtxt('Data/final_corr1_eta25.txt')
SamelayerSameNode2_ETA25 = np.loadtxt('Data/final_corr2_eta25.txt')
DiffLayerSameNode_ETA25 = np.loadtxt('Data/final_corr3_eta25.txt')
SameLayerDiffNode1_ETA25 = np.loadtxt('Data/final_corr4_eta25.txt')
DiffLayerDiffNode_ETA25 = np.loadtxt('Data/final_corr5_eta25.txt')
SamelayerDiffNode2_ETA25 = np.loadtxt('Data/final_corr6_eta25.txt')

# Load the data from the text files
SamelayerSameNode1_ETA001 = np.loadtxt('Data/final_corr1_eta0.01.txt')
SamelayerSameNode2_ETA001  = np.loadtxt('Data/final_corr2_eta0.01.txt')
DiffLayerSameNode_ETA001  = np.loadtxt('Data/final_corr3_eta0.01.txt')
SameLayerDiffNode1_ETA001  = np.loadtxt('Data/final_corr4_eta0.01.txt')
DiffLayerDiffNode_ETA001  = np.loadtxt('Data/final_corr5_eta0.01.txt')
SamelayerDiffNode2_ETA001  = np.loadtxt('Data/final_corr6_eta0.01.txt')

# Load the data from the text files
SamelayerSameNode1_BETA0034 = np.loadtxt('Data/final_corr1_beta0.034.txt')
SamelayerSameNode2_BETA0034 = np.loadtxt('Data/final_corr2_beta0.034.txt')
DiffLayerSameNode_BETA0034 = np.loadtxt('Data/final_corr3_beta0.034.txt')
SameLayerDiffNode1_BETA0034 = np.loadtxt('Data/final_corr4_beta0.034.txt')
DiffLayerDiffNode_BETA0034 = np.loadtxt('Data/final_corr5_beta0.034.txt')
SamelayerDiffNode2_BETA0034 = np.loadtxt('Data/final_corr6_beta0.034.txt')

h = np.linspace(0, 14, 15,dtype=int)

# Master equation solution
SS1MSTEQ_ETA1 = np.loadtxt('Data/final_corr1msteq_eta1.txt')
SS2MSTEQ_ETA1 = np.loadtxt('Data/final_corr2msteq_eta1.txt')

SS1MSTEQ_ETA25 = np.loadtxt('Data/final_corr1msteq_eta25.txt')
SS2MSTEQ_ETA25 = np.loadtxt('Data/final_corr2msteq_eta25.txt')

SS1MSTEQ_ETA001 = np.loadtxt('Data/final_corr1msteq_eta0.01.txt')
SS2MSTEQ_ETA001 = np.loadtxt('Data/final_corr2msteq_eta0.01.txt')

SS1MSTEQ_BETA0034 = np.loadtxt('Data/final_corr1msteq_beta0.034.txt')
SS2MSTEQ_BETA0034 = np.loadtxt('Data/final_corr2msteq_beta0.034.txt')

#%%
# for i in range(25):
#     plt.plot(h, SamelayerSameNode1[i,1:], label='SamelayerSameNode1', color='blue')
#     plt.plot(h, personcorr_aprox3(filename[i],h,1,0.5), label='SamelayerSameNode1 theoric', color='black')
#     plt.show()
#%%
plt.plot(h, SamelayerSameNode1_ETA001[6,1:], label='SamelayerSameNode1')
plt.plot(h, SamelayerSameNode2_ETA001[4,1:], label='SamelayerSameNode2')
plt.plot(h, DiffLayerSameNode_ETA001[6,1:], label='DiffLayerSameNode')
plt.plot(h, SameLayerDiffNode1_ETA001[6,1:], label='SameLayerDiffNode1')
plt.plot(h, SamelayerDiffNode2_ETA001[6,1:], label='SamelayerDiffNode2')
plt.plot(h, DiffLayerDiffNode_ETA001[6,1:], label='DiffLayerDiffNode')
plt.xlabel('time lag (h)',fontsize=14)
plt.ylabel('Correlation',fontsize=14)
plt.title(r'$\beta$=0.037',fontsize=14)
plt.legend()
plt.grid()
plt.savefig('CorrAll.png')
plt.show()
#%%

plt.title(f'Indent space {1}')
plt.plot(SamelayerSameNode1_ETA001[:,0],SamelayerSameNode1_ETA001[:,2],linestyle='',marker='o',label='SamelayerSameNode1')
plt.plot(SamelayerSameNode2_ETA001[:,0],SamelayerSameNode2_ETA001[:,2],linestyle='',marker='o',label='SamelayerSameNode2')
plt.plot(DiffLayerSameNode_ETA001[:,0],DiffLayerSameNode_ETA001[:,2],linestyle='',marker='o',label='DiffLayerSameNode') 
plt.plot(SameLayerDiffNode1_ETA001[:,0],SameLayerDiffNode1_ETA001[:,2],linestyle='',marker='o',label='SameLayerDiffNode1')
plt.plot(SamelayerDiffNode2_ETA001[:,0],SamelayerDiffNode2_ETA001[:,2],linestyle='',marker='o',label='SamelayerDiffNode2')
plt.plot(DiffLayerDiffNode_ETA001[:,0],DiffLayerDiffNode_ETA001[:,2],linestyle='',marker='o',label='DiffLayerDiffNode')
plt.plot(SS1MSTEQ_ETA001[:,0],SS1MSTEQ_ETA001[:,2],linestyle='-',label='SamelayerSameNode1 Mstereq')
plt.plot(SS2MSTEQ_ETA001[:,0],SS2MSTEQ_ETA001[:,2],linestyle='-',label='SamelayerSameNode2 Mstereq')
plt.xlabel(r'$beta$',fontsize=14)
plt.ylabel('Correlation',fontsize=14)
plt.legend()
plt.grid()
plt.savefig('CorrBeta_ETA001.png')
plt.show()

#%%
plt.title(f'Indent space {1}') 
plt.plot(SamelayerSameNode1_ETA25[:,0],SamelayerSameNode1_ETA25[:,2],linestyle='',marker='o',label='SamelayerSameNode1')
plt.plot(SamelayerSameNode2_ETA25[:,0],SamelayerSameNode2_ETA25[:,2],linestyle='',marker='o',label='SamelayerSameNode2')
plt.plot(DiffLayerSameNode_ETA25[:,0],DiffLayerSameNode_ETA25[:,2],linestyle='',marker='o',label='DiffLayerSameNode')
plt.plot(SameLayerDiffNode1_ETA25[:,0],SameLayerDiffNode1_ETA25[:,2],linestyle='',marker='o',label='SameLayerDiffNode1')
plt.plot(SamelayerDiffNode2_ETA25[:,0],SamelayerDiffNode2_ETA25[:,2],linestyle='',marker='o',label='SamelayerDiffNode2')
plt.plot(DiffLayerDiffNode_ETA25[:,0],DiffLayerDiffNode_ETA25[:,2],linestyle='',marker='o',label='DiffLayerDiffNode')
plt.plot(SS1MSTEQ_ETA25[:,0],SS1MSTEQ_ETA25[:,2],linestyle='-',label='SamelayerSameNode1 Mstereq')
plt.plot(SS2MSTEQ_ETA25[:,0],SS2MSTEQ_ETA25[:,2],linestyle='-',label='SamelayerSameNode2 Mstereq')
plt.xlabel(r'$beta$',fontsize=14)
plt.ylabel('Correlation',fontsize=14)
plt.legend()
plt.grid()
plt.savefig('CorrBeta_ETA25.png')
plt.show()
#%% 
plt.title(f'Indent space {1}')
plt.plot(SamelayerSameNode1_ETA1[:,0],SamelayerSameNode1_ETA1[:,2],linestyle='',marker='o',label='SamelayerSameNode1')  
plt.plot(SamelayerSameNode2_ETA1[:,0],SamelayerSameNode2_ETA1[:,2],linestyle='',marker='o',label='SamelayerSameNode2')
plt.plot(DiffLayerSameNode_ETA1[:,0],DiffLayerSameNode_ETA1[:,2],linestyle='',marker='o',label='DiffLayerSameNode')
plt.plot(SameLayerDiffNode1_ETA1[:,0],SameLayerDiffNode1_ETA1[:,2],linestyle='',marker='o',label='SameLayerDiffNode1')
plt.plot(SamelayerDiffNode2_ETA1[:,0],SamelayerDiffNode2_ETA1[:,2],linestyle='',marker='o',label='SamelayerDiffNode2')
plt.plot(DiffLayerDiffNode_ETA1[:,0],DiffLayerDiffNode_ETA1[:,2],linestyle='',marker='o',label='DiffLayerDiffNode')
plt.plot(SS1MSTEQ_ETA1[:,0],SS1MSTEQ_ETA1[:,2],linestyle='-',label='SamelayerSameNode1 Mstereq')
plt.plot(SS2MSTEQ_ETA1[:,0],SS2MSTEQ_ETA1[:,2],linestyle='-',label='SamelayerSameNode2 Mstereq')
plt.xlabel(r'$beta$',fontsize=14)
plt.ylabel('Correlation',fontsize=14)
plt.legend()
plt.grid()
plt.savefig('CorrBeta_ETA1.png')
plt.show()
#%%
plt.title(f'Indent space {1}')
plt.plot(SamelayerSameNode1_BETA0034[:,0],SamelayerSameNode1_BETA0034[:,2],linestyle='',marker='o',label='SamelayerSameNode1')
plt.plot(SamelayerSameNode2_BETA0034[:,0],SamelayerSameNode2_BETA0034[:,2],linestyle='',marker='o',label='SamelayerSameNode2')
plt.plot(DiffLayerSameNode_BETA0034[:,0],DiffLayerSameNode_BETA0034[:,2],linestyle='',marker='o',label='DiffLayerSameNode')
plt.plot(SameLayerDiffNode1_BETA0034[:,0],SameLayerDiffNode1_BETA0034[:,2],linestyle='',marker='o',label='SameLayerDiffNode1')
plt.plot(SamelayerDiffNode2_BETA0034[:,0],SamelayerDiffNode2_BETA0034[:,2],linestyle='',marker='o',label='SamelayerDiffNode2')
plt.plot(DiffLayerDiffNode_BETA0034[:,0],DiffLayerDiffNode_BETA0034[:,2],linestyle='',marker='o',label='DiffLayerDiffNode')
plt.plot(SS1MSTEQ_BETA0034[:,0],SS1MSTEQ_BETA0034[:,2],linestyle='-',label='SamelayerSameNode1 Mstereq')
plt.plot(SS2MSTEQ_BETA0034[:,0],SS2MSTEQ_BETA0034[:,2],linestyle='-',label='SamelayerSameNode2 Mstereq')
plt.xlabel(r'$eta$',fontsize=14)
plt.ylabel('Correlation',fontsize=14)
plt.legend()
plt.grid()
plt.savefig('CorrEta_BETA0034.png')
plt.show()
#%%
plt.title(f'Difeerence correlation between layers at indent 1')
plt.plot(SamelayerSameNode1_BETA0034[:,0],-SamelayerSameNode1_BETA0034[:,2]+SamelayerSameNode2_BETA0034[:,2],linestyle='',marker='o',label='SamelayerSameNode1')
plt.xlabel(r'$\eta$',fontsize=14)
plt.ylabel(r'$\Delta$Correlation',fontsize=14)
plt.legend()
max_diff_idx = np.argmax(-SamelayerSameNode1_BETA0034[:,2] + SamelayerSameNode2_BETA0034[:,2])
plt.axvline(SamelayerSameNode1_BETA0034[max_diff_idx,0], color='red', linestyle='--', label='Max Difference')
plt.grid()
plt.savefig('CorrDiffEta_BETA0034.png')
plt.show()
#%%
plt.title(f'lag h in different eta')
plt.plot(h, SamelayerSameNode1_ETA001[16,1:], label='SamelayerSameNode1 eta=001')
plt.plot(h, SamelayerSameNode2_ETA001[16,1:], label='SamelayerSameNode2 eta=001')
plt.plot(h, SS1MSTEQ_ETA001[24,1:], label='DiffLayerSameNode eta=001')
plt.plot(h, SS2MSTEQ_ETA001[24,1:], label='SameLayerDiffNode1 eta=001')
