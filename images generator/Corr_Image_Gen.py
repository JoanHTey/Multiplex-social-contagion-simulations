from pathlib import Path

import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT       = Path(__file__).parent          # repo root (this script's folder)
DATA_DIR   = ROOT / "Data"
IMAGES_DIR = ROOT / "Images"
IMAGES_DIR.mkdir(parents=True, exist_ok=True)


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
def groupby_mean_and_error(arr):
    arr = arr[np.argsort(arr[:, 0])]
    # Take only the first 6 digits to match the format in .txt files
    #x = np.logspace(-2, -1, 20)
    x = np.logspace(0,1.778,20) #BETA0034
    log_mids = 0.5 * (np.log10(x[1:]) + np.log10(x[:-1]))
    mids = 10**log_mids
    alls = np.sort(np.append(x, mids))
    space = np.true_divide(np.rint(alls * 10**6), 10**6)
    errors = np.zeros([len(arr[:,0]),16])
    errors[:,0] = arr[:,0]
    for i in space:
        values = arr[arr[:,0]==i]
        arr[arr[:,0]==i] = np.mean(values[:,:], axis=0)
        mask = errors[:,0]==i
        errors[mask,1:16] = np.std(values[:,1:16], axis=0)
    return arr,errors

#%%

density= np.load(str(DATA_DIR / 'density_eta1_N5000.npy'))
beta_density = np.logspace(-2,-0.3,50)

SamelayerSameNode1 = np.loadtxt(str(DATA_DIR / 'eta1_mu05_N1000/final_corr1.txt'))
SamelayerSameNode2 = np.loadtxt(str(DATA_DIR / 'eta1_mu05_N1000/final_corr2.txt'))
DiffLayerSameNode = np.loadtxt(str(DATA_DIR / 'eta1_mu05_N1000/final_corr3.txt'))
SameLayerDiffNode1 = np.loadtxt(str(DATA_DIR / 'eta1_mu05_N1000/final_corr4.txt'))
SamelayerDiffNode2 = np.loadtxt(str(DATA_DIR / 'eta1_mu05_N1000/final_corr6.txt'))
DiffLayerDiffNode = np.loadtxt(str(DATA_DIR / 'eta1_mu05_N1000/final_corr5.txt'))
SamelayerSameNode1,ESamelayerSameNode1 = groupby_mean_and_error(SamelayerSameNode1)
SamelayerSameNode2,ESamelayerSameNode2 = groupby_mean_and_error(SamelayerSameNode2)
DiffLayerSameNode,EDiffLayerSameNode = groupby_mean_and_error(DiffLayerSameNode)
SameLayerDiffNode1,ESameLayerDiffNode1 = groupby_mean_and_error(SameLayerDiffNode1)
DiffLayerDiffNode,EDiffLayerDiffNode = groupby_mean_and_error(DiffLayerDiffNode)
SamelayerDiffNode2,ESamelayerDiffNode2= groupby_mean_and_error(SamelayerDiffNode2)

SS1MSTEQ = np.loadtxt(str(DATA_DIR / 'eta1_mu05_N1000/final_corr1msteq.txt'))
SS2MSTEQ = np.loadtxt(str(DATA_DIR / 'eta1_mu05_N1000/final_corr2msteq.txt'))

den1T = np.loadtxt(str(DATA_DIR / 'eta1_mu05_N1000/den1msteq.txt'))
den2T = np.loadtxt(str(DATA_DIR / 'eta1_mu05_N1000/den2msteq.txt'))


fig, axs = plt.subplots(1, 2, figsize=(12, 5))

x,y,yerr = SamelayerSameNode1[:, 0],SamelayerSameNode1[:, 2],ESamelayerSameNode1[:, 2]
axs[0].plot(x, y, linestyle='-', marker='o', markersize=2, color='red', label=r'$\rho_{u,u}^{1,1}(1)$')
axs[0].fill_between(x, y - yerr, y + yerr, color='red', alpha=0.2)
x,y,yerr = SamelayerSameNode2[:, 0],SamelayerSameNode2[:, 2],ESamelayerSameNode2[:, 2]
axs[0].plot(x, y, linestyle='-', marker='o', markersize=2, color='green', label=r'$\rho_{u,u}^{2,2}(1)$')
axs[0].fill_between(x, y - yerr, y + yerr, color='green', alpha=0.2)
x,y,yerr = DiffLayerSameNode[:,0],DiffLayerSameNode[:,2],EDiffLayerSameNode[:,2]
axs[0].plot(x, y, linestyle='-', marker='o', markersize=2, color='purple', label=r'$\rho_{u,u}^{1,2}(1)$')
axs[0].fill_between(x, y - yerr, y + yerr, color='purple', alpha=0.2)
x,y,yerr = SameLayerDiffNode1[:,0],SameLayerDiffNode1[:,2],ESameLayerDiffNode1[:,2]
axs[0].plot(x, y, linestyle='-', marker='o', markersize=2, color='cyan', label=r'$\rho_{u,v}^{1,1}(1)$')
axs[0].fill_between(x, y - yerr, y + yerr, color='cyan', alpha=0.2)
x,y,yerr = SamelayerDiffNode2[:,0],SamelayerDiffNode2[:,2],ESamelayerDiffNode2[:,2]
axs[0].plot(x, y, linestyle='-', marker='o', markersize=2, color='orange', label=r'$\rho_{u,v}^{2,2}(1)$')
axs[0].fill_between(x, y - yerr, y + yerr, color='orange', alpha=0.2)
x,y,yerr = DiffLayerDiffNode[:,0],DiffLayerDiffNode[:,2],EDiffLayerDiffNode[:,2]
axs[0].plot(x, y, linestyle='-', marker='o', markersize=2, color='blue', label=r'$\rho_{u,v}^{1,2}(1)$')
axs[0].fill_between(x, y - yerr, y + yerr, color='blue', alpha=0.2)

mask = SS1MSTEQ[:,0]>1/30
axs[0].plot(SS1MSTEQ[mask,0],SS1MSTEQ[mask,2],linestyle='--',color='darkred',label=r'$\rho_{u,u}^{1,1}(1)$ T.')
mask = SS2MSTEQ[:,0]>1/30
axs[0].plot(SS2MSTEQ[mask,0],SS2MSTEQ[mask,2],linestyle='--',color='darkgreen',label=r'$\rho_{u,u}^{2,2}(1)$ T.')
axs[0].set_xlabel(r'$\beta/\mu$',fontsize=14)
axs[0].set_ylabel(r'$\rho(1)$',fontsize=14)
axs[0].set_xscale('log')
axs[0].axvline(1/30,linestyle='--',color='red',label = r'$\frac{1}{\Lambda_1}$')
axs[0].axvline(1/10,linestyle='--',color='green',label = r'$\frac{1}{\Lambda_2}$')
axs[0].legend(fontsize=12)

axs[1].plot(beta_density,density[:,0],marker='o',markersize='2',color='red',label=r'$\phi_1$')
axs[1].plot(beta_density,density[:,1],marker='o',markersize='2',color='green',label=r'$\phi_2$')
axs[1].plot(den1T[:,0],den1T[:,1],linestyle='--',color='darkred',label=r'$\phi_1$ T.')
axs[1].plot(den2T[:,0],den2T[:,1],linestyle='--',color='darkgreen',label=r'$\phi_2$ T.')
axs[1].set_xlabel(r'$\beta/\mu$',fontsize=14)
axs[1].set_ylabel(r'$\phi$',fontsize=14)
axs[1].set_xscale('log')
axs[1].axvline(1/30,linestyle='--',color='red',label = r'$\frac{1}{\Lambda_1}$')
axs[1].axvline(1/10,linestyle='--',color='green',label = r'$\frac{1}{\Lambda_2}$')
axs[1].legend(fontsize=14)

plt.tight_layout()
plt.savefig(str(IMAGES_DIR / 'CorrBeta_ETA1.pdf'))
plt.show()


#%%


# # Load the data from the text files
SamelayerSameNode1 = np.loadtxt(str(DATA_DIR / 'beta0034_mu05_N1000/final_corr1.txt'))
SamelayerSameNode2 = np.loadtxt(str(DATA_DIR / 'beta0034_mu05_N1000/final_corr2.txt'))
DiffLayerSameNode = np.loadtxt(str(DATA_DIR / 'beta0034_mu05_N1000/final_corr3.txt'))
SameLayerDiffNode1 = np.loadtxt(str(DATA_DIR / 'beta0034_mu05_N1000/final_corr4.txt'))
DiffLayerDiffNode = np.loadtxt(str(DATA_DIR / 'beta0034_mu05_N1000/final_corr5.txt'))
SamelayerDiffNode2= np.loadtxt(str(DATA_DIR / 'beta0034_mu05_N1000/final_corr6.txt'))
SamelayerSameNode1,ESamelayerSameNode1 = groupby_mean_and_error(SamelayerSameNode1)
SamelayerSameNode2,ESamelayerSameNode2 = groupby_mean_and_error(SamelayerSameNode2)
DiffLayerSameNode,EDiffLayerSameNode = groupby_mean_and_error(DiffLayerSameNode)
SameLayerDiffNode1,ESameLayerDiffNode1 = groupby_mean_and_error(SameLayerDiffNode1)
DiffLayerDiffNode,EDiffLayerDiffNode = groupby_mean_and_error(DiffLayerDiffNode)
SamelayerDiffNode2,ESamelayerDiffNode2= groupby_mean_and_error(SamelayerDiffNode2)


SS1MSTEQ = np.loadtxt(str(DATA_DIR / 'beta0034_mu05_N1000/final_corr1msteq.txt'))
SS2MSTEQ = np.loadtxt(str(DATA_DIR / 'beta0034_mu05_N1000/final_corr2msteq.txt'))


fig, axs = plt.subplots(1, 1, figsize=(6, 5))


x,y,yerr = SamelayerSameNode1[:, 0],SamelayerSameNode1[:, 2],ESamelayerSameNode1[:, 2]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='red', label=r'$\rho_{u,u}^{1,1}(1)$')
axs.fill_between(x, y - yerr, y + yerr, color='red', alpha=0.2)
x,y,yerr = SamelayerSameNode2[:, 0],SamelayerSameNode2[:, 2],ESamelayerSameNode2[:, 2]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='green', label=r'$\rho_{u,u}^{2,2}(1)$')
axs.fill_between(x, y - yerr, y + yerr, color='green', alpha=0.2)
x,y,yerr = DiffLayerSameNode[:,0],DiffLayerSameNode[:,2],EDiffLayerSameNode[:,2]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='purple', label=r'$\rho_{u,u}^{1,2}(1)$')
axs.fill_between(x, y - yerr, y + yerr, color='purple', alpha=0.2)
x,y,yerr = SameLayerDiffNode1[:,0],SameLayerDiffNode1[:,2],ESameLayerDiffNode1[:,2]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='cyan', label=r'$\rho_{u,v}^{1,1}(1)$')
axs.fill_between(x, y - yerr, y + yerr, color='cyan', alpha=0.2)
x,y,yerr = SamelayerDiffNode2[:,0],SamelayerDiffNode2[:,2],ESamelayerDiffNode2[:,2]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='orange', label=r'$\rho_{u,v}^{2,2}(1)$')
axs.fill_between(x, y - yerr, y + yerr, color='orange', alpha=0.2)
x,y,yerr = DiffLayerDiffNode[:,0],DiffLayerDiffNode[:,2],EDiffLayerDiffNode[:,2]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='blue', label=r'$\rho_{u,v}^{1,2}(1)$')
axs.fill_between(x, y - yerr, y + yerr, color='blue', alpha=0.2)


axs.plot(SS1MSTEQ[:,0],SS1MSTEQ[:,2],linestyle='--',color='darkred',label=r'$\rho_{u,u}^{1,1}(1)$ T.')
axs.plot(SS2MSTEQ[:,0],SS2MSTEQ[:,2],linestyle='--',color='darkgreen',label=r'$\rho_{u,u}^{2,2}(1)$ T.')
axs.set_xlabel(r'$\eta/\beta$',fontsize=14)
axs.set_ylabel(r'$\rho(1)$',fontsize=14)
axs.set_xscale('log')
axs.axvline(21.1271827928,linestyle='--',color='black', label=r'$\left(\frac{\eta}{\beta}\right)_{\mathrm{crit}}$')
axs.legend(fontsize=12,loc='upper left')

plt.tight_layout()
plt.savefig(str(IMAGES_DIR / 'CorrBeta_BETA0034.pdf'))
plt.show()

#%%
# Extract data and errors
x = SamelayerSameNode1[:, 0]  # eta/beta
y1 = SamelayerSameNode1[:, 2]
y2 = SamelayerSameNode2[:, 2]
yerr1 = ESamelayerSameNode1[:, 2]
yerr2 = ESamelayerSameNode2[:, 2]

# Compute delta correlation and its propagated error
delta_rho = -y1 + y2
delta_rho_err = np.sqrt(yerr1**2 + yerr2**2)

# Plot
plt.figure(figsize=(6, 5))
plt.plot(x, delta_rho, linestyle='-', marker='o', markersize=4, label=r'$\Delta\rho(1)$')
plt.fill_between(x, delta_rho - delta_rho_err, delta_rho + delta_rho_err, alpha=0.2)

# Labels and formatting
plt.xlabel(r'$\eta/\beta$', fontsize=14)
plt.ylabel(r'$\Delta\rho(1)$', fontsize=14)
plt.xscale('log')

# Add vertical critical line
plt.axvline(21.1271827928, color='black', linestyle='--', label=r'$\left(\frac{\eta}{\beta}\right)_{\mathrm{crit}}$')

# Legend and save
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig(str(IMAGES_DIR / 'CorrDiffEta_BETA0034.pdf'))
plt.show()


#%%
density= np.load(str(DATA_DIR / 'density_eta001_N5000.npy'))
beta_density = np.logspace(-2,-0.3,50)
h = np.linspace(0, 14, 15,dtype=int)

SamelayerSameNode1 = np.loadtxt(str(DATA_DIR / 'eta001_mu05_N1000/final_corr1.txt'))
SamelayerSameNode2 = np.loadtxt(str(DATA_DIR / 'eta001_mu05_N1000/final_corr2.txt'))
DiffLayerSameNode = np.loadtxt(str(DATA_DIR / 'eta001_mu05_N1000/final_corr3.txt'))
SameLayerDiffNode1 = np.loadtxt(str(DATA_DIR / 'eta001_mu05_N1000/final_corr4.txt'))
DiffLayerDiffNode = np.loadtxt(str(DATA_DIR / 'eta001_mu05_N1000/final_corr5.txt'))
SamelayerDiffNode2 = np.loadtxt(str(DATA_DIR / 'eta001_mu05_N1000/final_corr6.txt'))
SamelayerSameNode1,ESamelayerSameNode1 = groupby_mean_and_error(SamelayerSameNode1)
SamelayerSameNode2,ESamelayerSameNode2 = groupby_mean_and_error(SamelayerSameNode2)
DiffLayerSameNode,EDiffLayerSameNode = groupby_mean_and_error(DiffLayerSameNode)
SameLayerDiffNode1,ESameLayerDiffNode1 = groupby_mean_and_error(SameLayerDiffNode1)
DiffLayerDiffNode,EDiffLayerDiffNode = groupby_mean_and_error(DiffLayerDiffNode)
SamelayerDiffNode2,ESamelayerDiffNode2= groupby_mean_and_error(SamelayerDiffNode2)


SS1MSTEQ = np.loadtxt(str(DATA_DIR / 'eta001_mu05_N1000/final_corr1msteq.txt'))
SS2MSTEQ = np.loadtxt(str(DATA_DIR / 'eta001_mu05_N1000/final_corr2msteq.txt'))

den1T = np.loadtxt(str(DATA_DIR / 'eta001_mu05_N1000/den1msteq.txt'))
den2T = np.loadtxt(str(DATA_DIR / 'eta001_mu05_N1000/den2msteq.txt'))

fig, axs = plt.subplots(1, 2, figsize=(12, 5))

x,y,yerr = SamelayerSameNode1[:, 0],SamelayerSameNode1[:, 2],ESamelayerSameNode1[:, 2]
axs[0].plot(x, y, linestyle='-', marker='o', markersize=2, color='red', label=r'$\rho_{u,u}^{1,1}(1)$')
axs[0].fill_between(x, y - yerr, y + yerr, color='red', alpha=0.2)
x,y,yerr = SamelayerSameNode2[:, 0],SamelayerSameNode2[:, 2],ESamelayerSameNode2[:, 2]
axs[0].plot(x, y, linestyle='-', marker='o', markersize=2, color='green', label=r'$\rho_{u,u}^{2,2}(1)$')
axs[0].fill_between(x, y - yerr, y + yerr, color='green', alpha=0.2)
x,y,yerr = DiffLayerSameNode[:,0],DiffLayerSameNode[:,2],EDiffLayerSameNode[:,2]
axs[0].plot(x, y, linestyle='-', marker='o', markersize=2, color='purple', label=r'$\rho_{u,u}^{1,2}(1)$')
axs[0].fill_between(x, y - yerr, y + yerr, color='purple', alpha=0.2)
x,y,yerr = SameLayerDiffNode1[:,0],SameLayerDiffNode1[:,2],ESameLayerDiffNode1[:,2]
axs[0].plot(x, y, linestyle='-', marker='o', markersize=2, color='cyan', label=r'$\rho_{u,v}^{1,1}(1)$')
axs[0].fill_between(x, y - yerr, y + yerr, color='cyan', alpha=0.2)
x,y,yerr = SamelayerDiffNode2[:,0],SamelayerDiffNode2[:,2],ESamelayerDiffNode2[:,2]
axs[0].plot(x, y, linestyle='-', marker='o', markersize=2, color='orange', label=r'$\rho_{u,v}^{2,2}(1)$')
axs[0].fill_between(x, y - yerr, y + yerr, color='orange', alpha=0.2)
x,y,yerr = DiffLayerDiffNode[:,0],DiffLayerDiffNode[:,2],EDiffLayerDiffNode[:,2]
axs[0].plot(x, y, linestyle='-', marker='o', markersize=2, color='blue', label=r'$\rho_{u,v}^{1,2}(1)$')
axs[0].fill_between(x, y - yerr, y + yerr, color='blue', alpha=0.2)

mask = SS1MSTEQ[:,0]>1/30
axs[0].plot(SS1MSTEQ[mask,0],SS1MSTEQ[mask,2],linestyle='--',color='darkred',label=r'$\rho_{u,u}^{1,1}(1)$ T.')
mask = SS2MSTEQ[:,0]>1/30
axs[0].plot(SS2MSTEQ[mask,0],SS2MSTEQ[mask,2],linestyle='--',color='darkgreen',label=r'$\rho_{u,u}^{2,2}(1)$ T.')
axs[0].set_xlabel(r'$\beta/\mu$',fontsize=14)
axs[0].set_ylabel(r'$\rho(1)$',fontsize=14)
axs[0].set_xscale('log')
axs[0].axvline(1/30,linestyle='--',color='red',label = r'$\frac{1}{\Lambda_1}$')
axs[0].axvline(1/10,linestyle='--',color='green',label = r'$\frac{1}{\Lambda_2}$')
axs[0].legend(fontsize=12)

axs[1].plot(beta_density,density[:,0],marker='o',markersize='2',color='red',label=r'$\phi_1$')
axs[1].plot(beta_density,density[:,1],marker='o',markersize='2',color='green',label=r'$\phi_2$')
axs[1].plot(den1T[:,0],den1T[:,1],linestyle='--',color='darkred',label=r'$\phi_1$ T.')
axs[1].plot(den2T[:,0],den2T[:,1],linestyle='--',color='darkgreen',label=r'$\phi_2$ T.')
axs[1].set_xlabel(r'$\beta/\mu$',fontsize=14)
axs[1].set_ylabel(r'$\phi$',fontsize=14)
axs[1].set_xscale('log')
axs[1].axvline(1/30,linestyle='--',color='red',label = r'$\frac{1}{\Lambda_1}$')
axs[1].axvline(1/10,linestyle='--',color='green',label = r'$\frac{1}{\Lambda_2}$')
axs[1].legend(fontsize=14)

plt.tight_layout()
plt.savefig(str(IMAGES_DIR / 'CorrBeta_ETA001.pdf'))
plt.show()
#%%
plt.figure(figsize=(6, 5))

mask = 230

plt.plot(h[:7], SamelayerSameNode1[mask, 1:8], linestyle='-', marker='o', markersize=2, color='red', label=r'$\rho_{u,u}^{1,1}(h)$')
plt.plot(h[:7], SamelayerSameNode2[mask, 1:8], linestyle='-', marker='o', markersize=2, color='green', label=r'$\rho_{u,u}^{2,2}(h)$')
plt.plot(h[:7], DiffLayerSameNode[mask, 1:8], linestyle='-', marker='o', markersize=2, color='purple', label=r'$\rho_{u,u}^{1,2}(h)$')
plt.plot(h[:7], SameLayerDiffNode1[mask, 1:8], linestyle='-', marker='o', markersize=2, color='cyan', label=r'$\rho_{u,v}^{1,1}(h)$')
plt.plot(h[:7], SamelayerDiffNode2[mask, 1:8], linestyle='-', marker='o', markersize=2, color='orange', label=r'$\rho_{u,v}^{2,2}(h)$')
plt.plot(h[:7], DiffLayerDiffNode[mask, 1:8], linestyle='-', marker='o', markersize=2, color='blue', label=r'$\rho_{u,v}^{1,2}(h)$')
plt.plot(h[:7],SS1MSTEQ[49,1:8],linestyle='--',color='darkred',label=r'$\rho_{u,u}^{1,1}(1)$ T.' )
plt.plot(h[:7],SS2MSTEQ[49,1:8],linestyle='--',color='darkgreen',label=r'$\rho_{u,u}^{2,2}(1)$ T.' )
plt.xlabel('lag h',fontsize=14)
plt.ylabel(r'$\rho(h)$',fontsize=14)
plt.legend(fontsize=14)
plt.savefig(str(IMAGES_DIR / 'corrcorr.pdf'))
plt.show()

#%%
density= np.load(str(DATA_DIR / 'density_eta25_N5000.npy'))
beta_density = np.logspace(-2,-0.3,50)
h = np.linspace(0, 14, 15,dtype=int)

den1T = np.loadtxt(str(DATA_DIR / 'eta25_mu05_N1000/den1msteq.txt'))
den2T = np.loadtxt(str(DATA_DIR / 'eta25_mu05_N1000/den2msteq.txt'))

SamelayerSameNode1 = np.loadtxt(str(DATA_DIR / 'eta25_mu05_N1000/final_corr1.txt'))
SamelayerSameNode2 = np.loadtxt(str(DATA_DIR / 'eta25_mu05_N1000/final_corr2.txt'))
DiffLayerSameNode = np.loadtxt(str(DATA_DIR / 'eta25_mu05_N1000/final_corr3.txt'))
SameLayerDiffNode1 = np.loadtxt(str(DATA_DIR / 'eta25_mu05_N1000/final_corr4.txt'))
DiffLayerDiffNode = np.loadtxt(str(DATA_DIR / 'eta25_mu05_N1000/final_corr5.txt'))
SamelayerDiffNode2 = np.loadtxt(str(DATA_DIR / 'eta25_mu05_N1000/final_corr6.txt'))
SamelayerSameNode1,ESamelayerSameNode1 = groupby_mean_and_error(SamelayerSameNode1)
SamelayerSameNode2,ESamelayerSameNode2 = groupby_mean_and_error(SamelayerSameNode2)
DiffLayerSameNode,EDiffLayerSameNode = groupby_mean_and_error(DiffLayerSameNode)
SameLayerDiffNode1,ESameLayerDiffNode1 = groupby_mean_and_error(SameLayerDiffNode1)
DiffLayerDiffNode,EDiffLayerDiffNode = groupby_mean_and_error(DiffLayerDiffNode)
SamelayerDiffNode2,ESamelayerDiffNode2= groupby_mean_and_error(SamelayerDiffNode2)


SS1MSTEQ = np.loadtxt(str(DATA_DIR / 'eta25_mu05_N1000/final_corr1msteq.txt'))
SS2MSTEQ = np.loadtxt(str(DATA_DIR / 'eta25_mu05_N1000/final_corr2msteq.txt'))


fig, axs = plt.subplots(1, 2, figsize=(12, 5))

x,y,yerr = SamelayerSameNode1[:, 0],SamelayerSameNode1[:, 2],ESamelayerSameNode1[:, 2]
axs[0].plot(x, y, linestyle='-', marker='o', markersize=2, color='red', label=r'$\rho_{u,u}^{1,1}(1)$')
axs[0].fill_between(x, y - yerr, y + yerr, color='red', alpha=0.2)
x,y,yerr = SamelayerSameNode2[:, 0],SamelayerSameNode2[:, 2],ESamelayerSameNode2[:, 2]
axs[0].plot(x, y, linestyle='-', marker='o', markersize=2, color='green', label=r'$\rho_{u,u}^{2,2}(1)$')
axs[0].fill_between(x, y - yerr, y + yerr, color='green', alpha=0.2)
x,y,yerr = DiffLayerSameNode[:,0],DiffLayerSameNode[:,2],EDiffLayerSameNode[:,2]
axs[0].plot(x, y, linestyle='-', marker='o', markersize=2, color='purple', label=r'$\rho_{u,u}^{1,2}(1)$')
axs[0].fill_between(x, y - yerr, y + yerr, color='purple', alpha=0.2)
x,y,yerr = SameLayerDiffNode1[:,0],SameLayerDiffNode1[:,2],ESameLayerDiffNode1[:,2]
axs[0].plot(x, y, linestyle='-', marker='o', markersize=2, color='cyan', label=r'$\rho_{u,v}^{1,1}(1)$')
axs[0].fill_between(x, y - yerr, y + yerr, color='cyan', alpha=0.2)
x,y,yerr = SamelayerDiffNode2[:,0],SamelayerDiffNode2[:,2],ESamelayerDiffNode2[:,2]
axs[0].plot(x, y, linestyle='-', marker='o', markersize=2, color='orange', label=r'$\rho_{u,v}^{2,2}(1)$')
axs[0].fill_between(x, y - yerr, y + yerr, color='orange', alpha=0.2)
x,y,yerr = DiffLayerDiffNode[:,0],DiffLayerDiffNode[:,2],EDiffLayerDiffNode[:,2]
axs[0].plot(x, y, linestyle='-', marker='o', markersize=2, color='blue', label=r'$\rho_{u,v}^{1,2}(1)$')
axs[0].fill_between(x, y - yerr, y + yerr, color='blue', alpha=0.2)

mask = SS1MSTEQ[:,0]>1/50
axs[0].plot(SS1MSTEQ[mask,0],SS1MSTEQ[mask,2],linestyle='--',color='darkred',label=r'$\rho_{u,u}^{1,1}(1)$ T.')
mask = SS2MSTEQ[:,0]>1/50
axs[0].plot(SS2MSTEQ[mask,0],SS2MSTEQ[mask,2],linestyle='--',color='darkgreen',label=r'$\rho_{u,u}^{2,2}(1)$ T.')
axs[0].set_xlabel(r'$\beta/\mu$',fontsize=14)
axs[0].set_ylabel(r'$\rho(1)$',fontsize=14)
axs[0].set_xscale('log')
axs[0].axvline(1/46.93,linestyle='--',color='red',label = r'$\frac{1}{\Lambda_T}$')
axs[0].axvline(1/43,linestyle='--',color='green',label = r'$\frac{1}{\Lambda_S}$')
axs[0].legend(fontsize=12)

axs[1].plot(beta_density[:27],density[:27,0],marker='o',markersize='2',color='red',label=r'$\phi_1$')
axs[1].plot(beta_density[:27],density[:27,1],marker='o',markersize='2',color='green',label=r'$\phi_2$')
axs[1].plot(den1T[:,0],den1T[:,1],linestyle='--',color='darkred',label=r'$\phi_1$ T.')
axs[1].plot(den2T[:,0],den2T[:,1],linestyle='--',color='darkgreen',label=r'$\phi_2$ T.')
axs[1].set_xlabel(r'$\beta/\mu$',fontsize=14)
axs[1].set_ylabel(r'$\phi$',fontsize=14)
axs[1].set_xscale('log')
axs[1].axvline(1/46.93,linestyle='--',color='red',label = r'$\frac{1}{\Lambda_T}$')
axs[1].axvline(1/43,linestyle='--',color='green',label = r'$\frac{1}{\Lambda_S}$')
axs[1].legend(fontsize=14)

plt.tight_layout()
plt.savefig(str(IMAGES_DIR / 'CorrBeta_ETA25.pdf'))
plt.show()



#%%
h = np.linspace(0, 14, 15,dtype=int)

mu1lay1 = np.loadtxt(str(DATA_DIR / 'Muanalysis/mu1/final_corr1.txt'))
Emu1lay1 = np.std(mu1lay1,axis=0)
mu1lay1 = np.mean(mu1lay1,axis=0)
mu1lay2 = np.loadtxt(str(DATA_DIR / 'Muanalysis/mu1/final_corr2.txt'))
Emu1lay2 = np.std(mu1lay2,axis=0)
mu1lay2 = np.mean(mu1lay2,axis=0)
mu1lay12 = np.loadtxt(str(DATA_DIR / 'Muanalysis/mu1/final_corr3.txt'))
Emu1lay12 = np.std(mu1lay12,axis=0)
mu1lay12 = np.mean(mu1lay12,axis=0)
mu1lay11 = np.loadtxt(str(DATA_DIR / 'Muanalysis/mu1/final_corr4.txt'))
Emu1lay11 = np.std(mu1lay11,axis=0)
mu1lay11 = np.mean(mu1lay11,axis=0)
mu1lay22 = np.loadtxt(str(DATA_DIR / 'Muanalysis/mu1/final_corr5.txt'))
Emu1lay22 = np.std(mu1lay22,axis=0)
mu1lay22 = np.mean(mu1lay22,axis=0)
mu1lay122 = np.loadtxt(str(DATA_DIR / 'Muanalysis/mu1/final_corr6.txt'))
Emu1lay122 = np.std(mu1lay122,axis=0)
mu1lay122 = np.mean(mu1lay122,axis=0)

mu05lay1 = np.loadtxt(str(DATA_DIR / 'Muanalysis/mu05/final_corr1.txt'))
Emu05lay1 = np.std(mu05lay1,axis=0)
mu05lay1 = np.mean(mu05lay1,axis=0)
mu05lay2 = np.loadtxt(str(DATA_DIR / 'Muanalysis/mu05/final_corr2.txt'))
Emu05lay2 = np.std(mu05lay2,axis=0)
mu05lay2 = np.mean(mu05lay2,axis=0)
mu05lay12 = np.loadtxt(str(DATA_DIR / 'Muanalysis/mu05/final_corr3.txt'))
Emu05lay12 = np.std(mu05lay12,axis=0)
mu05lay12 = np.mean(mu05lay12,axis=0)
mu05lay11 = np.loadtxt(str(DATA_DIR / 'Muanalysis/mu05/final_corr4.txt'))
Emu05lay11 = np.std(mu05lay11,axis=0)
mu05lay11 = np.mean(mu05lay11,axis=0)
mu05lay22 = np.loadtxt(str(DATA_DIR / 'Muanalysis/mu05/final_corr5.txt'))
Emu05lay22 = np.std(mu05lay22,axis=0)
mu05lay22 = np.mean(mu05lay22,axis=0)
mu05lay122 = np.loadtxt(str(DATA_DIR / 'Muanalysis/mu05/final_corr6.txt'))
Emu05lay122 = np.std(mu05lay122,axis=0)
mu05lay122 = np.mean(mu05lay122,axis=0)

mu075lay1 = np.loadtxt(str(DATA_DIR / 'Muanalysis/mu075/final_corr1.txt'))
Emu075lay1 = np.std(mu075lay1,axis=0)
mu075lay1 = np.mean(mu075lay1,axis=0)
mu075lay2 = np.loadtxt(str(DATA_DIR / 'Muanalysis/mu075/final_corr2.txt'))
Emu075lay2 = np.std(mu075lay2,axis=0)
mu075lay2 = np.mean(mu075lay2,axis=0)
mu075lay12 = np.loadtxt(str(DATA_DIR / 'Muanalysis/mu075/final_corr3.txt'))
Emu075lay12 = np.std(mu075lay12,axis=0)
mu075lay12 = np.mean(mu075lay12,axis=0)
mu075lay11 = np.loadtxt(str(DATA_DIR / 'Muanalysis/mu075/final_corr4.txt'))
Emu075lay11 = np.std(mu075lay11,axis=0)
mu075lay11 = np.mean(mu075lay11,axis=0)
mu075lay22 = np.loadtxt(str(DATA_DIR / 'Muanalysis/mu075/final_corr5.txt'))
Emu075lay22 = np.std(mu075lay22,axis=0)
mu075lay22 = np.mean(mu075lay22,axis=0)
mu075lay122 = np.loadtxt(str(DATA_DIR / 'Muanalysis/mu075/final_corr6.txt'))
Emu075lay122 = np.std(mu075lay122,axis=0)
mu075lay122 = np.mean(mu075lay122,axis=0)

mu095lay1 = np.loadtxt(str(DATA_DIR / 'Muanalysis/mu095/final_corr1.txt'))
Emu095lay1 = np.std(mu095lay1,axis=0)
mu095lay1 = np.mean(mu095lay1,axis=0)
mu095lay2 = np.loadtxt(str(DATA_DIR / 'Muanalysis/mu095/final_corr2.txt'))
Emu095lay2 = np.std(mu095lay2,axis=0)
mu095lay2 = np.mean(mu095lay2,axis=0)
mu095lay12 = np.loadtxt(str(DATA_DIR / 'Muanalysis/mu095/final_corr3.txt'))
Emu095lay12 = np.std(mu095lay12,axis=0)
mu095lay12 = np.mean(mu095lay12,axis=0)
mu095lay11 = np.loadtxt(str(DATA_DIR / 'Muanalysis/mu095/final_corr4.txt'))
Emu095lay11 = np.std(mu095lay11,axis=0)
mu095lay11 = np.mean(mu095lay11,axis=0)
mu095lay22 = np.loadtxt(str(DATA_DIR / 'Muanalysis/mu095/final_corr5.txt'))
Emu095lay22 = np.std(mu095lay22,axis=0)
mu095lay22 = np.mean(mu095lay22,axis=0)
mu095lay122 = np.loadtxt(str(DATA_DIR / 'Muanalysis/mu095/final_corr6.txt'))
Emu095lay122 = np.std(mu095lay122,axis=0)
mu095lay122 = np.mean(mu095lay122,axis=0)

fig, axs = plt.subplots(1, 1, figsize=(6, 5))

x,y,yerr = h,mu095lay1[1:,],Emu095lay1[1:,]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='red', label=r'$\mu=0.95$')
axs.fill_between(x, y - yerr, y + yerr, color='red', alpha=0.2)

x,y,yerr = h,mu1lay1[1:,],Emu1lay1[1:,]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='blue', label=r'$\mu=1$')
axs.fill_between(x, y - yerr, y + yerr, color='red', alpha=0.2)

x,y,yerr = h,mu05lay1[1:,],Emu05lay1[1:,]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='green', label=r'$\mu=0.5$')
axs.fill_between(x, y - yerr, y + yerr, color='red', alpha=0.2)

x,y,yerr = h,mu075lay1[1:,],Emu075lay1[1:,]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='yellow', label=r'$\mu=0.75$')
axs.fill_between(x, y - yerr, y + yerr, color='red', alpha=0.2)

axs.set_xlabel(r'lag $h$',fontsize=14)
axs.set_ylabel(r'$\rho(h)$',fontsize=14)
axs.legend(fontsize=14)

plt.tight_layout()
plt.savefig(str(IMAGES_DIR / 'muAnalisis.pdf'))
plt.show()

#%%

SamelayerSameNode1 = np.loadtxt(str(DATA_DIR / 'EReta1_mu05/final_corr1.txt'))
SamelayerSameNode2 = np.loadtxt(str(DATA_DIR / 'EReta1_mu05/final_corr2.txt'))
DiffLayerSameNode = np.loadtxt(str(DATA_DIR / 'EReta1_mu05/final_corr3.txt'))
SameLayerDiffNode1 = np.loadtxt(str(DATA_DIR / 'EReta1_mu05/final_corr4.txt'))
DiffLayerDiffNode = np.loadtxt(str(DATA_DIR / 'EReta1_mu05/final_corr5.txt'))
SamelayerDiffNode2 = np.loadtxt(str(DATA_DIR / 'EReta1_mu05/final_corr6.txt'))
SamelayerSameNode1,ESamelayerSameNode1 = groupby_mean_and_error(SamelayerSameNode1)
SamelayerSameNode2,ESamelayerSameNode2 = groupby_mean_and_error(SamelayerSameNode2)
DiffLayerSameNode,EDiffLayerSameNode = groupby_mean_and_error(DiffLayerSameNode)
SameLayerDiffNode1,ESameLayerDiffNode1 = groupby_mean_and_error(SameLayerDiffNode1)
DiffLayerDiffNode,EDiffLayerDiffNode = groupby_mean_and_error(DiffLayerDiffNode)
SamelayerDiffNode2,ESamelayerDiffNode2= groupby_mean_and_error(SamelayerDiffNode2)


SS1MSTEQ = np.loadtxt(str(DATA_DIR / 'eta1_mu05_N1000/final_corr1msteq.txt'))
SS2MSTEQ = np.loadtxt(str(DATA_DIR / 'eta1_mu05_N1000/final_corr2msteq.txt'))


fig, axs = plt.subplots(1, 1, figsize=(6, 5))

x,y,yerr = SamelayerSameNode1[:, 0],SamelayerSameNode1[:, 2],ESamelayerSameNode1[:, 2]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='red', label=r'$\rho_{u,u}^{1,1}(1)$')
axs.fill_between(x, y - yerr, y + yerr, color='red', alpha=0.2)
x,y,yerr = SamelayerSameNode2[:, 0],SamelayerSameNode2[:, 2],ESamelayerSameNode2[:, 2]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='green', label=r'$\rho_{u,u}^{2,2}(1)$')
axs.fill_between(x, y - yerr, y + yerr, color='green', alpha=0.2)
x,y,yerr = DiffLayerSameNode[:,0],DiffLayerSameNode[:,2],EDiffLayerSameNode[:,2]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='purple', label=r'$\rho_{u,u}^{1,2}(1)$')
axs.fill_between(x, y - yerr, y + yerr, color='purple', alpha=0.2)
x,y,yerr = SameLayerDiffNode1[:,0],SameLayerDiffNode1[:,2],ESameLayerDiffNode1[:,2]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='cyan', label=r'$\rho_{u,v}^{1,1}(1)$')
axs.fill_between(x, y - yerr, y + yerr, color='cyan', alpha=0.2)
x,y,yerr = SamelayerDiffNode2[:,0],SamelayerDiffNode2[:,2],ESamelayerDiffNode2[:,2]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='orange', label=r'$\rho_{u,v}^{2,2}(1)$')
axs.fill_between(x, y - yerr, y + yerr, color='orange', alpha=0.2)
x,y,yerr = DiffLayerDiffNode[:,0],DiffLayerDiffNode[:,2],EDiffLayerDiffNode[:,2]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='blue', label=r'$\rho_{u,v}^{1,2}(1)$')
axs.fill_between(x, y - yerr, y + yerr, color='blue', alpha=0.2)

axs.axvline(1/31,linestyle='--',color='red',label = r'$\frac{1}{\Lambda_1}$')
axs.axvline(1/11,linestyle='--',color='green',label = r'$\frac{1}{\Lambda_2}$')
axs.set_xlabel(r'$\eta/\beta$',fontsize=14)
axs.set_ylabel(r'$\rho(1)$',fontsize=14)
axs.set_xscale('log')
axs.legend(fontsize=12)

plt.tight_layout()
plt.savefig(str(IMAGES_DIR / 'EReta1.pdf'))
plt.show()

#%%
SamelayerSameNode1 = np.loadtxt(str(DATA_DIR / 'EReta001_mu05/final_corr1.txt'))
SamelayerSameNode2 = np.loadtxt(str(DATA_DIR / 'EReta001_mu05/final_corr2.txt'))
DiffLayerSameNode = np.loadtxt(str(DATA_DIR / 'EReta001_mu05/final_corr3.txt'))
SameLayerDiffNode1 = np.loadtxt(str(DATA_DIR / 'EReta001_mu05/final_corr4.txt'))
DiffLayerDiffNode = np.loadtxt(str(DATA_DIR / 'EReta001_mu05/final_corr5.txt'))
SamelayerDiffNode2 = np.loadtxt(str(DATA_DIR / 'EReta001_mu05/final_corr6.txt'))
SamelayerSameNode1,ESamelayerSameNode1 = groupby_mean_and_error(SamelayerSameNode1)
SamelayerSameNode2,ESamelayerSameNode2 = groupby_mean_and_error(SamelayerSameNode2)
DiffLayerSameNode,EDiffLayerSameNode = groupby_mean_and_error(DiffLayerSameNode)
SameLayerDiffNode1,ESameLayerDiffNode1 = groupby_mean_and_error(SameLayerDiffNode1)
DiffLayerDiffNode,EDiffLayerDiffNode = groupby_mean_and_error(DiffLayerDiffNode)
SamelayerDiffNode2,ESamelayerDiffNode2= groupby_mean_and_error(SamelayerDiffNode2)


SS1MSTEQ = np.loadtxt(str(DATA_DIR / 'eta001_mu05_N1000/final_corr1msteq.txt'))
SS2MSTEQ = np.loadtxt(str(DATA_DIR / 'eta001_mu05_N1000/final_corr2msteq.txt'))


fig, axs = plt.subplots(1, 1, figsize=(6, 5))

x,y,yerr = SamelayerSameNode1[:, 0],SamelayerSameNode1[:, 2],ESamelayerSameNode1[:, 2]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='red', label=r'$\rho_{u,u}^{1,1}(1)$')
axs.fill_between(x, y - yerr, y + yerr, color='red', alpha=0.2)
x,y,yerr = SamelayerSameNode2[:, 0],SamelayerSameNode2[:, 2],ESamelayerSameNode2[:, 2]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='green', label=r'$\rho_{u,u}^{2,2}(1)$')
axs.fill_between(x, y - yerr, y + yerr, color='green', alpha=0.2)
x,y,yerr = DiffLayerSameNode[:,0],DiffLayerSameNode[:,2],EDiffLayerSameNode[:,2]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='purple', label=r'$\rho_{u,u}^{1,2}(1)$')
axs.fill_between(x, y - yerr, y + yerr, color='purple', alpha=0.2)
x,y,yerr = SameLayerDiffNode1[:,0],SameLayerDiffNode1[:,2],ESameLayerDiffNode1[:,2]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='cyan', label=r'$\rho_{u,v}^{1,1}(1)$')
axs.fill_between(x, y - yerr, y + yerr, color='cyan', alpha=0.2)
x,y,yerr = SamelayerDiffNode2[:,0],SamelayerDiffNode2[:,2],ESamelayerDiffNode2[:,2]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='orange', label=r'$\rho_{u,v}^{2,2}(1)$')
axs.fill_between(x, y - yerr, y + yerr, color='orange', alpha=0.2)
x,y,yerr = DiffLayerDiffNode[:,0],DiffLayerDiffNode[:,2],EDiffLayerDiffNode[:,2]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='blue', label=r'$\rho_{u,v}^{1,2}(1)$')
axs.fill_between(x, y - yerr, y + yerr, color='blue', alpha=0.2)

axs.axvline(1/31,linestyle='--',color='red',label = r'$\frac{1}{\Lambda_1}$')
axs.axvline(1/11,linestyle='--',color='green',label = r'$\frac{1}{\Lambda_2}$')
axs.set_xlabel(r'$\eta/\beta$',fontsize=14)
axs.set_ylabel(r'$\rho(1)$',fontsize=14)
axs.set_xscale('log')
axs.legend(fontsize=12)

plt.tight_layout()
plt.savefig(str(IMAGES_DIR / 'EReta001.pdf'))
plt.show()

#%%
SamelayerSameNode1 = np.loadtxt(str(DATA_DIR / 'EReta25_mu05/final_corr1.txt'))
SamelayerSameNode2 = np.loadtxt(str(DATA_DIR / 'EReta25_mu05/final_corr2.txt'))
DiffLayerSameNode = np.loadtxt(str(DATA_DIR / 'EReta25_mu05/final_corr3.txt'))
SameLayerDiffNode1 = np.loadtxt(str(DATA_DIR / 'EReta25_mu05/final_corr4.txt'))
DiffLayerDiffNode = np.loadtxt(str(DATA_DIR / 'EReta25_mu05/final_corr5.txt'))
SamelayerDiffNode2 = np.loadtxt(str(DATA_DIR / 'EReta25_mu05/final_corr6.txt'))
SamelayerSameNode1,ESamelayerSameNode1 = groupby_mean_and_error(SamelayerSameNode1)
SamelayerSameNode2,ESamelayerSameNode2 = groupby_mean_and_error(SamelayerSameNode2)
DiffLayerSameNode,EDiffLayerSameNode = groupby_mean_and_error(DiffLayerSameNode)
SameLayerDiffNode1,ESameLayerDiffNode1 = groupby_mean_and_error(SameLayerDiffNode1)
DiffLayerDiffNode,EDiffLayerDiffNode = groupby_mean_and_error(DiffLayerDiffNode)
SamelayerDiffNode2,ESamelayerDiffNode2= groupby_mean_and_error(SamelayerDiffNode2)


SS1MSTEQ = np.loadtxt(str(DATA_DIR / 'eta25_mu05_N1000/final_corr1msteq.txt'))
SS2MSTEQ = np.loadtxt(str(DATA_DIR / 'eta25_mu05_N1000/final_corr2msteq.txt'))


fig, axs = plt.subplots(1, 1, figsize=(6, 5))

x,y,yerr = SamelayerSameNode1[:, 0],SamelayerSameNode1[:, 2],ESamelayerSameNode1[:, 2]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='red', label=r'$\rho_{u,u}^{1,1}(1)$')
axs.fill_between(x, y - yerr, y + yerr, color='red', alpha=0.2)
x,y,yerr = SamelayerSameNode2[:, 0],SamelayerSameNode2[:, 2],ESamelayerSameNode2[:, 2]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='green', label=r'$\rho_{u,u}^{2,2}(1)$')
axs.fill_between(x, y - yerr, y + yerr, color='green', alpha=0.2)
x,y,yerr = DiffLayerSameNode[:,0],DiffLayerSameNode[:,2],EDiffLayerSameNode[:,2]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='purple', label=r'$\rho_{u,u}^{1,2}(1)$')
axs.fill_between(x, y - yerr, y + yerr, color='purple', alpha=0.2)
x,y,yerr = SameLayerDiffNode1[:,0],SameLayerDiffNode1[:,2],ESameLayerDiffNode1[:,2]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='cyan', label=r'$\rho_{u,v}^{1,1}(1)$')
axs.fill_between(x, y - yerr, y + yerr, color='cyan', alpha=0.2)
x,y,yerr = SamelayerDiffNode2[:,0],SamelayerDiffNode2[:,2],ESamelayerDiffNode2[:,2]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='orange', label=r'$\rho_{u,v}^{2,2}(1)$')
axs.fill_between(x, y - yerr, y + yerr, color='orange', alpha=0.2)
x,y,yerr = DiffLayerDiffNode[:,0],DiffLayerDiffNode[:,2],EDiffLayerDiffNode[:,2]
axs.plot(x, y, linestyle='-', marker='o', markersize=2, color='blue', label=r'$\rho_{u,v}^{1,2}(1)$')
axs.fill_between(x, y - yerr, y + yerr, color='blue', alpha=0.2)

axs.axvline(1/47.93,linestyle='--',color='black',label = r'$\frac{1}{\Lambda_T}$')
axs.set_xlabel(r'$\eta/\beta$',fontsize=14)
axs.set_ylabel(r'$\rho(1)$',fontsize=14)
axs.set_xscale('log')
axs.legend(fontsize=12)

plt.tight_layout()
plt.savefig(str(IMAGES_DIR / 'EReta25.pdf'))
plt.show()
