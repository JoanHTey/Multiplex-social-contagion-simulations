import numpy as np
import matplotlib.pyplot as plt


    

# Time steps where beta changes
timestep = np.linspace(0, 25000, 6)

# Load your correlation data
c1 = np.load('corr1.npy')
c2 = np.load('corr2.npy')
c3 = np.load('corr3.npy')
c4 = np.load('corr4.npy')
c5 = np.load('corr5.npy')  # You had c4 twice, assuming this was a typo
c6 = np.load('corr6.npy')
ec1 = np.load('ecorr1.npy')
ec2 = np.load('ecorr2.npy')
ec3 = np.load('ecorr3.npy')
ec4 = np.load('ecorr4.npy')
ec5 = np.load('ecorr5.npy')
ec6 = np.load('ecorr6.npy')
d1 = np.load('den1.npy')
d2 = np.load('den2.npy')


# Time axis
t = np.linspace(0, 24500, 50)
print(t)
# Initialize figure with two subplots
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# First subplot: correlations
ax1 = axes[0]
ax1.plot(t, c1[:, 1], label=r'$\rho_{uu}^{11}(1)$', color='tab:blue',marker='o',markersize=5, linewidth=1.5)
ax1.fill_between(t, c1[:, 1] - ec1[:, 1], c1[:, 1] + ec1[:, 1], color='blue', alpha=0.2)
ax1.plot(t, c2[:, 1], label=r'$\rho_{uu}^{22}(1)$', color='tab:orange',marker='o',markersize=5, linewidth=1.5)
ax1.fill_between(t, c2[:, 1] - ec2[:, 1], c2[:, 1] + ec2[:, 1], color='orange', alpha=0.2)
ax1.plot(t, c3[:, 1], label=r'$\rho_{uu}^{12}(1)$', color='tab:green',marker='o',markersize=5, linewidth=1.5)
ax1.fill_between(t, c3[:, 1] - ec3[:, 1], c3[:, 1] + ec3[:, 1], color='green', alpha=0.2)
ax1.plot(t, c4[:, 1], label=r'$\rho_{uv}^{11}(1)$', color='tab:red',marker='o',markersize=5, linewidth=1.5)
ax1.fill_between(t, c4[:, 1] - ec4[:, 1], c4[:, 1] + ec4[:, 1], color='red', alpha=0.2)
ax1.plot(t, c5[:, 1], label=r'$\rho_{uv}^{22}(1)$', color='tab:purple',marker='o',markersize=5, linewidth=1.5)
ax1.fill_between(t, c5[:, 1] - ec5[:, 1], c5[:, 1] + ec5[:, 1], color='purple', alpha=0.2)
ax1.plot(t, c6[:, 1], label=r'$\rho_{uv}^{12}(1)$', color='tab:brown',marker='o',markersize=5, linewidth=1.5)

# Add vertical lines for beta change points
ax1.axvline(0, linestyle='--', color='red',linewidth=1, alpha=0.7,label=fr'$\beta={np.linspace(0.01,0.15,5)[0]:.2f}$')
ax1.axvline(5000, linestyle='--',color='green', linewidth=1, alpha=0.7,label=fr'$\beta={np.linspace(0.01,0.15,5)[1]:.2f}$')
ax1.axvline(10000, linestyle='--',color='purple', linewidth=1, alpha=0.7,label=fr'$\beta={np.linspace(0.01,0.15,5)[2]:.2f}$')
ax1.axvline(15000, linestyle='--',color='blue', linewidth=1, alpha=0.7,label=fr'$\beta={np.linspace(0.01,0.15,5)[3]:.2f}$')
ax1.axvline(20000, linestyle='--',color='yellow', linewidth=1, alpha=0.7,label=fr'$\beta={np.linspace(0.01,0.15,5)[4]:.2f}$')

# Axis labels and title
ax1.set_xlabel('Time steps', fontsize=14)
ax1.set_ylabel(r'$\rho_{uv}^{\alpha\kappa}(1)$', fontsize=14)

# Legend and layout
ax1.legend(fontsize=10, loc='best')
ax1.set_ylim(0, 1.05)

# Second subplot: densities
ax2 = axes[1]
ax2.plot(t, d1, label=r'$\langle X \rangle_1$', color='tab:blue',marker='o',markersize=5, linewidth=1.5)
ax2.plot(t, d2, label=r'$\langle X \rangle_2$', color='tab:orange',marker='o',markersize=5, linewidth=1.5)
# Add vertical lines for beta change points
ax2.axvline(0, linestyle='--', color='red',linewidth=1, alpha=0.7,label=fr'$\beta={np.linspace(0.01,0.15,5)[0]:.2f}$')
ax2.axvline(5000, linestyle='--',color='green', linewidth=1, alpha=0.7,label=fr'$\beta={np.linspace(0.01,0.15,5)[1]:.2f}$')
ax2.axvline(10000, linestyle='--',color='purple', linewidth=1, alpha=0.7,label=fr'$\beta={np.linspace(0.01,0.15,5)[2]:.2f}$')
ax2.axvline(15000, linestyle='--',color='blue', linewidth=1, alpha=0.7,label=fr'$\beta={np.linspace(0.01,0.15,5)[3]:.2f}$')
ax2.axvline(20000, linestyle='--',color='yellow', linewidth=1, alpha=0.7,label=fr'$\beta={np.linspace(0.01,0.15,5)[4]:.2f}$')

ax2.legend(fontsize=10, loc='best')

ax2.set_xlabel('Time steps', fontsize=14)
ax2.set_ylabel(r'$\langle X \rangle_i$', fontsize=14)

plt.tight_layout()


plt.savefig('expe.pdf')
# Show figure
plt.show()
