import matplotlib.pyplot as plt
import numpy as np

# Read the data
clttPath = '/mnt/c/Users/gdzhao/projects/unwise_sz/Planck_temperature/namaster_cltt.txt'

# Read the data
ell, cltt = np.loadtxt(clttPath, unpack=True,skiprows=1)
pll = ell*ell*(ell+1)/(2*np.pi)

fig,ax = plt.subplots(figsize=(8,6))

plotmask = ell < 2e3
ax.plot(ell[plotmask], cltt[plotmask]*pll[plotmask], label='namaster', color='blue', lw=2)

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlabel(r'$\ell$', fontsize=16)

ax.set_ylabel(r'$\ell(\ell+1)C_{\ell}^{TT}/(2\pi)$', fontsize=16)

ax.legend(fontsize=14)

plt.show()

# Save the plot
plotPath = '/mnt/c/Users/gdzhao/projects/unwise_sz/Planck_temperature/namaster_cltt.png'
fig.savefig(plotPath, dpi=300, bbox_inches='tight')