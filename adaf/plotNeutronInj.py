import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from matplotlib import rc        # TO MANAGE MATPLOTLIB PARAMETERS"
rc('font',family='serif')
rc('text',usetex = True)

from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
viridis = cm.get_cmap('viridis', 12)

logeVp,r,lognp = np.loadtxt('protonDis.txt',unpack=True)
logeVe,r,logne = np.loadtxt('electronDis.txt',unpack=True)
logerge = np.log10(np.power(10,logeVe)*1.6e-12)
logergp = np.log10(np.power(10,logeVp)*1.6e-12)

x_eVe = [logeVe[0]+.5,logeVe[-1]]
x_eVp = [logeVp[0]+.5,logeVp[-1]]
y_p = [0,15]
y_e = [-10,10]

nR = 20
f = 1
nE = 50
colors = np.arange(nR)/nR



fig, ax1 = plt.subplots()

ax1.tick_params(axis='both',labelsize=14)
ax1.set_xlim(x_eVp)
ax1.set_ylim(y_p)
ax1.set_xlabel(r'$\mathrm{Log}(E ~ [\mathrm{eV}])$',fontsize=14)
ax1.set_ylabel(r'$\mathrm{Log}(N_\mathrm{p} ~ [\mathrm{erg~cm^{-3}}])$',fontsize=14)

npTot = np.zeros(nE)
for r1 in np.arange(nR//f):
    ax1.plot(logeVp[f*r1*nE:(f*r1+1)*nE],2.0*logergp[f*r1*nE:(f*r1+1)*nE] +
             lognp[f*r1*nE:(f*r1+1)*nE],color=viridis(colors[r1*f]))

for e in np.arange(nE):
    for r1 in np.arange(nR):
        npTot[e] = npTot[e] + np.power(10,lognp[r1*nE+e])

ax1.plot(logeVp[:nE],2.0*logergp[:nE]+np.log10(npTot),ls='--',lw=3,label='Total', \
            color = 'k')
ax1.legend(loc='lower left',fontsize=15)
fig.savefig('neutronInj.pdf')
