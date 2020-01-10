import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from matplotlib import rc        # TO MANAGE MATPLOTLIB PARAMETERS"
rc('font',family='serif')
rc('text',usetex = True)

from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
viridis = cm.get_cmap('viridis', 12)

logeVp,r,lognp = np.loadtxt('protonDistribution.dat',unpack=True)
logeVe,r,logne = np.loadtxt('electronDistribution.dat',unpack=True)
logerge = np.log10(np.power(10,logeVe)*1.6e-12)
logergp = np.log10(np.power(10,logeVp)*1.6e-12)

log_gp = logeVp - 8.973
log_ge = logeVe - 5.71

x_ge = [log_ge[0]+.5,log_ge[-1]]
x_gp = [log_gp[0]+.5,log_gp[-1]]
y_p = [40,50]
y_e = [30,45]

nR = 30
f = 1
nE = 100
colors = np.arange(nR)/nR

fig, ax1 = plt.subplots()

ax1.tick_params(axis='both',labelsize=15)
ax1.set_xlim(x_ge)
ax1.set_ylim(y_e)

ax1.set_xlabel(r'$\mathrm{Log}~\gamma_e$',fontsize=19)
ax1.set_ylabel(r'$\mathrm{Log}(E^2 N_\mathrm{e}(E) ~ [\mathrm{erg~cm^{-3}}])$',fontsize=19)

neTot = np.zeros(nE)
for r1 in np.arange(nR//f):
    ax1.plot(log_ge[f*r1*nE:(f*r1+1)*nE],2.0*logerge[f*r1*nE:(f*r1+1)*nE] +
             logne[f*r1*nE:(f*r1+1)*nE],color=viridis(colors[r1*f]))

for r1 in np.arange(nR):
     for e in np.arange(nE):
        neTot[e] = neTot[e] + np.power(10,logne[r1*nE+e])

    
ax1.plot(log_ge[:nE],2.0*logerge[:nE]+np.log10(neTot),ls='-',lw=4,color='red',label='Total')
ax1.legend(loc='best',fontsize=15)
plt.tight_layout()

fig.savefig('electronDist.pdf')

fig, ax1 = plt.subplots()

ax1.tick_params(axis='both',labelsize=15)
ax1.set_xlim(x_gp)
ax1.set_ylim(y_p)
ax1.set_xlabel(r'$\mathrm{Log}~\gamma_p$',fontsize=19)
ax1.set_ylabel(r'$\mathrm{Log}(E^2 N_\mathrm{p}(E) ~ [\mathrm{erg~cm^{-3}}])$',fontsize=19)

npTot = np.zeros(nE)

for r1 in np.arange(nR//f):
    ax1.plot(log_gp[f*r1*nE:(f*r1+1)*nE],2.0*logergp[f*r1*nE:(f*r1+1)*nE] +
             lognp[f*r1*nE:(f*r1+1)*nE],color=viridis(colors[r1*f]))

for r1 in np.arange(nR):
     for e in np.arange(nE):
        npTot[e] = npTot[e] + np.power(10,lognp[r1*nE+e])

ax1.plot(log_gp[:nE],2.0*logergp[:nE]+np.log10(npTot),ls='-',lw=4,color='red',label='Total')
ax1.legend(loc='best',fontsize=15)
plt.tight_layout()

fig.savefig('protonDist.pdf')
