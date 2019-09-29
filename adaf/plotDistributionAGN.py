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
y_p = [-15,3]
y_e = [-3,0]

nR = 200
f = 1
nE = 60
colors = np.arange(nR)/nR


fig, ax1 = plt.subplots()

ax1.tick_params(axis='both',labelsize=14)
ax1.set_xlim(x_eVe)
ax1.set_ylim(y_e)

ax1.set_xlabel(r'$\mathrm{Log}(E ~ [\mathrm{eV}])$',fontsize=14)
ax1.set_ylabel(r'$\mathrm{Log}(E^2 N_\mathrm{e}(E) ~ [\mathrm{erg~cm^{-3}}])$',fontsize=14)

neTot = np.zeros(nE)
for r1 in np.arange(nR//f):
    ax1.plot(logeVe[f*r1*nE:(f*r1+1)*nE],2.0*logerge[f*r1*nE:(f*r1+1)*nE] +
             logne[f*r1*nE:(f*r1+1)*nE],color=viridis(colors[r1*f]))

for r1 in np.arange(nR):
     for e in np.arange(nE):
        neTot[e] = neTot[e] + np.power(10,logne[r1*nE+e])

    
ax1.plot(logeVe[:nE],2.0*logerge[:nE]+np.log10(neTot),ls='--',lw=3,label='Total', \
            color = 'k')
ax1.legend(loc='lower left',fontsize=15)
fig.savefig('electronDist.pdf')

fig, ax1 = plt.subplots()

ax1.tick_params(axis='both',labelsize=14)
ax1.set_xlim(x_eVp)
ax1.set_ylim(y_p)
ax1.set_xlabel(r'$\mathrm{Log}(E ~ [\mathrm{eV}])$',fontsize=14)
ax1.set_ylabel(r'$\mathrm{Log}(E^2 N_\mathrm{p}(E) ~ [\mathrm{erg~cm^{-3}}])$',fontsize=14)

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
fig.savefig('protonDist.pdf')


logeVp,r,lognp_le = np.loadtxt("linearEmitter_p.txt",unpack=True)

fig, ax1 = plt.subplots()

ax1.tick_params(axis='both',labelsize=14)
ax1.set_xlim(x_eVp)
ax1.set_ylim(y_p)

ax1.set_xlabel(r'$\mathrm{Log}(E ~ [\mathrm{eV}])$',fontsize=14)
ax1.set_ylabel(r'$\mathrm{Log}(E^2 N_\mathrm{p}(E) ~ [\mathrm{erg~cm^{-3}}])$',fontsize=14)

npTot_le = np.zeros(nE)
for r1 in np.arange(nR//f):
    ax1.plot(logeVp[f*r1*nE:(f*r1+1)*nE],2.0*logergp[f*r1*nE:(f*r1+1)*nE] +
             lognp_le[f*r1*nE:(f*r1+1)*nE],color=viridis(colors[r1*f]))

for r1 in np.arange(nR):
     for e in np.arange(nE):
        npTot_le[e] = npTot_le[e] + np.power(10,lognp_le[r1*nE+e])

ax1.plot(logeVp[:nE],2.0*logergp[:nE]+np.log10(npTot_le),ls='--',lw=3,label='Total', \
            color = 'k')
ax1.legend(loc='lower left',fontsize=15)
fig.savefig('protonLinearEmitter.pdf')



logeVe,r,logne_le = np.loadtxt("linearEmitter_e.txt",unpack=True)

fig, ax1 = plt.subplots()

ax1.tick_params(axis='both',labelsize=14)
ax1.set_xlim(x_eVe)
ax1.set_ylim(y_e)

ax1.set_xlabel(r'$\mathrm{Log}(E ~ [\mathrm{eV}])$',fontsize=14)
ax1.set_ylabel(r'$\mathrm{Log}(E^2 N_\mathrm{e}(E) ~ [\mathrm{erg~cm^{-3}}])$',fontsize=14)

neTot_le = np.zeros(nE)
for r1 in np.arange(nR//f):
    ax1.plot(logeVe[f*r1*nE:(f*r1+1)*nE],2.0*logerge[f*r1*nE:(f*r1+1)*nE] +
             logne_le[f*r1*nE:(f*r1+1)*nE],color=viridis(colors[r1*f]))

for r1 in np.arange(nR):
     for e in np.arange(nE):
        neTot_le[e] = neTot_le[e] + np.power(10,logne_le[r1*nE+e])

ax1.plot(logeVe[:nE],2.0*logerge[:nE]+np.log10(neTot_le),ls='--',lw=3,label='Total', \
            color = 'k')
ax1.legend(loc='lower left',fontsize=15)
fig.savefig('electronLinearEmitter.pdf')
