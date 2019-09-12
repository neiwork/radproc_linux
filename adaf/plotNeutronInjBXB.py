import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from matplotlib import rc        # TO MANAGE MATPLOTLIB PARAMETERS"
rc('font',family='serif')
rc('text',usetex = True)

from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
viridis = cm.get_cmap('viridis', 12)
magma = cm.get_cmap('magma',12)

logeV,r,logpp,logpY = np.loadtxt('neutronInj.txt',unpack=True,skiprows=1)
logerg = np.log10(np.power(10,logeV)*1.6e-12)

x_eV = [logeV[0],logeV[-1]]
y = [-5,16]

nR = 30
f = 2
nE = 100
colors = np.arange(nR)/nR

fig, ax1 = plt.subplots()

ax1.tick_params(axis='both',labelsize=14)
ax1.set_xlim(x_eV)
ax1.set_ylim(y)
ax1.set_xlabel(r'$\mathrm{Log}(E ~ [\mathrm{eV}])$',fontsize=14)
ax1.set_ylabel(r'$\mathrm{Log}(E^2Q_\mathrm{n}(E) ~ [\mathrm{erg~cm^{-3}s^{-1}}])$',fontsize=14)

qTot = np.zeros(nE)
for r1 in np.arange(nR//f):
    ax1.plot(logeV[f*r1*nE:(f*r1+1)*nE],2.0*logerg[f*r1*nE:(f*r1+1)*nE] +
             logpp[f*r1*nE:(f*r1+1)*nE],color=viridis(colors[r1*f]))
    ax1.plot(logeV[f*r1*nE:(f*r1+1)*nE],2.0*logerg[f*r1*nE:(f*r1+1)*nE] +
             logpY[f*r1*nE:(f*r1+1)*nE],color=magma(colors[r1*f]))

for e in np.arange(nE):
    for r1 in np.arange(nR):
        qTot[e] = qTot[e] + np.power(10,logpp[r1*nE+e]) + np.power(10,logpY[r1*nE+e])

ax1.plot(logeV[:nE],2.0*logerg[:nE]+np.log10(qTot),ls='--',lw=3,label='Total', \
            color = 'k')
ax1.legend(loc='best',fontsize=15)
fig.savefig('neutronInj_energy.pdf')





fig, ax1 = plt.subplots()

ax1.tick_params(axis='both',labelsize=14)
ax1.set_xlim(x_eV)
ax1.set_ylim(y)
ax1.set_xlabel(r'$\mathrm{Log}(E ~ [\mathrm{eV}])$',fontsize=14)
ax1.set_ylabel(r'$\mathrm{Log}(Q_\mathrm{n}(E) ~ [\mathrm{erg^{-1}cm^{-3}s^{-1}}])$',fontsize=14)

qTot = np.zeros(nE)
for r1 in np.arange(nR//f):
    ax1.plot(logeV[f*r1*nE:(f*r1+1)*nE], logpp[f*r1*nE:(f*r1+1)*nE],color=viridis(colors[r1*f]))
    ax1.plot(logeV[f*r1*nE:(f*r1+1)*nE], logpY[f*r1*nE:(f*r1+1)*nE],color=magma(colors[r1*f]))

for e in np.arange(nE):
    for r1 in np.arange(nR):
        qTot[e] = qTot[e] + np.power(10,logpp[r1*nE+e]) + np.power(10,logpY[r1*nE+e])

ax1.plot(logeV[:nE],np.log10(qTot),ls='--',lw=3,label='Total', \
            color = 'k')
ax1.legend(loc='best',fontsize=15)
fig.savefig('neutronInj.pdf')
