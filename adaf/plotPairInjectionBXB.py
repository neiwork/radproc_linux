import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from matplotlib import rc        # TO MANAGE MATPLOTLIB PARAMETERS"
rc('font',family='serif')
rc('text',usetex = True)

from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
viridis = cm.get_cmap('viridis', 12)

logeV1,ra,logQ1 = np.loadtxt('secondaryPairsInj.txt',unpack=True)
logeV2,rb,logQ2 = np.loadtxt('secondaryPairsInj2.txt',unpack=True)
logerg1 = np.log10(np.power(10,logeV1)*1.6e-12)
logerg2 = np.log10(np.power(10,logeV2)*1.6e-12)

x_eV = [logeV1[0],logeV1[-1]]
y = [0,24]

nR = 30
f = 4
nE = 80
colors = np.arange(nR)/nR

fig, ax1 = plt.subplots()

ax1.tick_params(axis='both',labelsize=14)
ax1.set_xlim(x_eV)
ax1.set_ylim(y)
ax1.set_xlabel(r'$\mathrm{Log}(E ~ [\mathrm{eV}])$',fontsize=14)
ax1.set_ylabel(r'$\mathrm{Log}(E^2Q_\mathrm{e^+e^-}(E) ~ [\mathrm{erg~cm^{-3}s^{-1}}])$',fontsize=14)

for r1 in np.arange(nR//f):
    ax1.plot(logeV1[f*r1*nE:(f*r1+1)*nE],2.0*logerg1[f*r1*nE:(f*r1+1)*nE] +
             logQ1[f*r1*nE:(f*r1+1)*nE],color=viridis(colors[r1*f]))
    ax1.plot(logeV2[f*r1*nE:(f*r1+1)*nE],2.0*logerg2[f*r1*nE:(f*r1+1)*nE] +
             logQ2[f*r1*nE:(f*r1+1)*nE],color=viridis(colors[r1*f]))

fig.savefig('pairInj_energy.pdf')


fig, ax1 = plt.subplots()

ax1.tick_params(axis='both',labelsize=14)
ax1.set_xlim(x_eV)
ax1.set_ylim(y)
ax1.set_xlabel(r'$\mathrm{Log}(E ~ [\mathrm{eV}])$',fontsize=14)
ax1.set_ylabel(r'$\mathrm{Log}(Q_\mathrm{e^+e^-}(E) ~ [\mathrm{erg^{-1}cm^{-3}s^{-1}}])$',fontsize=14)

for r1 in np.arange(nR//f):
    ax1.plot(logeV1[f*r1*nE:(f*r1+1)*nE], logQ1[f*r1*nE:(f*r1+1)*nE],color=viridis(colors[r1*f]))
    ax1.plot(logeV2[f*r1*nE:(f*r1+1)*nE], logQ2[f*r1*nE:(f*r1+1)*nE],color=viridis(colors[r1*f]))

fig.savefig('pairInj.pdf')
