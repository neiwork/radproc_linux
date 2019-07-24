import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from matplotlib import rc        # TO MANAGE MATPLOTLIB PARAMETERS"
rc('font',family='serif')
rc('text',usetex = True)

nu,eV,Sy,Br,IC,pp,CD,Refl,Tot = np.loadtxt('Release/lum.txt',unpack=True)
AGN_nu = [1.0e10,1.0e22]
AGN_eV = np.log10(np.array(AGN_nu)/2.418e14)
AGN_nu = np.log10(AGN_nu)
BXB_nu = [1.0e14,1.0e21]
BXB_eV = np.log10(np.array(BXB_nu)/2.418e14)
BXB_nu = np.log10(BXB_nu)
AGN = [37,44]
BXB = [32,38]

fig, ax1 = plt.subplots()
ax2 = ax1.twiny()

ax1.tick_params(axis='both', labelsize=12)
ax1.set_xlim(AGN_nu)
ax1.set_ylim(AGN)
ax1.set_xlabel(r'$\mathrm{Log}(\nu/\mathrm{Hz})$',fontsize=13)
ax1.xaxis.set_label_position('bottom')

ax2.set_xlabel(r'$\mathrm{Log}(E/\mathrm{eV})$',fontsize=13)
ax2.xaxis.set_label_position('top')
ax2.set_xlim(AGN_eV)

ax1.set_ylabel(r'$\mathrm{Log}(\nu L_\nu / \mathrm{erg~s}^{-1})$',fontsize=13)

ax1.plot(np.log10(nu),np.log10(Sy),label='Sy')
ax1.plot(np.log10(nu),np.log10(Br),label='Br')
ax1.plot(np.log10(nu),np.log10(IC),label='IC')
ax1.plot(np.log10(nu),np.log10(CD),label='Cold Disk')
ax1.plot(np.log10(nu),np.log10(Refl),label='Reflection')
ax1.plot(np.log10(nu),np.log10(Tot+CD+Refl),lw=3,label='Tot')

fig.savefig('lum.pdf')

r,Sy,Br,IC,pp,Tot = np.loadtxt('lumRadius.txt',unpack=True)

fig, ax = plt.subplots()

ax.tick_params(axis='both', labelsize=12)
ax.set_xlim([0.2,1.9])
ax.set_ylim(BXB)
ax.set_xlabel(r'$\mathrm{Log}(r/\mathrm{M})$',fontsize=13)
ax.xaxis.set_label_position('bottom')

ax.set_ylabel(r'$\mathrm{Log}(L / \mathrm{erg~s}^{-1})$',fontsize=13)

rS = 3*1.0e5*2.0e8
ax.plot(np.log10(r/rS),np.log10(Sy),c='r')
ax.plot(np.log10(r/rS),np.log10(Br),label='Br')
ax.plot(np.log10(r/rS),np.log10(IC),label='IC')
ax.plot(np.log10(r/rS),np.log10(Sy+Br+IC),label='Tot')

fig.savefig('lumRadius.pdf')

J,Tot = np.loadtxt('Release/photonDensity1.txt',unpack=True)

fig, ax1 = plt.subplots()
ax2 = ax1.twiny()

BXB_J = [-20,-12]
AGN_J = [-18,-12]
BXB_n = [5,15]
AGN_n = [0,15]
ax1.tick_params(axis='both', labelsize=12)
ax1.set_xlim(AGN_J)
ax1.set_ylim(AGN_n)
ax1.set_xlabel(r'$\mathrm{Log}(E/\mathrm{J})$',fontsize=13)
ax1.xaxis.set_label_position('bottom')

ax2.set_xlabel(r'$\mathrm{Log}(E/\mathrm{eV})$',fontsize=13)
ax2.xaxis.set_label_position('top')
ax2.set_xlim(AGN_eV)

ax1.set_ylabel(r'$\mathrm{Log}(E n_E / \mathrm{m}^{-3})$',fontsize=13)

ax1.plot(np.log10(J),np.log10(J*Tot))

fig.savefig('photonDensity.pdf')
