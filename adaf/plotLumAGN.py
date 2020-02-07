import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from matplotlib import rc        # TO MANAGE MATPLOTLIB PARAMETERS"
rc('font',family='serif')
rc('text',usetex = True)

import scipy.optimize as optimization

nu,eV,Sy,Br,IC,pp,CD,Refl,Tot,ICin = np.loadtxt('lumThermal.dat',unpack=True)
NT_logeV,NT_logSyp,NT_logIC,NT_logpp,NT_logpg,NT_logAbs = np.loadtxt('lumNonThermal.dat',unpack=True,skiprows=1)
NT_logeVSy,NT_logSye = np.loadtxt('lumSy.dat',unpack=True,skiprows=1)
NT_logeVs,NT_logSys,NT_logSymu,NT_logSypi,NT_logICp,NT_logpip,NT_logpig,NT_logAbss = np.loadtxt('secondariesLum.dat',unpack=True,skiprows=1)

x_eV = [-5,17]
y_axis = [37,45]

fig, ax1 = plt.subplots(figsize=(15,7))

ax1.tick_params(axis='both',labelsize=30)
ax1.set_xlim(x_eV)
ax1.set_ylim(y_axis)

ax1.set_xlabel(r'$\mathrm{Log}(E/\mathrm{eV})$',fontsize=30)
ax1.set_ylabel(r'$\mathrm{Log}(\nu L_\nu / \mathrm{erg~s}^{-1})$',fontsize=30)

ax1.plot(np.log10(eV),np.log10(Tot),label='Thermal',lw=4,ls='solid')
ax1.plot(NT_logeVSy,NT_logSye,label='eSy',lw=3,ls='--',marker='s',markevery=4,markersize=5)
ax1.plot(NT_logeV,NT_logSyp,label='pSy',lw=3,ls='--',marker='o',markevery=4,markersize=5)
ax1.plot(NT_logeV,NT_logSys,label='pairSy',lw=3,ls='--',marker='^',markevery=4,markersize=5)
ax1.plot(NT_logeV,NT_logSymu,label='muSy',lw=3,ls='--',marker='^',markevery=4,markersize=5)
ax1.plot(NT_logeV,NT_logICp,label='pairIC',lw=3,ls='--',marker='o',markevery=4,markersize=5)
ax1.plot(NT_logeV,NT_logIC,label='IC',lw=3,ls='-.',marker='*',markevery=4,markersize=5)
ax1.plot(NT_logeV,NT_logpp,label='pp',lw=3,ls=':',marker='v',markevery=4,markersize=5)
ax1.plot(NT_logeV,NT_logpg,label=r'p$\gamma$',lw=3,ls=':',marker='^',markevery=4,markersize=5)
logTotSoft = np.log10( np.power(10,NT_logSye) + Tot)
logTotHard = np.log10(np.power(10,NT_logAbs)+np.power(10,NT_logAbss))
ax1.plot(NT_logeV[45:],logTotHard[45:],lw=8,label='Total Abs',color='k',ls='solid')
ax1.plot(NT_logeVSy[:90],logTotSoft[:90],lw=8,label='Total Abs',color='k',ls='solid')
ax1.legend(loc='best',fontsize=20)

plt.tight_layout()
fig.savefig('nonThermalLum.pdf')
