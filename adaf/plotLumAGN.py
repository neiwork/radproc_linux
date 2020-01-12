import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from matplotlib import rc        # TO MANAGE MATPLOTLIB PARAMETERS"
rc('font',family='serif')
rc('text',usetex = True)

import scipy.optimize as optimization

nu,eV,Sy,Br,IC,pp,CD,Refl,Tot,ICin = np.loadtxt('lumThermal.dat',unpack=True)
NT_logeV,NT_logSye,NT_logSyp,NT_logIC,NT_logpp,NT_logpg,NT_logAbs = np.loadtxt('lumNonThermal.dat',unpack=True,skiprows=1)
x_eV = [-5,17]
y_axis = [37,45]

fig, ax1 = plt.subplots()

ax1.tick_params(axis='both',labelsize=12)
ax1.set_xlim(x_eV)
ax1.set_ylim(y_axis)

ax1.set_xlabel(r'$\mathrm{Log}(E/\mathrm{eV})$',fontsize=13)
ax1.set_ylabel(r'$\mathrm{Log}(\nu L_\nu / \mathrm{erg~s}^{-1})$',fontsize=13)

ax1.plot(np.log10(eV),np.log10(Tot),label='Thermal',lw=2,ls='solid')
ax1.plot(NT_logeV,NT_logSye,label='eSy',lw=1.3,ls='--',marker='s',markevery=4,markersize=5)
ax1.plot(NT_logeV,NT_logSyp,label='pSy',lw=1.3,ls='--',marker='o',markevery=4,markersize=5)
ax1.plot(NT_logeV,NT_logIC,label='IC',lw=1.3,ls='-.',marker='*',markevery=4,markersize=5)
ax1.plot(NT_logeV,NT_logpp,label='pp',lw=1.3,ls=':',marker='v',markevery=4,markersize=5)
ax1.plot(NT_logeV,NT_logpg,label=r'p$\gamma$',lw=1.3,ls=':',marker='^',markevery=4,markersize=5)
ax1.plot(NT_logeV,NT_logAbs,lw=3,label='Abs',color='k',ls='solid')
ax1.legend(loc='best',fontsize=8)

fig.savefig('nonThermalLum.pdf')
