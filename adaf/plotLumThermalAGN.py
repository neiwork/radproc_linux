import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from matplotlib import rc        # TO MANAGE MATPLOTLIB PARAMETERS

rc('font',family='serif')
rc('text',usetex = True)

import scipy.optimize as optimization

nu,eV,Sy,Br,IC,pp,CD,Refl,Tot,ICin = np.loadtxt('lumThermal.dat',unpack=True)

x_nu = [8,22]
y_axis = [35,42]

fig, ax1 = plt.subplots()

ax1.tick_params(axis='both',labelsize=12)
ax1.set_xlim(x_nu)
ax1.set_ylim(y_axis)

ax1.set_xlabel(r'$\mathrm{Log}(\nu/\mathrm{Hz})$',fontsize=13)
ax1.set_ylabel(r'$\mathrm{Log}(\nu L_\nu / \mathrm{erg~s}^{-1})$',fontsize=13)

ax1.plot(np.log10(nu),np.log10(Tot))

plt.savefig('thermalLum.eps')
