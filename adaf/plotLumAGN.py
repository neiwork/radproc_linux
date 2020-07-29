#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
plt.style.use('mystyle')

import scipy.interpolate as interpol


nu,eV,Sy,Br,IC,pp,ColdDisk,Refl,Tot,ICin = np.loadtxt('lumThermal.dat',unpack=True)
NT_logeV,NT_logSyp,NT_logIC,NT_logpp,NT_logpg,NT_logNotAbs, NT_logAbs = np.loadtxt('lumNonThermal.dat',unpack=True,skiprows=1)
NT_logeVSy,NT_logSye = np.loadtxt('lumSy.dat',unpack=True,skiprows=1)
NT_logeVs,NT_logSys,NT_logSymu,NT_logSypi,NT_logICp,NT_logpip,NT_logpig,a,NT_logAbss = np.loadtxt('secondariesLum.dat',unpack=True,skiprows=1)

NT_logNotAbs_fun = interpol.interp1d(NT_logeV,NT_logNotAbs,bounds_error=False,fill_value=-1000.0)
NT_logAbss_fun = interpol.interp1d(NT_logeVs,NT_logAbss,bounds_error=False,fill_value=-1000.0)
NT_logSye_fun = interpol.interp1d(NT_logeVSy,NT_logSye,bounds_error=False,fill_value=-1000.0)
Th_logTot_fun = interpol.interp1d(np.log10(eV),np.log10(Tot),bounds_error=False,fill_value=-1000.0)

x_eV = [-6,17]
y_axis = [36,42]

fig, ax1 = plt.subplots(figsize=(15,7))

ax1.tick_params(axis='both',labelsize=30)
ax1.set_xlim(x_eV)
ax1.set_ylim(y_axis)

ax1.set_xlabel(r'$\mathrm{Log}(E/\mathrm{eV})$',fontsize=30)
ax1.set_ylabel(r'$\mathrm{Log}(\nu L_\nu / \mathrm{erg~s}^{-1})$',fontsize=30)

markersize=10


ax1.plot(np.log10(eV),np.log10(Tot),label='Thermal',lw=4,ls='solid')
#ax1.plot(np.log10(eV),np.log10(ColdDisk),label='Thin disk',lw=4,ls='solid')
ax1.plot(NT_logeVSy,NT_logSye,label='eSy',lw=2,ls='-',marker='v',markevery=4,markersize=markersize)
ax1.plot(NT_logeV,NT_logSyp,label='pSy',lw=2,ls='-.',marker='v',markevery=4,markersize=markersize)
ax1.plot(NT_logeV,NT_logSys,label='pairSy',lw=2,ls='--',marker='v',markevery=4,markersize=markersize)
ax1.plot(NT_logeV,NT_logSymu,label='muSy',lw=2,ls=':',marker='v',markevery=4,markersize=markersize)
ax1.plot(NT_logeV,NT_logICp,label='pairIC',lw=2,ls='--',marker='o',markevery=4,markersize=5)
ax1.plot(NT_logeV,NT_logSypi,label='piSy',lw=2,ls='--',marker='v',markevery=4,markersize=5)
ax1.plot(NT_logeV,NT_logpip,label='pi-p',lw=2,ls='--',marker='v',markevery=4,markersize=5)
ax1.plot(NT_logeV,NT_logpig,label='pi-g',lw=2,ls='--',marker='v',markevery=4,markersize=5)
ax1.plot(NT_logeV,NT_logIC,label='IC',lw=2,ls='--',marker='*',markevery=4,markersize=markersize)
ax1.plot(NT_logeV,NT_logpp,label='pp',lw=2,ls=':',marker='s',markevery=4,markersize=markersize)
ax1.plot(NT_logeV,NT_logpg,label=r'p$\gamma$',lw=2,ls='-.',marker='o',markevery=4,markersize=markersize)

logeV_tot = np.log10(np.logspace(-6,16,500))
logTotNotAbs = np.log10(np.power(10,NT_logNotAbs_fun(logeV_tot)) + np.power(10,NT_logAbss_fun(logeV_tot)) +
                        np.power(10,Th_logTot_fun(logeV_tot)) + np.power(10,NT_logSye_fun(logeV_tot)))
#logTotNotAbs = np.log10(np.power(10,NT_logNotAbs_fun(logeV_tot)) + np.power(10,Th_logTot_fun(logeV_tot)) + np.power(10,NT_logSye_fun(logeV_tot)))

ax1.plot(logeV_tot,logTotNotAbs,lw=5,label='Total Abs',color='k',ls='solid')
ax1.legend(loc='best',fontsize=20)

plt.tight_layout()
fig.savefig('nonThermalLum.pdf')
