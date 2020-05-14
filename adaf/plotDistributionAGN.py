#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from matplotlib import rc        # TO MANAGE MATPLOTLIB PARAMETERS"
rc('font',family='serif')
rc('text',usetex = True)

from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
viridis = cm.get_cmap('viridis', 12)

logeVp,iR,logr,lognp = np.loadtxt('protonDistribution.dat',unpack=True)
logeVe,iR,logr,logne = np.loadtxt('electronDistribution.dat',unpack=True)
logeVp,iR,logr,logqp = np.loadtxt('protonInjection.dat',unpack=True)
logeVe,iR,logr,logqe = np.loadtxt('electronInjection.dat',unpack=True)
logeVpion,iR,logr,logqpion = np.loadtxt('pionInjection.dat',unpack=True)
logeVpion,iR,logr,lognpion = np.loadtxt('pionDistribution.dat',unpack=True)
logeVmuon,iR,logr,logqmuon = np.loadtxt('muonInjection.dat',unpack=True)
logeVmuon,iR,logr,lognmuon = np.loadtxt('muonDistribution.dat',unpack=True)
logeVpair,iR,logr,logqpair = np.loadtxt('secondaryPairInjection.dat',unpack=True)
NT_logeV,NT_logSyp,NT_logIC,NT_logpp,NT_logpg,NT_logAbs,NT_logAbsorbido = np.loadtxt('lumNonThermal.dat',unpack=True,skiprows=1)
NT_logeVs,NT_logSys,NT_logSymu,NT_logSypi,NT_logICs,NT_logpip,NT_logpig,NT_logAbs = np.loadtxt('secondariesLum.dat',unpack=True,skiprows=1)

logerge = np.log10(np.power(10,logeVe)*1.6e-12)
logergp = np.log10(np.power(10,logeVp)*1.6e-12)
logergpion = np.log10(np.power(10,logeVpion)*1.6e-12)
logergmuon = np.log10(np.power(10,logeVmuon)*1.6e-12)
logergpair = np.log10(np.power(10,logeVpair)*1.6e-12)

log_gp = logeVp - 8.973
log_ge = logeVe - 5.71
log_gpion = logeVpion - 8.145

x_ge = [log_ge[0]+.5,log_ge[-1]]
x_gp = [0,log_gp[-1]]
y_p = [33,50]
y_pq = [30,45]
y_e = [30,45]
y_eq = [35,44]

nR = 30
f = 1
nE = 100
colors = np.arange(nR)/nR

fig, ax1 = plt.subplots()

ax1.tick_params(axis='both',labelsize=15)
ax1.set_xlim(x_gp)
ax1.set_ylim(y_p)
ax1.set_xlabel(r'$\mathrm{Log}~\gamma_p$',fontsize=19)
ax1.set_ylabel(r'$\mathrm{Log}(E^2 N_\mathrm{p}(E) ~ [\mathrm{erg~cm^{-3}}])$',fontsize=19)

npTot = np.zeros(nE)
npionTot = np.zeros(nE)
nmuonTot = np.zeros(nE)

for r1 in np.arange(nR//f):
    ax1.plot(log_gp[f*r1*nE:(f*r1+1)*nE],2.0*logergp[f*r1*nE:(f*r1+1)*nE] +
             lognp[f*r1*nE:(f*r1+1)*nE],color=viridis(colors[r1*f]))

for r1 in np.arange(nR):
     for e in np.arange(nE):
        npionTot[e] = npionTot[e] + np.power(10,lognpion[r1*nE+e])
        npTot[e] = npTot[e] + np.power(10,lognp[r1*nE+e])
        nmuonTot[e] = nmuonTot[e] + np.power(10,lognmuon[r1*nE+e])

ax1.plot(log_gp[:nE],2.0*logergp[:nE]+np.log10(npTot),ls='-',lw=4,color='red',label='Total')
ax1.legend(loc='best',fontsize=15)
plt.tight_layout()

fig.savefig('protonDist.pdf')

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
ax1.set_xlim(x_ge)
ax1.set_ylim(y_eq)

ax1.set_xlabel(r'$\mathrm{Log}~\gamma_e$',fontsize=19)
ax1.set_ylabel(r'$\mathrm{Log}(E^2 Q_\mathrm{e}(E) ~ [\mathrm{erg~cm^{-3}s^{-1}}])$',fontsize=19)

qeTot = np.zeros(nE)
for r1 in np.arange(nR//f):
    ax1.plot(log_ge[f*r1*nE:(f*r1+1)*nE],2.0*logerge[f*r1*nE:(f*r1+1)*nE] +
             logqe[f*r1*nE:(f*r1+1)*nE],color=viridis(colors[r1*f]))

for r1 in np.arange(nR):
     for e in np.arange(nE):
        qeTot[e] = qeTot[e] + np.power(10,logqe[r1*nE+e])

    
ax1.plot(log_ge[:nE],2.0*logerge[:nE]+np.log10(qeTot),ls='-',lw=4,color='red',label='Total')
ax1.legend(loc='best',fontsize=15)
plt.tight_layout()

fig.savefig('electronInj.pdf')

fig, ax1 = plt.subplots()

ax1.tick_params(axis='both',labelsize=15)
ax1.set_xlim(x_gp)
ax1.set_ylim(y_pq)
ax1.set_xlabel(r'$\mathrm{Log}~\gamma_p$',fontsize=19)
ax1.set_ylabel(r'$\mathrm{Log}(E^2 Q_\mathrm{p}(E) ~ [\mathrm{erg~cm^{-3}s^{-1}}])$',fontsize=19)

qpTot = np.zeros(nE)
qpionTot = np.zeros(nE)
qmuonTot = np.zeros(nE)
qpairTot = np.zeros(nE)

for r1 in np.arange(nR//f):
    ax1.plot(log_gp[f*r1*nE:(f*r1+1)*nE],2.0*logergp[f*r1*nE:(f*r1+1)*nE] +
             logqp[f*r1*nE:(f*r1+1)*nE],color=viridis(colors[r1*f]))

for r1 in np.arange(nR):
     for e in np.arange(nE):
        qpTot[e] = qpTot[e] + np.power(10,logqp[r1*nE+e])
        qpionTot[e] = qpionTot[e] + np.power(10,logqpion[r1*nE+e])
        qmuonTot[e] = qmuonTot[e] + np.power(10,logqmuon[r1*nE+e])
        qpairTot[e] = qpairTot[e] + np.power(10,logqpair[r1*nE+e])

ax1.plot(log_gp[:nE],2.0*logergp[:nE]+np.log10(qpTot),ls='-',lw=4,color='red',label='Total')

ax1.legend(loc='best',fontsize=15)
plt.tight_layout()

fig.savefig('protonInj.pdf')


fig, ax1 = plt.subplots()

ax1.tick_params(axis='both',labelsize=15)

logr,logQe_r = np.loadtxt('electronInjection_R.dat',unpack=True)
logr,logQp_r = np.loadtxt('protonInjection_R.dat',unpack=True)
logr,logNe_r = np.loadtxt('electronDistribution_R.dat',unpack=True)
logr,logNp_r = np.loadtxt('protonDistribution_R.dat',unpack=True)

ax1.set_xlabel(r'${\rm Log}(r/r_{\rm schw})$')
ax1.set_ylabel(r'${\rm Log}(L_{\rm inj}/{\rm erg~cm^{-3}s^{-1}})$')

ax1.plot(logr,logQe_r,ls='--',lw=2,color='r',label='Electrons')
ax1.plot(logr,logQp_r,ls='-.',lw=2,color='b',label='Protons')
ax1.legend()

fig.savefig('nonThermalInjections_R.pdf')

fig, ax1 = plt.subplots()

ax1.tick_params(axis='both',labelsize=15)

ax1.set_ylim([-6,6])
ax1.set_xlabel(r'${\rm Log}(r/r_{\rm schw})$')
ax1.set_ylabel(r'${\rm Log}(E_{\rm NT}/{\rm erg~cm^{-3}})$')

ax1.plot(logr,logNe_r,ls='--',lw=2,color='r',label='Electrons')
ax1.plot(logr,logNp_r,ls='-.',lw=2,color='b',label='Protons')
ax1.legend()

fig.savefig('nonThermalDistributions_R.pdf')










fig, ax1 = plt.subplots()

ax1.tick_params(axis='both',labelsize=15)
ax1.set_xlim([8,17])
ax1.set_ylim(y_pq)
ax1.set_xlabel(r'$\mathrm{Log}(E/{\rm eV})$',fontsize=19)
ax1.set_ylabel(r'$\mathrm{Log}(E^2 Q_\mathrm{p}(E) ~ [\mathrm{erg~cm^{-3}s^{-1}}])$',fontsize=19)

#ax1.plot(logeVp[:nE],2.0*logergp[:nE]+np.log10(qpTot),ls='-',lw=4,color='red',label='Primary protons')
#ax1.plot(logeVe[:nE],2.0*logerge[:nE]+np.log10(qeTot),ls='-',lw=4,color='cyan',label='Primary electrons')
ax1.plot(logeVpion[:nE],2.0*logergpion[:nE]+np.log10(qpionTot),ls='-.',lw=3,color='blue',label='Pions')
ax1.plot(logeVmuon[:nE],2.0*logergmuon[:nE]+np.log10(qmuonTot),ls=':',lw=3,color='green',label='Muons')
ax1.plot(logeVpair[:nE],2.0*logergpair[:nE]+np.log10(qpairTot),ls='-',lw=4,color='yellow',label='Pairs')
ax1.plot(NT_logeV,NT_logpp,ls='--',lw=2,color='purple',label='pp')
ax1.plot(NT_logeVs,NT_logSys,ls='--',lw=2,color='red',label='SyPairs')
ax1.plot(NT_logeVs,NT_logICs,ls='--',lw=2,color='orange',label='ICPairs')
ax1.plot(NT_logeV,NT_logpg,ls=':',lw=2,color='cyan',label='pg')
ax1.legend(loc='best',fontsize=10)
plt.tight_layout()

fig.savefig('protonCascades.pdf')

