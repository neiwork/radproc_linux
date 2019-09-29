import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from matplotlib import rc        # TO MANAGE MATPLOTLIB PARAMETERS"
rc('font',family='serif')
rc('text',usetex = True)

logeVe,logr,logAcce,logtcelle,logAdve,logDiffe,logEmaxHe,logSye,logIC,logBr = np.loadtxt('neutronLosses.txt',unpack=True,skiprows=1)

x_eVe = [logeVe[0]+.5,12]
x_eVp = [logeVp[0]+.5,logeVp[-1]]
y_p = [-10,15]
y_e = [-10,15]

nR = 30
f = 5
nE = 100
colors = np.arange(nR)/nR


Cell = np.arange(nR//f)
for r1 in np.arange(nR//f):
    fig, ax1 = plt.subplots()
    ax1.tick_params(axis='both',labelsize=12)
    ax1.set_xlim(x_eVe)
    ax1.set_ylim(y_e)

    ax1.set_xlabel(r'$\mathrm{Log}(E/\mathrm{eV})$',fontsize=13)
    ax1.set_ylabel(r'$\mathrm{Log}(t^{-1} ~ [\mathrm{s}^{-1}])$',fontsize=13)
    
    ax1.plot(logeVe[f*r1*nE:(f*r1+1)*nE],logtcelle[f*r1*nE:(f*r1+1)*nE],label='tEscape')
    ax1.plot(logeVe[f*r1*nE:(f*r1+1)*nE],logIC[f*r1*nE:(f*r1+1)*nE],label='pp')
    ax1.plot(logeVe[f*r1*nE:(f*r1+1)*nE],logBr[f*r1*nE:(f*r1+1)*nE],label='p$\gamma$')
    ax1.legend(loc='best',fontsize=8)
    fig.savefig('neutronLosses_'+str(Cell[r1]+1)+'.pdf')

for r1 in np.arange(nR//f):
    
    fig, ax1 = plt.subplots()
    ax1.tick_params(axis='both',labelsize=12)
    ax1.set_xlim(x_eVp)
    ax1.set_ylim(y_p)

    ax1.set_xlabel(r'$\mathrm{Log}(E/\mathrm{eV})$',fontsize=13)
    ax1.set_ylabel(r'$\mathrm{Log}(t^{-1} ~ [\mathrm{s}^{-1}])$',fontsize=13)
    
    ax1.plot(logeVp[f*r1*nE:(f*r1+1)*nE],logSyp[f*r1*nE:(f*r1+1)*nE],label='Sy')
    ax1.plot(logeVp[f*r1*nE:(f*r1+1)*nE],logAdvp[f*r1*nE:(f*r1+1)*nE],label='Adv')
    ax1.plot(logeVp[f*r1*nE:(f*r1+1)*nE],logtcellp[f*r1*nE:(f*r1+1)*nE],label='tCell')
    ax1.plot(logeVp[f*r1*nE:(f*r1+1)*nE],logDiffp[f*r1*nE:(f*r1+1)*nE],label='Diff')
    ax1.axvline(logEmaxHp[0],label='Hillas')
    ax1.plot(logeVp[f*r1*nE:(f*r1+1)*nE],logAccp[f*r1*nE:(f*r1+1)*nE],label='Acc')
    ax1.plot(logeVp[f*r1*nE:(f*r1+1)*nE],logpp[f*r1*nE:(f*r1+1)*nE],label='pp')
    ax1.plot(logeVp[f*r1*nE:(f*r1+1)*nE],logpg[f*r1*nE:(f*r1+1)*nE],label=r'p$\gamma$')
    ax1.legend(loc='best',fontsize=8)
    fig.savefig('protonLosses_'+str(Cell[r1]+1)+'.pdf')
