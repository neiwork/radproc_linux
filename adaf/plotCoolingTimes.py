import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from matplotlib.backends.backend_pdf import *
from matplotlib import rc        # TO MANAGE MATPLOTLIB PARAMETERS"
rc('font',family='serif')
rc('text',usetex = True)

import scipy.optimize as optimization

r,log_ge,logAcce,logAdve,logDiffe,logEmaxHe,logSye,logIC,logBr = np.loadtxt('electronCoolingTimes.dat',unpack=True,skiprows=1)
r,log_gp,logAccp,logAdvp,logDiffp,logEmaxHp,logSyp,logICp_Th,logICp,logpp,logpg,logBH = np.loadtxt('protonCoolingTimes.dat',unpack=True,skiprows=1)

x_ge = [log_ge[0]+.5,log_ge[-1]]
x_gp = [log_gp[0]+.5,log_gp[-1]]
y_e = [-15,10]
y_p = [-10,15]

aux = r[0]
nE = 0
while (aux == r[0]):
    nE += 1
    aux = r[nE]

with PdfPages('electronCoolingTimes.pdf') as pdf:
    for r1 in range(np.size(r)//nE):
        fig, ax1 = plt.subplots()
        ax1.tick_params(axis='both',labelsize=12)
        ax1.set_xlim(x_ge)
        ax1.set_ylim(y_e)

        ax1.set_title(r'$R=$'+str(np.int8(r[r1*nE]))+r'$R_{\rm schw}$',fontsize=15)
        ax1.set_xlabel(r'$\mathrm{Log}(\gamma_e)$',fontsize=13)
        ax1.set_ylabel(r'$\mathrm{Log}(t ~ [\mathrm{s}])$',fontsize=13)
    
        logTot = -np.log10(np.power(10,-logSye) + np.power(10,-logBr) + np.power(10,-logIC))

        ax1.plot(log_ge[r1*nE:(r1+1)*nE],logSye[r1*nE:(r1+1)*nE],lw=1.5,ls='dashed',c='blue',label='Sy')
        ax1.plot(log_ge[r1*nE:(r1+1)*nE],logAdve[r1*nE:(r1+1)*nE],lw=1.5,ls='dashdot',c='g',label='Adv')
        ax1.plot(log_ge[r1*nE:(r1+1)*nE],logDiffe[r1*nE:(r1+1)*nE],lw=1.5,ls='dotted',c='orange',label='Diff')
        ax1.plot(log_ge[r1*nE:(r1+1)*nE],logAcce[r1*nE:(r1+1)*nE],lw=2.0,ls='solid',c='k',label='Acc')
        ax1.plot(log_ge[r1*nE:(r1+1)*nE],logBr[r1*nE:(r1+1)*nE],lw=1.5,ls=(0,(3,1,1,1)),c='purple',label='Br')
        ax1.plot(log_ge[r1*nE:(r1+1)*nE],logIC[r1*nE:(r1+1)*nE],lw=2,ls=(0,(5,1)),c='magenta',label='IC')
        ax1.plot(log_ge[r1*nE:(r1+1)*nE],logTot[r1*nE:(r1+1)*nE],lw=4,ls='solid',c='r',label='Total')
        ax1.axvline(logEmaxHe[0])
        ax1.legend(loc='best',fontsize=8)
        plt.tight_layout()
        pdf.savefig(fig)


with PdfPages('protonCoolingTimes.pdf') as pdf:
    for r1 in range(np.size(r)//nE):
        fig, ax1 = plt.subplots()
        ax1.tick_params(axis='both',labelsize=12)
        ax1.set_xlim(x_gp)
        ax1.set_ylim(y_p)

        ax1.set_title(r'$R=$'+str(np.int8(r[r1*nE]))+r'$R_{\rm schw}$',fontsize=15)
        ax1.set_xlabel(r'$\mathrm{Log}(\gamma_p)$',fontsize=13)
        ax1.set_ylabel(r'$\mathrm{Log}(t ~ [\mathrm{s}])$',fontsize=13)
    
        logph = -np.log10(np.power(10,-logpg)+np.power(10,-logBH))
        logTot = -np.log10(np.power(10,-logSyp) + np.power(10,-logpp) + np.power(10,-logph))
        ax1.plot(log_gp[r1*nE:(r1+1)*nE],logSyp[r1*nE:(r1+1)*nE],lw=1.5,ls='dashed',c='blue',label='Sy')
        ax1.plot(log_gp[r1*nE:(r1+1)*nE],logICp_Th[r1*nE:(r1+1)*nE],lw=1.5,ls='dashed',c='k',label='ICth')
        ax1.plot(log_gp[r1*nE:(r1+1)*nE],logICp[r1*nE:(r1+1)*nE],lw=1.5,ls='dashed',c='purple',label='IC')
        ax1.plot(log_gp[r1*nE:(r1+1)*nE],logAdvp[r1*nE:(r1+1)*nE],lw=1.5,ls='dashdot',c='g',label='Adv')
        ax1.plot(log_gp[r1*nE:(r1+1)*nE],logDiffp[r1*nE:(r1+1)*nE],lw=1.5,ls='dotted',c='orange',label='Diff')
        ax1.plot(log_gp[r1*nE:(r1+1)*nE],logAccp[r1*nE:(r1+1)*nE],lw=2.0,ls='solid',c='k',label='Acc')
        ax1.plot(log_gp[r1*nE:(r1+1)*nE],logpp[r1*nE:(r1+1)*nE],lw=1.5,ls=(0,(3,1,1,1)),c='purple',label='pp')
        ax1.plot(log_gp[r1*nE:(r1+1)*nE],logph[r1*nE:(r1+1)*nE],lw=2,ls=(0,(5,1)),c='magenta',label=r'p$\gamma$')
        ax1.plot(log_gp[r1*nE:(r1+1)*nE],logTot[r1*nE:(r1+1)*nE],lw=4,ls='solid',c='r',label='Total')
        ax1.axvline(logEmaxHe[0])
        ax1.legend(loc='best',fontsize=8)
        plt.tight_layout()
        pdf.savefig(fig)
