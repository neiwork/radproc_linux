import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from matplotlib.backends.backend_pdf import *
plt.style.use("mystyle")

r,log_g,logDiff,logAdv,ratio = np.loadtxt('diffusion_adv.dat',unpack=True,skiprows=1)

x_ge = [log_g[0],log_g[-1]]

aux = r[0]
nE = 0
while (aux == r[0]):
    nE += 1
    aux = r[nE]

with PdfPages('diff_to_Adv.pdf') as pdf:
    for r1 in range(np.size(r)//nE):
        fig, ax1 = plt.subplots()
        ax1.tick_params(axis='both',labelsize=12)
        ax1.set_xlim(x_ge)
        ax1.set_ylim([-5,7])
        ax1.set_title(r'$R=$'+str(np.int32(r[r1*nE]))+r'$R_{\rm schw}$',fontsize=15)
        ax1.set_xlabel(r'$\mathrm{Log}(\gamma_p)$',fontsize=13)
        ax1.set_ylabel(r'$\mathrm{Log}(l ~ [\mathrm{cm}])$',fontsize=13)
    

        ax1.plot(log_g[r1*nE:(r1+1)*nE],logDiff[r1*nE:(r1+1)*nE],lw=1.5,ls='-',c='blue',label='Diff')
        ax1.plot(log_g[r1*nE:(r1+1)*nE],logAdv[r1*nE:(r1+1)*nE],lw=1.5,ls='dashed',c='g',label='Adv')
        ax1.legend(loc='best',fontsize=8)
        plt.tight_layout()
        pdf.savefig(fig)

