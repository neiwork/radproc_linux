#!/home/emgutierrez/anaconda3/bin/python

import numpy as np
import astropy.constants as const
import matplotlib.pyplot as plt
import matplotlib
from constants_cgs import *
import scipy.interpolate as interp
import scipy.integrate as integ
import scipy.special as sp

eMMW = 1.14
iMMW = 1.23

# READ THE ADAF PARAMETERS
logr,logTi,logTe,logv = np.loadtxt('adafFile.txt',skiprows=1,unpack=True)
massBH,accRateOut,sWind,beta,alpha,jAngMom,delta = np.loadtxt('adafParameters.txt')

accRateOut *= 1.39e18 * massBH
massBH *= solarMass
schwRadius = 2.0*gravConstant*massBH/cLight2

logTifun = interp.interp1d(logr,logTi,kind='quadratic',fill_value='extrapolate')
logTefun = interp.interp1d(logr,logTe,kind='quadratic',fill_value='extrapolate')
logvfun = interp.interp1d(logr,logv,kind='quadratic',fill_value='extrapolate')

def electronTemp(r):
    logrr = np.log(r/schwRadius)
    return np.exp(logTefun(logrr))*eMMW

def ionTemp(r):
    logrr = np.log(r/schwRadius)
    return np.exp(logTifun(logrr))*iMMW

def radialVel(r):
    logrr = np.log(r/schwRadius)
    rg = schwRadius/2.0
    r_rg = r/rg
    return -np.exp(logvfun(logrr)) / np.where(r_rg <= 30.0, 0.93*np.exp(2.13/r_rg), 1.0)

def integrand(r):
    return 1.0/radialVel(r)

def accretionTime(r):
    return np.vectorize(integ.quad)(integrand,r,schwRadius)[0]
    
def keplAngVel(r):
    return np.sqrt(gravConstant*massBH/r) / (r-schwRadius)

def sqrdSoundVel(r):
    return boltzmann / (beta*atomicMassUnit) * (ionTemp(r)/iMMW + electronTemp(r)/eMMW)

def height_fun(r):
    return np.sqrt(sqrdSoundVel(r))/keplAngVel(r)

def angVel(r):
    return alpha*sqrdSoundVel(r)/(-radialVel(r)*r) + (jAngMom*schwRadius*cLight)/(r*r)

def accRate(r):
    rOut = np.exp(logr[-1])*schwRadius
    return accRateOut * np.power(r/rOut,sWind)

def massDensity(r):
    return accRate(r) / (4.0*np.pi*r*height_fun(r)*np.sqrt(0.5*np.pi)*(-radialVel(r)))

def magField(r):
    return np.sqrt(8.0*np.pi*(1.0-beta)*massDensity(r)*sqrdSoundVel(r))

def eDens(r):
    return massDensity(r)/(atomicMassUnit*eMMW)

def iDens(r):
    return massDensity(r)/(atomicMassUnit*iMMW)

def tRelax_e(r):
    theta = boltzmann*electronTemp(r)/electronRestEnergy
    n_e = eDens(r)
    return 4.0*np.sqrt(np.pi)/(20.0*n_e*thomson*cLight) * np.power(theta,1.5)

def tRelax_i(r):
    theta = boltzmann*ionTemp(r)/protonRestEnergy
    n_p = iDens(r)
    return 4.0*np.sqrt(np.pi)/(20.0*n_p*thomson*cLight) * np.square(protonMass/electronMass) * np.power(theta,1.5)

def tDiss(r):
    return 1.0/angVel(r)/alpha

def tRelax_e_rel(E,r):
    n_e = eDens(r)
    theta = boltzmann*electronTemp(r)/electronRestEnergy
    k1 = sp.k1(1.0/theta)
    k2 = sp.kn(2,1.0/theta)
    g = E/electronRestEnergy
    factor = np.power(np.abs(k1/k2-1.0/g),-1)
    return 2.0/3.0 * g /(n_e*thomson*cLight*(20.0+9.0/16.0-np.log(np.sqrt(2.0)))) * factor

def integrand2(log_gg,g):
    theta = boltzmann*electronTemp(10*schwRadius)/electronRestEnergy
    gg = np.exp(log_gg)
    uee = (g/gg + gg/g)/(2.0*theta)
    return gg * np.exp(-uee)*(theta*(1.0+2.0*uee)-g)

def tRelax_e_rel_exact(E,r):
    n_e = eDens(r)
    theta = boltzmann*electronTemp(r)/electronRestEnergy
    k1 = sp.k1(1.0/theta)
    k2 = sp.kn(2,1.0/theta)
    g = E/electronRestEnergy
    factor = np.power(np.abs(np.vectorize(integ.quad)(integrand2,0.0,20.0,args=(g))[0]),-1)
    return 4.0/3.0 * k2 * g*g*g /(n_e*thomson*cLight*(20.0+9.0/16.0-np.log(np.sqrt(2.0)))) * factor

def tAcc(E,r,etaAcc):
    return etaAcc*E/(electronCharge*magField(r)*cLight)

def tDSA(E,r,etaAcc,etag):
    vshock = -radialVel(r)
    return etaAcc*etag*cLight*E/(3.0*electronCharge*magField(r)*vshock*vshock)

def tMR_e(E,r):
    g = E / electronRestEnergy
    vA = magField(r)/np.sqrt(4.0*np.pi*massDensity(r))
    return 7.7e-6 * height_fun(r)/vA * np.power(electronMass/protonMass,0.9) * np.power(g,0.4)

def tMR_p(E,r):
    vA = magField(r)/np.sqrt(4.0*np.pi*massDensity(r))
    g = E / protonRestEnergy
    return 7.7e-6 * height_fun(r)/vA * np.power(g,0.4)

def tSDA_e(E,r,zeda,q):
    g = E / electronRestEnergy
    rL = E/(electronCharge*magField(r))
    vA = magField(r)/np.sqrt(4.0*np.pi*massDensity(r))
    return 1.0/zeda * np.power(vA/cLight,-2) * height_fun(r)/cLight * np.power(rL*g/height_fun(r),2-q)

def tSDA_p(E,r,zeda,q):
    g = E / protonRestEnergy
    rL = E/(electronCharge*magField(r))
    vA = magField(r)/np.sqrt(4.0*np.pi*massDensity(r))
    return 1.0/zeda * np.power(vA/cLight,-2) * height_fun(r)/cLight * np.power(rL*g/height_fun(r),2-q)

def tDiff_p(E,r,zeda,q):
    g = E / protonRestEnergy
    rL = E/(electronCharge*magField(r))
    vA = magField(r)/np.sqrt(4.0*np.pi*massDensity(r))
    return zeda * height_fun(r)/cLight * np.power(rL*g/height_fun(r),q-2)

def tDiff_e(E,r,zeda,q):
    g = E / electronRestEnergy
    rL = E/(electronCharge*magField(r))
    vA = magField(r)/np.sqrt(4.0*np.pi*massDensity(r))
    return zeda * height_fun(r)/cLight * np.power(rL*g/height_fun(r),q-2)

def tSync(E,r):
    g = E / electronRestEnergy
    Umag = magField(r)*magField(r)/(8.0*np.pi)
    dEdt = 4.0/3.0 * thomson * cLight * Umag * g*g
    return E/dEdt

rVector = np.logspace(0,2,100)*schwRadius

plt.xlabel('Log10(r/rS)')
plt.ylabel('Log10(t/s)')
plt.plot(np.log10(rVector/schwRadius),np.log10(tDiss(rVector)),color='b',label='tdiss')
plt.plot(np.log10(rVector/schwRadius),np.log10(tRelax_i(rVector)),color='r',label='tRel_i')
plt.plot(np.log10(rVector/schwRadius),np.log10(tRelax_e(rVector)),color='g',label='tRel_e')
plt.plot(np.log10(rVector/schwRadius),np.log10(accretionTime(rVector)),color='purple',label='tAcc_exact')
plt.plot(np.log10(rVector/schwRadius),np.log10(rVector/(-radialVel(rVector))),color='k',label='tAcc_approx')

plt.legend()
plt.savefig('timescales.pdf')
