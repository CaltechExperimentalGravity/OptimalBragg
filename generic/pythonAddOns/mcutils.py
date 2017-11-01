#Collection of functions that are useful for doing some MC analysis of coating performance

import numpy as np
import scipy.io as scio
from scipy.interpolate import interp1d, PchipInterpolator
import sys
#Import the pygwinc functions
sys.path.append('/ligo/svncommon/pygwinc')
from gwinc import *

#Some function definitions
def multidiel1(n,L,lamb,theta=0,pol='te'):
    M = len(n)-2                # number of slabs  
    if M==0:
        L = np.array([])
    theta = theta * np.pi / 180.
    costh = np.conj(np.sqrt(np.conj(1 - (n[0] * np.sin(theta) / n)**2)))
    if (pol=='te' or pol=='TE'):
        nT = n * costh
    else:
        nT = n / costh;
    if M>0:
        L = L * costh[1:M+1]
    r = -np.diff(nT) / (np.diff(nT) + 2*nT[0:M+1])
    Gamma1 = r[M] * np.ones(len(np.array([lamb])))
    for i in range(M-1,-1,-1):
        delta = 2*np.pi*L[i]/lamb                  
        z = np.exp(-2*1j*delta)
        Gamma1 = (r[i] + Gamma1*z) / (1 + r[i]*Gamma1*z)
    Z1 = (1 + Gamma1) / (1 - Gamma1)
    return Gamma1,Z1

def op2phys(L,n):
    phys = L/n
    return phys

def lnprob(x, mu, icov):
    diff = x-mu
    return -np.dot(diff,np.dot(icov,diff))/2.0

def surfaceField(gamm,Ei=27.46):
	sField = Ei * np.abs(1+gamm)
	return sField

def specREFL(matFileName, dispFileName, lambda_0=1064e-9, lam=np.linspace(0.4,1.6,2200), aoi=0., pol='tm'):
    data = scio.loadmat(matFileName, struct_as_record=False, squeeze_me=True)
    dispersion=scio.loadmat(dispFileName, struct_as_record=False, squeeze_me=True)
    aoi = float(np.copy(data['costOut'].aoi))
    n_IR = np.copy(data['costOut'].n_IR)
    L = np.copy(data['costOut'].L)
    no_of_stacks = len(L)/2
    Rp = np.ones(len(lam))
    Rs = np.ones(len(lam))
    Tp = np.ones(len(lam))
    Ts = np.ones(len(lam))
    nSiO2 = PchipInterpolator(dispersion['SiO2'][:,0], dispersion['SiO2'][:,1],
                        extrapolate=True)
    nTa2O5 = PchipInterpolator(dispersion['Ta2O5'][:,0], dispersion['Ta2O5'][:,1],
                    extrapolate=True)
    n1_IR = nSiO2(1064)
    n2_IR = nTa2O5(1064)

    for ii in range(len(lam)):
        n1 = nSiO2(lam[ii]*1064.)
        n2 = nTa2O5(lam[ii]*1064.)
        nb = 1.
        n_c = np.copy(n_IR)
        n_c[1:-1:2] = n1
        n_c[2::2] = n2
        n_c[-1] = n1
        L_temp = np.copy(L)
        L_temp[1::2] *= (n2/n2_IR)
        L_temp[2::2] *= (n1/n1_IR)
        [Gammap, Z1] = multidiel1(n_c, L_temp, lam[ii],aoi,pol)
        Rp[ii] = np.abs(Gammap)**2
        Tp[ii] = 1 - Rp[ii]
    return Rp, Tp

