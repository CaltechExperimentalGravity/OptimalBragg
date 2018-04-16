#Collection of functions that are useful for doing some MC analysis of coating performance

import numpy as np
import scipy.io as scio
from scipy.interpolate import interp1d, PchipInterpolator
import sys
#Import the pygwinc functions
#sys.path.append('/ligo/svncommon/pygwinc')
#from gwinc import *
from itertools import izip

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

def fieldDepth(L, n, lam=1064e-9, theta=0., pol='s',nPts=30):
	'''
	Function that calculates "Normalized" E-Field strength squared
	as a function of penetration depth in a dielectric coating.
	Input argument L is an array of PHYSICAL thickness of coating layers.
	Input argument n is an array of refractive indices.

	Following derivation set out in Arnon and Baumeister, 1980
	https://www.osapublishing.org/ao/abstract.cfm?uri=ao-19-11-1853
	'''
	def arrayTheta(n,theta0):
		alpha=[]
		alpha.append(np.deg2rad(theta0))
		for ii in range(len(n)-1):
			t_r = np.arcsin(n[ii]*np.sin(alpha[ii])/n[ii+1])
			alpha.append(t_r)
		return np.array(alpha) #Note that angles are all in radians

	#Calculate the array of angles, in radians
	angles = arrayTheta(n,theta)
	def M_i(b_i,qq_i):
		#Eqn 2 from paper
		return np.matrix([[np.cos(b_i), 1j*np.sin(b_i)/qq_i],
			[1j*np.sin(b_i)*qq_i, np.cos(b_i)]])
	def q_i(n_i, theta_i):
		if pol == 's':
			return n_i*np.cos(theta_i) #Eqn 5
		elif pol == 'p':
			return n_i / np.cos(theta_i) #Eqn 6
	def beta_i(tt_i, nn_i, hh_i):
		return 2*np.pi*np.cos(tt_i)*nn_i*hh_i/lam #Eqn 3

	#Calculate the total matrix, as per Eqn 7.
	Mtot = np.eye(2)
	for n_i, h_i, theta_i in izip(n[1:-1], L, angles):
		Mtot = Mtot * M_i(beta_i(theta_i,n_i,h_i), q_i(n_i,theta_i))

	Mtotz = Mtot

	def E0pk(Mtot):
		q0 = q_i(n[0],theta)
		qSub = q_i(n[-1],angles[-1])
		#Eqn 10
		return 0.25*(np.abs(Mtot[0,0] + Mtot[1,1]*qSub/q0)**2 +
				np.abs(Mtot[1,0]/q0/1j + Mtot[0,1]*qSub/1j)**2)
	def delta_h(bb_i, qq_i):
		#Equation 11
		return M_i(bb_i, -qq_i)

	#Initialize some arrays to store the calculated E field profile
	E_profile = np.zeros(len(L)*nPts)
	z = np.zeros(len(L)*nPts)
	Z=0.
	correction=1.

	#Initialize the q-parameter at the rightmost interface
	qSub = q_i(n[-1], angles[-1])

	for ii in range(len(L)):
		n_i = n[ii+1]
		dL = L[ii] / nPts
		theta_i = angles[ii]

		if pol=='p':
			#correction=(np.cos(theta[0])/np.cos(theta_i))**2
			correction=(np.cos(theta)/np.cos(theta_i))**2

		for jj in range(0,nPts):
			Z += dL
			z[ii*nPts + jj] = Z
			Mtotz = delta_h(beta_i(theta_i, n_i, dL),q_i(n_i, theta_i)) * Mtotz
			E_profile[ii*nPts+jj] = correction * np.abs(Mtotz[0,0])**2 + np.abs(qSub*Mtotz[0,1]/1j)**2

	return z, E_profile/E0pk(Mtot)


