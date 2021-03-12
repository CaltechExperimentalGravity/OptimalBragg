# Collection of functions that are useful for doing some
# MC analysis of coating performance

import numpy as np
import scipy.io as scio
from scipy.interpolate import interp1d, PchipInterpolator
import yaml
import gwinc
import sys

#Some function definitions
def multidiel1(n, L, lamb, theta=0, pol='te'):
    '''
    Calculates reflectivity and compelx impedance of a dielectric stack.

    Parameters:
    -----------
    n: array_like
        Array of refractive indices, including the incident and 
        transmitted media. Ordered from incident medium to 
        transmitted medium.
    L: array_like
        Array of optical thicknesses comprising the dielectric 
        stack, ordered from incident medium to transmitted medium.
        Should have 2 fewer elements than n.
    lamb: float or array_like
        Wavelength(s) at which the reflectivity is to be evaluated, 
        in units of some central (design) wavelength. 
    theta: float
        Angle of incidence in degrees. 
        Defaults to 0 degrees (normal incidence)
    pol: str, 'te' or 'tm'
        Polarization at which reflectivity is to be evaluated. 
        Defaults to 'te' (s-polarization)

    Returns:
    --------
    Gamma1: complex
        Amplitude reflectivity for the input dielectric stack.
    Z1:
        Complex impedance at the interface (only correct if incident medium
        is vacuum).
    Example usage:
    --------------
    r_p, _ = multidiel1(n, L, [1.0, 0.5], 45.3, 'te')
    evaluates the amplitude reflectivity for the dielectric stack 
    specified by n and L (which are wavelength dependent in general), 
    at a design wavelength and the second harmonic wavelength, 
    at an angle of incidence of 45.3 degrees for 'te' polarized (s-pol) light.

    References:
    -----------
    [1]: http://eceweb1.rutgers.edu/~orfanidi/ewa/
    '''
    M = len(n) - 2                # number of slabs
    if M == 0:
        L = np.array([])
    theta = theta * np.pi / 180
    costh = np.conj(np.sqrt(np.conj(1 - (n[0] * np.sin(theta) / n)**2)))
    if (pol=='te' or pol=='TE'):
        nT = n * costh
    else:
        nT = n / costh
    if M > 0:
        L = L * costh[1:M+1]
    r = -np.diff(nT) / (np.diff(nT) + 2*nT[0:M+1])
    Gamma1 = r[M] * np.ones(len(np.array([lamb])))
    for i in range(M-1,-1,-1):
        delta = 2*np.pi*L[i]/lamb
        z = np.exp(-2*1j*delta)
        Gamma1 = (r[i] + Gamma1*z) / (1 + r[i]*Gamma1*z)
    Z1 = (1 + Gamma1) / (1 - Gamma1)
    return Gamma1,Z1

def op2phys(L, n):
    '''
    Converts optical lengths to physical lengths.

    Parameters:
    -----------
    L: array_like
        Array of optical thicknesses
    n: array_like
        Array of refractive indices, of the same length as L.

    Returns:
    --------
    phys: array_like
        Array of physical thicknesses for the dielectric stack 
        specified by L and n.
    '''
    if len(L) != len(n):
        raise ValueError(f'L (dim {len(L)}) and n (dim {len(n)}) must have the same dimension.')
        sys.exit()
    phys = L / n
    return phys

def lnprob(x, mu, icov):
    diff = x - mu
    return -np.dot(diff, np.dot(icov, diff)) / 2.0

def surfaceField(gamm,Ei=27.46):
    '''
    Calculates the surface electric field for a dielectric coating 
    with given amplitude reflectivity at the interface of incidence, 
    for an incident electric field.
    Parameters:
    -----------
    gamm: float
        Amplitude reflectivity of coating at interface of incidence.
    Ei: float
        Incident electric field, in V/m. Defaults to 27.46 V/m, 
        corresponding to an intensity of 1 W/m^2.
    '''
    sField = Ei * np.abs(1+gamm)
    return sField

def specREFL(layers, dispFileName, lambda_0=1064e-9,
                 lam=np.linspace(0.4,1.6,2200), aoi=0., pol='tm'):
    '''
    Computes the spectral (power) reflectivity for coating output 
    from the optimization code.

    Parameters:
    -----------
    layers: str or array_like
        Path to the MATLAB output file from the coating optimization code.
        Alternatively, this can be a 1D array of optical thicknesses.
    dispFileName: str
        Path to a .mat file containing the dispersion data for the coating.
    lambda_0: float
        Design (central) wavelength for the coating optimization, in meters. 
        Defaults to 1064nm.
    lam: array_like
        Array of wavelengths at whihc to evaluate the reflectivity. 
        Defaults to [400nm 1600nm].
    aoi: float
        Angle of incidence at which to evaluate reflectivity. 
        Defaults to 0 (normal incidence)
    pol: str
        Polarization to evaluate reflectivity. 
        Defaults to 'tm' (p-polarization).

    Returns:
    --------
    Rp: array_like
        Reflectivity of coating 
    Tp: array_like
        Transmissivity of coating
    ll: array_like
        Array of wavelengths at which the reflectivity was evaluated.
    '''
    if type(layers)==str:
        data = scio.loadmat(layers, struct_as_record=False, squeeze_me=True)
        aoi = float(np.copy(data['costOut'].aoi))
        #n_IR = np.copy(data['costOut'].n_IR)
        L = np.copy(data['costOut'].L)
    elif type(layers)==np.ndarray:
        L = np.copy(layers)
    else:
        print('{} is an unknown type of layer specifications. Please provide \
                a path to a .mat file or a 1D array of optical thicknesses'.format(layers))
        return
    dispersion=scio.loadmat(dispFileName, struct_as_record=False, squeeze_me=True)
    nSiO2 = PchipInterpolator(dispersion['SiO2'][:,0], dispersion['SiO2'][:,1],
                        extrapolate=True)
    nTa2O5 = PchipInterpolator(dispersion['Ta2O5'][:,0], dispersion['Ta2O5'][:,1],
                    extrapolate=True)
    no_of_stacks = len(L)/2
    Rp = np.ones(len(lam))
    Rs = np.ones(len(lam))
    Tp = np.ones(len(lam))
    Ts = np.ones(len(lam))
    n1_IR = nSiO2(1064)
    n2_IR = nTa2O5(1064)

    for ii in range(len(lam)):
        n1 = nSiO2(lam[ii]*1064.)
        n2 = nTa2O5(lam[ii]*1064.)
        nb = 1.
        #n_c = np.copy(n_IR)
        n_c = np.ones(len(L)+2)
        n_c[1:-1:2] = n1
        n_c[2::2] = n2
        n_c[-1] = n1
        L_temp = np.copy(L)
        L_temp[1::2] *= (n2/n2_IR)
        L_temp[2::2] *= (n1/n1_IR)
        [Gammap, Z1] = multidiel1(n_c, L_temp, lam[ii],aoi,pol)
        Rp[ii] = np.abs(Gammap)**2
        Tp[ii] = 1 - Rp[ii]
    return Rp, Tp, lambda_0*lam

def fieldDepth(L, n, lam=1064e-9, theta=0, pol='s', nPts=30):
    '''
    Function that calculates "Normalized" E-Field strength squared
    as a function of penetration depth in a dielectric coating.

    Parameters:
    ------------
    L: array_like
        Array of PHYSICAL thickness of coating layers.
    n: array_like 
        Array of refractive indices.

    Returns:
    --------
    z: array_like
        Array of penetration depths at which E-field is evaluated.
    Enorm: array_like
        Electric field SQUARED normalized to that at the interface of incidence.

    Following derivation set out in Arnon and Baumeister, 1980
    https://www.osapublishing.org/ao/abstract.cfm?uri=ao-19-11-1853
    '''
    # check to see that the lengths of L and n match in some way

    
    def arrayTheta(n, theta0):
        alpha = []
        alpha.append(np.deg2rad(theta0))
        for ii in range(len(n) - 1):
            t_r = np.arcsin(n[ii] * np.sin(alpha[ii]) / n[ii+1])
            alpha.append(t_r)
        return np.array(alpha[1:]) #Note that angles are all in radians
    # Calculate the array of angles, in radians
    angles = arrayTheta(n, theta)
    qAngle = angles[-1]
    angles = angles[0:-1]
    def M_i(b_i, qq_i):
        # Eqn 2 from paper
        out = np.matrix([[np.cos(b_i), 1j*np.sin(b_i)/qq_i], [1j*np.sin(b_i)*qq_i, np.cos(b_i)]])
        return out
    def q_i(n_i, theta_i):
        if pol == 's':
            return n_i * np.cos(theta_i) # Eqn 5
        elif pol == 'p':
            return n_i / np.cos(theta_i) # Eqn 6
    def beta_i(tt_i, nn_i, hh_i):
        return 2*np.pi * np.cos(tt_i) * nn_i * hh_i / lam #Eqn 3

    # Calculate the total matrix, as per Eqn 7.
    Mtot = np.eye(2)
    for n_i, h_i, theta_i in zip(n[1:-1], L, angles):
        Mtot = Mtot * M_i(beta_i(theta_i, n_i, h_i), q_i(n_i, theta_i))

    Mtotz = Mtot
    def E0pk(Mtot):
        q0   = q_i(n[0],  theta)
        qSub = q_i(n[-1], qAngle)
        #Eqn 10
        return 0.25*(np.abs(Mtot[0,0] + Mtot[1,1]*qSub/q0)**2 +
                np.abs(Mtot[1,0]/q0/1j + Mtot[0,1]*qSub/1j)**2)
    def delta_h(bb_i, qq_i):
        #Equation 11
        return M_i(bb_i, -qq_i)

    #Initialize some arrays to store the calculated E field profile
    E_profile = np.zeros(len(L)*nPts)
    z = np.zeros(len(L)*nPts)
    Z = 0

    #Initialize the q-parameter at the rightmost interface
    qSub = q_i(n[-1], qAngle)

    for ii in range(len(L)):
        n_i = n[ii+1]
        dL = L[ii] / nPts
        theta_i = angles[ii]

        if pol=='p':
            #correction=(np.cos(theta[0])/np.cos(theta_i))**2
            correction = (np.cos(theta) / np.cos(theta_i))**2
        elif pol=='s':
            correction = 1

        for jj in range(0,nPts):
            Z += dL
            z[ii*nPts + jj] = Z
            Mtotz = delta_h(beta_i(theta_i, n_i, dL),q_i(n_i, theta_i)) * Mtotz
            E_profile[ii*nPts+jj] = correction * (np.abs(Mtotz[0,0])**2 +
                                            np.abs(qSub*Mtotz[0,1]/1j)**2)

    return z, E_profile/E0pk(Mtot)

def importParams(paramFile):
    '''
    Function to load some parameters from a yaml config file to run an optimizer.
    Parameters:
    -----------
    paramFile: str
        Path to parameter file
    Returns:
    --------
    pars: dict
        A dict from which we can access various params to set up the optimizer
    '''
    with open(paramFile,'r') as f:
        params = yaml.safe_load(f)
    return(params)

def calcAbsorption(Esq, L, nPts, alphaOdd, alphaEven):
    '''
    Function for calculating the Absorption given an electric field profile
    Made to work together with fieldDepth.
    Parameters:
    -----------
    Esq: array_like
        NORMALIZED electric field squared as a function of distance in a coating
    L: array_like         
        PHYSICAL thickness of coating layers
    nPts: int 
        # of points inside each layer at which field is to be evaluated
    alphaOdd: float
        Absorption coefficient of all odd layers [m ^-1] (top layer is assumed layer #1)
    alphaEven: float
        Absorption coefficient of all even layers [m ^-1] (top layer is assumed layer #1)
    Returns:
    -------------
    absorp: float 
        Absorption of coating [ppm]
    '''
    if len(Esq) != int(len(L)*nPts):
        raise ValueError(f'The input electric field vector length, {len(Esq)} is not consistent with the number of points requested per layer, {nPts}, and the number of layers, {len(L)}.')
    dL = L/nPts
    dz = []
    #Define the grid for rectangular integration
    for ii in range(len(L)):
        temp = np.tile(dL[ii], nPts)# Pythonic rep of repmat(dL(ii),nPts,1)
        dz = np.hstack((dz,temp)) # Pythonic implementation of vertcat
    alph = np.ones(nPts*len(L))
    alph[::2] = alphaOdd
    alph[1::2] = alphaEven
    absorp = 1e6 * np.sum(alph * dz * Esq)
    return(absorp)

def sellmeier(B=[0.696166300, 0.407942600, 0.897479400], C=[4.67914826e-3, 1.35120631e-2, 97.9340025], lam=1064e-9):
    '''
    Function to calculate dispersion using Sellmeier coefficients

    Parameters:
    ------------
    B: list or array_like
        "B" coefficients in the Sellmeier equation.
        Defaults to first 3 coefficients for Fused Silica.
    C: list or array_like
        "C" coefficients in the Sellmeier equation [um^2].
        Defaults to first 3 coefficients for Fused Silica.
    lam: float or array_like
        Value or array of wavelengths [m]. Conversion to 
        microns is done internally to the function.
        Defaults to 1064nm.

    Returns:
    ----------
    n: float or array_like
        Refractive index (same shape as lam)

    Ref:
    ----
    https://en.wikipedia.org/wiki/Sellmeier_equation#:~:text=The%20Sellmeier%20equation%20is%20an,of%20light%20in%20the%20medium.
    '''
    if len(B) != len(C):
        raise ValueError(f'The Sellmeier coefficients A (len {len(B)}) and B (len {len(C)}) must have the same number of elements.')
        sys.exit()
    n = 1
    ll = lam*1e6 # Sellmeier coefficients are quoted for wavelength in microns
    for bb,cc in zip(B,C):
        n += bb * ll**2 / (ll**2 - cc)
    n = np.sqrt(n)
    return(n)