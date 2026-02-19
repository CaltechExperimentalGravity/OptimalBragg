# Collection of functions that are useful for doing some
# MC analysis of coating performance

import numpy as np
import numba
import scipy.io as scio
from scipy.interpolate import interp1d, PchipInterpolator
import yaml
import sys


@numba.njit(cache=True)
def _transfer_matrix_loop(r, L_adj, lamb, M):
    """Numba-accelerated inner loop of the transfer matrix method.

    Parameters
    ----------
    r : ndarray of complex128
        Fresnel reflection coefficients at each interface.
    L_adj : ndarray of complex128
        Angle-adjusted optical thicknesses of each layer.
    lamb : ndarray of float64
        Wavelength(s) normalized to design wavelength.
    M : int
        Number of slabs.

    Returns
    -------
    Gamma1 : ndarray of complex128
        Reflection coefficient at each wavelength.
    """
    N_wl = lamb.shape[0]
    Gamma1 = np.empty(N_wl, dtype=np.complex128)
    for j in range(N_wl):
        Gamma1[j] = r[M]
    for i in range(M - 1, -1, -1):
        for j in range(N_wl):
            delta = 2.0 * np.pi * L_adj[i] / lamb[j]
            z = np.exp(-2.0j * delta)
            Gamma1[j] = (r[i] + Gamma1[j] * z) / (1.0 + r[i] * Gamma1[j] * z)
    return Gamma1


def multidiel1(n, L, lamb, theta=0, pol='te'):
    """Compute reflection coefficient of a multilayer dielectric stack.

    Uses the transfer matrix method (TMM) to calculate the complex
    reflection coefficient for a stack of dielectric layers. The inner
    loop is Numba JIT-compiled for performance.

    Parameters
    ----------
    n : array_like
        Refractive indices, including the incident and transmitted media.
        Ordered from incident medium to substrate. Length ``M + 2`` where
        ``M`` is the number of layers.
    L : array_like
        Optical thicknesses of each layer, in units of the design
        wavelength. Length ``M``.
    lamb : float or array_like
        Wavelength(s) at which to evaluate, normalized to design
        wavelength. For example, ``1.0`` is the design wavelength.
    theta : float, optional
        Angle of incidence in degrees. Default is 0 (normal incidence).
    pol : {'te', 'tm'}, optional
        Polarization. ``'te'`` is s-polarization, ``'tm'`` is
        p-polarization. Default is ``'te'``.

    Returns
    -------
    Gamma1 : complex or ndarray
        Complex reflection coefficient at each wavelength.
    Z1 : complex or ndarray
        Complex impedance at the first interface (correct only when the
        incident medium is vacuum).

    Examples
    --------
    >>> r, _ = multidiel1(n, L, [1.0, 0.5], 45.3, 'te')

    Evaluates the amplitude reflectivity at the design wavelength and the
    second harmonic, at 45.3 degrees for TE polarization.

    References
    ----------
    .. [1] http://eceweb1.rutgers.edu/~orfanidi/ewa/
    """
    # Number of slabs
    M = len(n) - 2

    # AoI in radians
    theta = theta * np.pi / 180

    # Cosine of theta; projection for oblique incidence
    costh = np.conj(np.sqrt(np.conj(1 - (n[0] * np.sin(theta) / n)**2)))

    # Polarization projection
    if (pol=='te' or pol=='TE'):
        nT = n * costh
    else:
        nT = n / costh

    # Thicknesses
    if M == 0:
        L = np.array([])
    else:
        L = L * costh[1:M+1]

    # Reflectivity
    r = -np.diff(nT) / (np.diff(nT) + 2*nT[0:M+1])

    # Ensure lamb is a 1-D float64 array for the JIT kernel
    lamb_arr = np.atleast_1d(np.asarray(lamb, dtype=np.float64))

    # Ensure r and L are complex128 arrays for the JIT kernel
    r_c = np.asarray(r, dtype=np.complex128)
    L_c = np.asarray(L, dtype=np.complex128)

    Gamma1 = _transfer_matrix_loop(r_c, L_c, lamb_arr, M)
    Z1 = (1 + Gamma1) / (1 - Gamma1)
    return Gamma1, Z1

def op2phys(L, n):
    """Convert optical thicknesses to physical thicknesses.

    Parameters
    ----------
    L : array_like
        Optical thicknesses (in units of design wavelength).
    n : array_like
        Refractive indices. Must have the same length as ``L``.

    Returns
    -------
    phys : ndarray
        Physical thicknesses.

    Raises
    ------
    ValueError
        If ``L`` and ``n`` have different lengths.
    """
    if len(L) != len(n):
        raise ValueError(f'L (dim {len(L)}) and n (dim {len(n)}) must have the same dimension.')
    phys = L / n
    return phys

def lnprob(x, mu, icov):
    """Log-probability of a multivariate Gaussian.

    Used as the target distribution for the ``emcee`` ensemble sampler
    in Monte Carlo sensitivity analysis.

    Parameters
    ----------
    x : array_like
        Sample point.
    mu : array_like
        Mean of the Gaussian.
    icov : ndarray
        Inverse covariance matrix.

    Returns
    -------
    logp : float
        Log-probability at ``x``.
    """
    diff = x - mu
    return -np.dot(diff, np.dot(icov, diff)) / 2.0

def surfaceField(gamm, Ei=27.46):
    """Compute the surface electric field of a dielectric coating.

    Parameters
    ----------
    gamm : complex
        Amplitude reflectivity at the interface of incidence.
    Ei : float, optional
        Incident electric field in V/m. Default is 27.46 V/m,
        corresponding to an intensity of 1 W/m^2.

    Returns
    -------
    sField : float
        Surface electric field magnitude in V/m.
    """
    sField = Ei * np.abs(1+gamm)
    return sField

def specREFL(layers, dispFileName, lambda_0=1064e-9,
                 lam=np.linspace(0.4,1.6,2200), aoi=0., pol='tm'):
    """Compute spectral power reflectivity with dispersion.

    Evaluates reflectivity across a wavelength range, accounting for the
    wavelength dependence of the refractive indices via Pchip
    interpolation of measured dispersion data.

    Parameters
    ----------
    layers : str or array_like
        Path to a MATLAB ``.mat`` output file, or a 1-D array of optical
        thicknesses.
    dispFileName : str
        Path to a ``.mat`` file containing dispersion data with keys
        ``'SiO2'`` and ``'Ta2O5'``.
    lambda_0 : float, optional
        Design wavelength in meters. Default is 1064 nm.
    lam : array_like, optional
        Wavelengths at which to evaluate, normalized to ``lambda_0``.
        Default spans 400--1600 nm.
    aoi : float, optional
        Angle of incidence in degrees. Default is 0.
    pol : {'te', 'tm'}, optional
        Polarization. Default is ``'tm'``.

    Returns
    -------
    Rp : ndarray
        Power reflectivity at each wavelength.
    Tp : ndarray
        Power transmissivity at each wavelength.
    ll : ndarray
        Wavelengths in meters.
    """
    if isinstance(layers, str):
        data = scio.loadmat(layers, struct_as_record=False, squeeze_me=True)
        aoi = float(np.copy(data['costOut'].aoi))
        #n_IR = np.copy(data['costOut'].n_IR)
        L = np.copy(data['costOut'].L)
    elif isinstance(layers, np.ndarray):
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
    """Compute normalized E-field squared as a function of coating depth.

    Follows the derivation in Arnon and Baumeister (1980) [1]_.

    Parameters
    ----------
    L : array_like
        Physical thickness of each coating layer in meters.
    n : array_like
        Refractive indices, including the incident and substrate media.
        Length is ``len(L) + 2``.
    lam : float, optional
        Wavelength in meters. Default is 1064 nm.
    theta : float, optional
        Angle of incidence in degrees. Default is 0.
    pol : {'s', 'p'}, optional
        Polarization. Default is ``'s'``.
    nPts : int, optional
        Number of sample points per layer. Default is 30.

    Returns
    -------
    z : ndarray
        Depth positions in meters (length ``len(L) * nPts``).
    Enorm : ndarray
        Electric field squared, normalized to the surface value
        (length ``len(L) * nPts``).

    References
    ----------
    .. [AB1980] Arnon and Baumeister, "Electric field distribution and the
       reduction of laser damage in dielectrics," Appl. Opt. 19,
       1853--1855 (1980).
    """
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
    """Load optimization parameters from a YAML config file.

    Parameters
    ----------
    paramFile : str
        Path to the YAML parameter file.

    Returns
    -------
    pars : dict
        Parsed parameters with keys ``'costs'`` and ``'misc'``.
    """
    with open(paramFile,'r') as f:
        params = yaml.safe_load(f)
    return(params)

def calcAbsorption(Esq, L, nPts, alphaOdd, alphaEven):
    """Compute integrated absorption from an E-field profile.

    Designed to work with the output of :func:`fieldDepth`.

    Parameters
    ----------
    Esq : array_like
        Normalized electric field squared as a function of depth.
        Length must equal ``len(L) * nPts``.
    L : array_like
        Physical thickness of each coating layer in meters.
    nPts : int
        Number of sample points per layer (must match ``Esq``).
    alphaOdd : float
        Absorption coefficient of odd-numbered layers (1st, 3rd, ...)
        in m^-1.
    alphaEven : float
        Absorption coefficient of even-numbered layers (2nd, 4th, ...)
        in m^-1.

    Returns
    -------
    absorp : float
        Total coating absorption in ppm.

    Raises
    ------
    ValueError
        If ``len(Esq) != len(L) * nPts``.
    """
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
    """Compute refractive index from Sellmeier dispersion coefficients.

    Parameters
    ----------
    B : list or array_like, optional
        Sellmeier B coefficients. Default is fused silica.
    C : list or array_like, optional
        Sellmeier C coefficients in um^2. Default is fused silica.
    lam : float or array_like, optional
        Wavelength(s) in meters. Default is 1064 nm.

    Returns
    -------
    n : float or ndarray
        Refractive index at each wavelength.

    Raises
    ------
    ValueError
        If ``B`` and ``C`` have different lengths.

    References
    ----------
    .. [Sellmeier] https://en.wikipedia.org/wiki/Sellmeier_equation
    """
    if len(B) != len(C):
        raise ValueError(f'The Sellmeier coefficients A (len {len(B)}) and B (len {len(C)}) must have the same number of elements.')
    n = 1
    ll = lam*1e6 # Sellmeier coefficients are quoted for wavelength in microns
    for bb,cc in zip(B,C):
        n += bb * ll**2 / (ll**2 - cc)
    n = np.sqrt(n)
    return(n)
