"""Transfer matrix method and E-field calculations for dielectric stacks.

Two APIs for the transfer matrix:

- ``multidiel1(n, L, lamb, theta, pol)`` — Numba JIT core.
  Uses optical thicknesses (fraction of design wavelength),
  normalized wavelengths, and theta in **degrees**.

- ``multilayer_diel(ns, Ls, lamb, aoi, pol)`` — Convenience wrapper.
  Uses physical thicknesses [m], absolute wavelength [m],
  and aoi in **radians**.  Converts internally and calls the JIT core.

References
----------
.. [Orfanidi] http://eceweb1.rutgers.edu/~orfanidi/ewa/ Chapter 8
.. [AB1980] Arnon & Baumeister, Appl. Opt. 19, 1853 (1980)
"""

import numpy as np
import numba

# NumPy 2.0 renamed trapz → trapezoid
_trapz = getattr(np, 'trapezoid', np.trapz)


# ── Numba JIT transfer matrix core ───────────────────────────────────────

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
        Wavelength(s) in the same units as *L_adj*.
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
    """Reflection coefficient of a multilayer dielectric stack (JIT).

    Parameters
    ----------
    n : array_like
        Refractive indices [incident, layer1, ..., layerM, substrate].
    L : array_like
        Optical thicknesses in units of the design wavelength.
    lamb : float or array_like
        Wavelength(s) normalized to the design wavelength.
    theta : float, optional
        Angle of incidence in **degrees**. Default 0.
    pol : {'te', 'tm'}, optional
        Polarization. Default ``'te'``.

    Returns
    -------
    Gamma1 : complex or ndarray
        Complex amplitude reflection coefficient.
    Z1 : complex or ndarray
        Complex impedance at the first interface.
    """
    n = np.asarray(n, dtype=complex)
    L = np.asarray(L, dtype=complex)
    M = len(n) - 2

    theta_rad = theta * np.pi / 180
    costh = np.conj(np.sqrt(np.conj(1 - (n[0] * np.sin(theta_rad) / n) ** 2)))

    if pol in ('te', 'TE'):
        nT = n * costh
    else:
        nT = n / costh

    if M == 0:
        L = np.array([], dtype=complex)
    else:
        L = L * costh[1:M + 1]

    r = -np.diff(nT) / (np.diff(nT) + 2 * nT[0:M + 1])

    lamb_arr = np.atleast_1d(np.asarray(lamb, dtype=np.float64))
    r_c = np.asarray(r, dtype=np.complex128)
    L_c = np.asarray(L, dtype=np.complex128)

    Gamma1 = _transfer_matrix_loop(r_c, L_c, lamb_arr, M)
    Z1 = (1 + Gamma1) / (1 - Gamma1)
    return Gamma1, Z1


# ── Physical-thickness API ───────────────────────────────────────────────

def multilayer_diel(ns, Ls, lamb, aoi=0, pol='te'):
    """Reflection coefficient using physical thicknesses and wavelength.

    Thin wrapper around :func:`multidiel1`; converts physical → optical
    thicknesses internally so the Numba JIT core handles the heavy work.

    Parameters
    ----------
    ns : array_like
        Refractive indices [incident, layer1, ..., layerM, substrate].
    Ls : array_like
        Physical layer thicknesses in meters.
    lamb : float or array_like
        Wavelength(s) in meters.
    aoi : float, optional
        Angle of incidence in **radians**. Default 0.
    pol : {'te', 'tm'}, optional
        Polarization. Default ``'te'``.

    Returns
    -------
    Gamma : complex or ndarray
        Complex amplitude reflection coefficient.
    Z : complex or ndarray
        Complex impedance at the first interface.
    """
    ns = np.asarray(ns, dtype=complex)
    Ls = np.asarray(Ls, dtype=float)
    # Optical path lengths in meters (same units as lamb)
    L_opt = np.real(ns[1:len(Ls) + 1]) * Ls
    theta_deg = np.degrees(aoi)
    Gamma, Z = multidiel1(ns, L_opt, lamb, theta=theta_deg, pol=pol)
    # Return scalar for scalar input
    if np.ndim(lamb) == 0:
        return Gamma.item(), Z.item()
    return Gamma, Z


# ── Stack-based convenience functions ────────────────────────────────────

def amp_refl(wavelengths, stack, **kwargs):
    """Spectral amplitude reflectivity of a stack.

    Parameters
    ----------
    wavelengths : array_like
        Wavelength(s) in meters.
    stack : dict
        Stack dict (must contain ``ns`` and ``Ls``).
    **kwargs
        Forwarded to :func:`multilayer_diel` (e.g. ``aoi``, ``pol``).

    Returns
    -------
    ndarray
        Complex amplitude reflectivity at each wavelength.
    """
    ns, Ls = stack['ns'], stack['Ls']
    wavelengths = np.atleast_1d(np.asarray(wavelengths, dtype=float))
    rr = np.empty_like(wavelengths, dtype=complex)
    for j, wl in enumerate(wavelengths):
        rr[j], _ = multilayer_diel(ns, Ls, float(wl), **kwargs)
    return rr


def refl(wavelengths, stack, **kwargs):
    """Spectral power reflectivity R = |r|^2.

    Parameters
    ----------
    wavelengths : array_like
        Wavelength(s) in meters.
    stack : dict
        Stack dict.
    **kwargs
        Forwarded to :func:`multilayer_diel`.

    Returns
    -------
    ndarray
        Power reflectivity at each wavelength.
    """
    return np.abs(amp_refl(wavelengths, stack, **kwargs)) ** 2


def trans(wavelengths, stack, **kwargs):
    """Spectral power transmissivity T = 1 - R.

    Parameters
    ----------
    wavelengths : array_like
        Wavelength(s) in meters.
    stack : dict
        Stack dict.
    **kwargs
        Forwarded to :func:`multilayer_diel`.

    Returns
    -------
    ndarray
        Power transmissivity at each wavelength.
    """
    return 1 - refl(wavelengths, stack, **kwargs)


# ── Surface E-field ──────────────────────────────────────────────────────

def surfield(rr, Ei=27.46, normalized=False):
    """Surface electric field of a dielectric coating.

    Parameters
    ----------
    rr : complex
        Amplitude reflectivity at the input interface.
    Ei : float, optional
        Incident E-field amplitude in V/m.  Default 27.46 (1 W/m^2).
    normalized : bool, optional
        If True, return ``|1 + r|`` (units of Ei). Default False.

    Returns
    -------
    float or ndarray
        Surface E-field magnitude.
    """
    if normalized:
        return np.abs(1 + rr)
    return Ei * np.abs(1 + rr)


# ── E-field depth profile ────────────────────────────────────────────────

def field_zmag(ns, Ls, lam, aoi=0, pol='s', n_pts=30):
    """Normalized E-field squared vs depth (Arnon & Baumeister 1980).

    Parameters
    ----------
    ns : array_like
        Refractive indices [incident, layer1, ..., substrate].
    Ls : array_like
        Physical thicknesses in meters.
    lam : float
        Reference wavelength in meters.
    aoi : float, optional
        Angle of incidence in **radians**. Default 0.
    pol : {'s', 'p', 'te', 'tm'}, optional
        Polarization. Default ``'s'``.
    n_pts : int, optional
        Sample points per layer. Default 30.

    Returns
    -------
    z_prof : ndarray
        Depth positions in meters.
    E_prof : ndarray
        |E|^2 normalized to the surface value.
    """
    ns = np.asarray(ns, dtype=complex)
    Ls = np.asarray(Ls, dtype=float)

    # Array of incidence angles via Snell's law
    alpha = [aoi]
    for ii in range(len(ns) - 1):
        t_r = np.arcsin(ns[ii] * np.sin(alpha[ii]) / ns[ii + 1])
        alpha.append(t_r)
    q_angle = alpha[-1]
    angles = np.array(alpha[1:-1])

    is_p = pol in ('tm', 'TM', 'p', 'P')

    def M_i(b_i, qq_i):
        return np.array([
            [np.cos(b_i), 1j * np.sin(b_i) / qq_i],
            [1j * np.sin(b_i) * qq_i, np.cos(b_i)],
        ])

    def q_i(n_i, theta_i):
        if is_p:
            return n_i / np.cos(theta_i)
        return n_i * np.cos(theta_i)

    def beta_i(tt_i, nn_i, hh_i):
        return 2 * np.pi * np.cos(tt_i) * nn_i * hh_i / lam

    # Total transfer matrix (Eq. 7)
    Mtot = np.eye(2, dtype=complex)
    for n_i, L_i, a_i in zip(ns[1:-1], Ls, angles):
        Mtot = Mtot @ M_i(beta_i(a_i, n_i, L_i), q_i(n_i, a_i))

    # Eq. 10 — |E|^2 at the surface
    q_0 = q_i(ns[0], aoi)
    q_sub = q_i(ns[-1], q_angle)
    Epeak_0 = 0.25 * (
        np.abs(Mtot[0, 0] + Mtot[1, 1] * q_sub / q_0) ** 2
        + np.abs(Mtot[1, 0] / q_0 / 1j + Mtot[0, 1] * q_sub / 1j) ** 2
    )

    def delta_h(bb_i, qq_i):
        return M_i(bb_i, -qq_i)

    E_prof = np.zeros(len(Ls) * n_pts)
    z_prof = np.zeros(len(Ls) * n_pts)
    Z = 0.0
    Mtotz = Mtot.copy()
    q_sub = q_i(ns[-1], q_angle)

    for ii in range(len(Ls)):
        n_i = ns[ii + 1]
        dL = Ls[ii] / n_pts
        a_i = angles[ii]
        corr = (np.cos(aoi) / np.cos(a_i)) ** 2 if is_p else 1.0

        for jj in range(n_pts):
            Z += dL
            z_prof[ii * n_pts + jj] = Z
            Mtotz = delta_h(beta_i(a_i, n_i, dL), q_i(n_i, a_i)) @ Mtotz
            E_prof[ii * n_pts + jj] = np.real(corr * (
                np.abs(Mtotz[0, 0]) ** 2
                + np.abs(q_sub * Mtotz[0, 1] / 1j) ** 2
            ))

    return z_prof, E_prof / Epeak_0


# ── Absorption ───────────────────────────────────────────────────────────

def calc_abs(Esq, Ls, alphas, n_pts=None):
    """Integrated absorption from an E-field profile.

    Parameters
    ----------
    Esq : ndarray
        |E|^2 inside each layer (from :func:`field_zmag`),
        ordered as ``n_layers * n_pts`` samples.
    Ls : array_like
        Physical thicknesses [m].
    alphas : array_like
        Absorption coefficients [1/m].
    n_pts : int, optional
        Samples per layer.  Default: ``len(Esq) // len(Ls)``.

    Returns
    -------
    float
        Integrated stack absorption [W/W].
    """
    if n_pts is None:
        n_pts = len(Esq) // len(Ls)
    absorp = 0.0
    for i, (alpha_i, Li) in enumerate(zip(alphas, Ls)):
        E_layer = Esq[i * n_pts : (i + 1) * n_pts]
        absorp += 2 * _trapz(E_layer * alpha_i, dx=Li / n_pts)
    return absorp


# ── Utilities ────────────────────────────────────────────────────────────

def op2phys(L, n):
    """Convert optical thicknesses to physical thicknesses.

    Parameters
    ----------
    L : array_like
        Optical thicknesses (fraction of design wavelength).
    n : array_like
        Refractive indices (same length as *L*).

    Returns
    -------
    ndarray
        Physical thicknesses (fraction of design wavelength / n).

    Raises
    ------
    ValueError
        If lengths don't match.
    """
    if len(L) != len(n):
        raise ValueError(
            f'L (dim {len(L)}) and n (dim {len(n)}) must have the same dimension.'
        )
    return np.asarray(L) / np.asarray(n)


def sellmeier(B=None, C=None, lam=1064e-9):
    """Refractive index from Sellmeier dispersion coefficients.

    Parameters
    ----------
    B : list or array_like, optional
        Sellmeier B coefficients. Default is fused silica.
    C : list or array_like, optional
        Sellmeier C coefficients in um^2. Default is fused silica.
    lam : float or array_like, optional
        Wavelength(s) in meters. Default 1064 nm.

    Returns
    -------
    float or ndarray
        Refractive index.
    """
    if B is None:
        B = [0.696166300, 0.407942600, 0.897479400]
    if C is None:
        C = [4.67914826e-3, 1.35120631e-2, 97.9340025]
    if len(B) != len(C):
        raise ValueError(
            f'B (len {len(B)}) and C (len {len(C)}) must have the same length.'
        )
    ll = np.asarray(lam) * 1e6  # wavelength in microns
    n_sq = 1.0
    for b, c in zip(B, C):
        n_sq = n_sq + b * ll ** 2 / (ll ** 2 - c)
    return np.sqrt(n_sq)


def qwbandedges(stack):
    """Quarter-wave HR band edges (Orfanidi Eq. 6.3.18).

    Parameters
    ----------
    stack : dict
        Stack dict with ``ns`` and ``Ls``.

    Returns
    -------
    lam1 : float
        Short-wavelength band edge [m].
    lam2 : float
        Long-wavelength band edge [m].
    """
    ns = np.real(stack['ns'][1:-1])
    Ls = stack['Ls']
    nH = np.max(ns)
    nL = np.min(ns)
    LH = Ls[np.argmax(ns)]
    LL = Ls[np.argmin(ns)]
    rho = (nH - nL) / (nH + nL)
    lam1 = np.pi * (nH * LH + nL * LL) / np.arccos(-rho)
    lam2 = np.pi * (nH * LH + nL * LL) / np.arccos(rho)
    return lam1, lam2
