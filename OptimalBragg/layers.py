import numpy as np
import scipy.io as scio

from scipy.interpolate import interp1d, PchipInterpolator


def multilayer_diel(ns, Ls, lamb, aoi=0, pol="te"):
    """Calculates amplitude reflectivity and complex
    impedance of a dielectric stack according to Chapter 8
    of http://eceweb1.rutgers.edu/~orfanidi/ewa/

    Args:
        ns (arr): Array of refractive indices, including the
                 incident and transmitted media. Ordered from
                 incident to transmitted medium.
        Ls (arr): Array of physical thicknesses comprising
                 the dielectric stack, ordered from incident
                 to transmitted medium. Should have 2 fewer
                 elements than n.
        lamb (float, arr): Wavelength(s) for which to evaluate
                           the stack reflectivity.
        aoi (float, optional): Angle of incidence in rad,
                               default=0.0 (normal inc)
        pol (str, optional): Polarization at which to evaluate
                             reflectivity, defaults to 'te' or
                             's' polarization.

    Returns:
        Gamma_0 (arr): Amplitude reflectivity (complex)
        z_0 (arr): Complex impedance at the interface, only for
                   incident medium with n=1.0, i.e. vacuum.

    """
    # Wavenumber(s)
    ki = 2 * np.pi / lamb

    # Oblique incidence projection (Eq. 8.2.2)
    proj_aoi = np.conj(1 - (ns[0] * np.sin(aoi) / ns) ** 2)

    # Final conj needed when n_inc > n(i) and aoi > aoi_crit
    costh = np.conj(np.sqrt(proj_aoi))

    # Transverse index projection (Eq. 8.1.4)
    if pol in ["te", "s", "TE", "S"]:
        nT = ns * costh
    else:
        nT = ns / costh

    # Oblique, lossy optical thicknesses
    opt_L = (ns[1:-1] * Ls) * costh[1 : len(Ls) + 1]

    # Per-layer amplitude reflectivity (Eq. 8.1.3)
    r = -np.diff(nT) / (np.diff(nT) + 2 * nT[0 : len(Ls) + 1])

    # Initialize complex amplitude reflectivities with incident layer
    Gamma_0 = r[len(Ls)] * np.ones_like(lamb)

    # Recursion relation (Eq. 8.1.7)
    for i in range(len(Ls) - 1, -1, -1):
        delta_i = ki * opt_L[i]
        z_i = np.exp(2j * delta_i)
        Gamma_0 = (r[i] + Gamma_0 * z_i) / (1 + r[i] * Gamma_0 * z_i)

    # Incident layer complex impedance
    Z_0 = (1 + Gamma_0) / (1 - Gamma_0)
    return Gamma_0, Z_0


def surfield(rr, Ei=27.46, normalized=False):
    """Surface electric field for a dielectric coating

    Args:
        rr (complex): Amplitude reflectivity at the input interface
        Ei (float): Incident electric field amplitude in V/m.
                    default = 27.46 V/m corresponding to 1 W/m^2.
        normalized (bool, optional): Return in units of Ei

    Returns:
        (float, arr): Surface E-field amplitude
    """
    if normalized:
        return np.abs(1 + rr)
    else:
        return Ei * np.abs(1 + rr)


def field_zmag(ns, Ls, lam, aoi=0, pol="s", n_pts=30):
    """Normalized longitudinal E-field strength squared following
    derivation set out in Arnon and Baumeister, 1980
    https://www.osapublishing.org/ao/abstract.cfm?uri=ao-19-11-1853

    Args:
        ns (arr): Refractive index
        Ls (arr): Physical thickness [m]
        lam (float): Reference wavelength [m]
        aoi (float, optional): Angle of incidence (rad),
                               default = 0
        pol (str, optional): Polarization, default = 's'.
        n_pts (int, optional): Number of points

    Returns:
        z_prof (arr): Array of penetration depths at which
                     field magnitude is evaluated.
        (arr): Magnitude squared electric field normalized
               to value at the interface of incidence.
    """

    # Calculate array of incidence angles
    alpha = [aoi]
    for ii in range(len(ns) - 1):
        t_r = np.arcsin(ns[ii] * np.sin(alpha[ii]) / ns[ii + 1])
        alpha.append(t_r)
    q_angle = alpha[-1]
    angles = np.array(alpha[1:-1])

    def M_i(b_i, qq_i):
        # Eq (2) from Arnon and Baumeister, 1980
        out = np.matrix(
            [
                [np.cos(b_i), 1j * np.sin(b_i) / qq_i],
                [1j * np.sin(b_i) * qq_i, np.cos(b_i)],
            ]
        )
        return out

    def q_i(n_i, theta_i):
        # Eqs (5)-(6) from Arnon and Baumeister, 1980
        if pol in ["te", "TE", "s", "S"]:
            return n_i * np.cos(theta_i)
        elif pol in ["tm", "TM", "p", "P"]:
            return n_i / np.cos(theta_i)

    def beta_i(tt_i, nn_i, hh_i):
        # Eq (3) from Arnon and Baumeister, 1980
        return 2 * np.pi * np.cos(tt_i) * nn_i * hh_i / lam

    # Calculate the total matrix as per Eq (7)
    Mtot = np.eye(2)
    for n_i, L_i, a_i in zip(ns[1:-1], Ls, angles):
        Mtot = Mtot * M_i(beta_i(a_i, n_i, L_i), q_i(n_i, a_i))

    # Eq (10) from Arnon and Baumeister, 1980
    q_0 = q_i(ns[0], aoi)
    q_sub = q_i(ns[-1], q_angle)
    Epeak_0 = 0.25 * (
        np.abs(Mtot[0, 0] + Mtot[1, 1] * q_sub / q_0) ** 2
        + np.abs(Mtot[1, 0] / q_0 / 1j + Mtot[0, 1] * q_sub / 1j) ** 2
    )

    def delta_h(bb_i, qq_i):
        # Eq (11) from Arnon and Baumeister, 1980
        return M_i(bb_i, -qq_i)

    # Initialize some arrays to store the calculated E field profile
    E_prof = np.zeros(len(Ls) * n_pts)
    z_prof = np.zeros(len(Ls) * n_pts)
    Z = 0
    Mtotz = Mtot

    # Initialize the q-parameter at the rightmost interface
    q_sub = q_i(ns[-1], q_angle)

    for ii in range(len(Ls)):
        n_i = ns[ii + 1]
        dL = Ls[ii] / n_pts
        a_i = angles[ii]

        if pol in ["tm", "TM", "p", "P"]:
            corr = (np.cos(aoi) / np.cos(a_i)) ** 2
        else:
            corr = 1

        for jj in range(0, n_pts):
            Z += dL
            z_prof[ii * n_pts + jj] = Z
            Mtotz = delta_h(beta_i(a_i, n_i, dL), q_i(n_i, a_i)) * Mtotz
            E_prof[ii * n_pts + jj] = corr * (
                np.abs(Mtotz[0, 0]) ** 2 + np.abs(q_sub * Mtotz[0, 1] / 1j) ** 2
            )

    return z_prof, E_prof / Epeak_0


def calc_abs(Esq, Ls, alphas):
    """
    Function for calculating the Absorption given an electric field profile

    Args:
        Esq (arr): Magnitude squared E-field inside layers
        Ls (arr): Physical thicknesses [m]
        alphas (arr): Absorption coefficients [m]

    Returns:
        (float): Integrated stack absorption [W/W]
    """

    # Remember, absorption coefficient = 4 * pi * k / lambda,
    # where k is the extinction coefficient (aka Im(n))
    # We should let the transfer matrix method handle this for speedup
    # Note the factor of 2, accounting for the two field passes in the stack
    absorp = 0
    for alpha_i, Li in zip(alphas, Ls):
        # Unclear if trapz is faster, probably not... but also not slower.
        # absorp += np.sum(Esq * alpha_i * np.ones_like(Esq) * Li / len(Esq))
        absorp += 2 * np.trapz(
            Esq * alpha_i * np.ones_like(Esq), dx=Li / (len(Esq))
        )
    return absorp


def amp_refl(wavelengths, stack, **multilayer_diel_pars):
    """Compute spectral amplitude reflectivity (r) of a given stack

    Args:
        wavelengths (arr): Wavelenght(s) [m]
        stack (dict): Container of stack attributes
        **multilayer_diel_pars: Kwargs for multilayer_diel

    Returns:
        (arr): Spectral amplitude reflectivities
    """
    ns, Ls = stack["ns"], stack["Ls"]
    try:
        iter(wavelengths)
        rrj = np.zeros_like(wavelengths).astype(complex)
        for j, wavelength in enumerate(wavelengths):
            rr, _ = multilayer_diel(ns, Ls, wavelength, **multilayer_diel_pars)
            rrj[j] = rr
    except TypeError:
        rrj, _ = multilayer_diel(ns, Ls, wavelengths, **multilayer_diel_pars)
    finally:
        return rrj


def refl(wavelengths, stack, **multilayer_diel_pars):
    """Compute spectral reflection (R) of a given stack

    Args:
        wavelengths (arr): Wavelength(s) [m]
        stack (dict): Container of stack attributes
        **multilayer_diel_pars: Kwargs for multilayer_diel

    Returns:
        (arr): Spectral reflection
    """
    return np.abs(amp_refl(wavelengths, stack, **multilayer_diel_pars)) ** 2


def trans(wavelengths, stack, **multilayer_diel_pars):
    """Compute spectral transmission (T) of a given stack

    Args:
        wavelengths (arr): Wavelength(s) [m]
        stack (dict): Container of stack attributes
        **multilayer_diel_pars: Kwargs for multilayer_diel

    Returns:
        (arr): Spectral transmission
    """
    return 1 - refl(wavelengths, stack, **multilayer_diel_pars)


def qwbandedges(stack):
    """Based on Eq (6.3.18) in orfanidi

    Args:
        stack (dict): Container of stack attributes

    Returns:
        (float, float): Wavelengths at which HR band edges occur [m]
    """
    ns = stack["ns"][1:-1]
    nH = max(ns)
    nL = min(ns)
    LH = stack["Ls"][np.argmax(ns)]
    LL = stack["Ls"][np.argmin(ns)]
    rho = (nH - nL) / (nH + nL)
    lam1 = np.pi * (nH * LH + nL * LL) / (np.arccos(-rho))
    lam2 = np.pi * (nH * LH + nL * LL) / (np.arccos(rho))
    return lam1, lam2
