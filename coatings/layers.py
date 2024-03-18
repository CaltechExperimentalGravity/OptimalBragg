import numpy as np
import scipy.io as scio

from scipy.interpolate import interp1d, PchipInterpolator


def multilayer_diel(ns, Ls, lamb, lamb_0, aoi=0, pol="te"):
    """Calculates amplitude reflectivity and complex
    impedance of a dielectric stack according to Eqs. (6.1.3)
    and (7.7.1) of http://eceweb1.rutgers.edu/~orfanidi/ewa/

    Args:
        ns (arr): Array of refractive indices, including the
                 incident and transmitted media. Ordered from
                 incident to transmitted medium.
        Ls (arr): Array of physical thicknesses comprising
                 the dielectric stack, ordered from incident
                 to transmitted medium. Should have 2 fewer
                 elements than n.
        lamb (float, arr): Wavelengths for which to evaluate
                           the stack reflectivity.
        lamb_0 (float): Center design wavelength for which to
                        evaluate the stack reflectivity.
        aoi (float, optional): Angle of incidence in rad,
                               default=0.0 (normal inc)
        pol (str, optional): Polarization at which to evaluate
                             reflectivity, defaults to 'te' or
                             s-pol.

    Returns:
        Gamma_0 (arr): Amplitude reflectivity (complex)
        z_0 (arr): Complex impedance at the interface, only for
                   incident medium with n=1.0, i.e. vacuum.

    """
    # Use center design wavelength as reference unit
    rel_lamb = lamb / lamb_0

    # Cosine of theta; projection for oblique incidence
    proj_aoi = np.conj(1 - (ns[0] * np.sin(aoi) / ns) ** 2)

    # Final conj needed when n_inc > n(i) and aoi > aoi_crit
    costh = np.conj(np.sqrt(proj_aoi))

    # Polarization projection in transmission
    if pol in ["te", "s", "TE", "S"]:
        nT = ns * costh
    else:
        nT = ns / costh

    # Optical thicknesses
    opt_L = (ns[1:-1] * Ls / lamb) * costh[1 : len(Ls) + 1]

    # Per-layer amplitude reflectivity
    r = -np.diff(nT) / (np.diff(nT) + 2 * nT[0 : len(Ls) + 1])

    # Initialize to incident layer amplitude reflectivity
    Gamma_0 = r[len(Ls)] * np.ones_like(rel_lamb)

    # Recursion relation
    for i in range(len(Ls) - 1, -1, -1):
        delta = 2 * np.pi * opt_L[i] / rel_lamb
        z = np.exp(-2 * 1j * delta)
        Gamma_0 = (r[i] + Gamma_0 * z) / (1 + r[i] * Gamma_0 * z)

    # Incident layer impedance
    Z_0 = (1 + Gamma_0) / (1 - Gamma_0)

    return Gamma_0, Z_0


def surfield(rr, Ei=27.46, normalized=False):
    """Surface electric field for a dielectric coating

    Args:
        rr (complex): Amplitude reflectivity at the interface
        Ei (float): Incident electric field amplitude in V/m.
                    default = 27.46 V/m corresponding to 1 W/m^2.
        normalized (bool, optional): Return in units of Ei

    Returns:
        (float): Surface E-field amplitude
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
        Ls (arr): Physical thickness
        lam (float): Wavelength
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
    Parameters:
    -----------
    Esq: array_like
        NORMALIZED electric field squared as a function of distance in a coating
    Ls: array_like
        PHYSICAL thickness of coating layers
    alphas: arr
        Absorption coefficients of all layers in 1/m
    Returns:
    -------------
    absorp: float
        Absorption of coating (W/W)
    """

    # Remember, absorption coefficient = 4 * pi * k / lambda,
    # where k is the extinction coefficient (aka Im(n))
    # Should the transfer matrix method include this?

    absorp = 0
    for alpha_i, Li in zip(alphas, Ls):
        # Unclear if trapz is faster, probably not... but also not slower.
        # absorp += np.sum(Esq * alpha_i * np.ones_like(Esq) * Li / len(Esq))
        absorp += np.trapz(
            Esq * alpha_i * np.ones_like(Esq), dx=Li / (len(Esq))
        )
    return absorp


def stack_refl(wavelengths, stack):
    lam_ref = stack["lam_ref"]
    ns, Ls = stack["ns"], stack["Ls"]
    try:
        iter(wavelengths)
        refl = np.zeros_like(wavelengths)
        for ii, wavelength in enumerate(wavelengths):
            rr, _ = multilayer_diel(ns, Ls, wavelength, lam_ref)
            refl[ii] = rr
    except TypeError:
        rr, _ = multilayer_diel(ns, Ls, wavelengths, lam_ref)
        refl = rr
    finally:
        return refl


def stack_R(wavelengths, stack):
    lam_ref = stack["lam_ref"]
    ns, Ls = stack["ns"], stack["Ls"]
    try:
        iter(wavelengths)
        RR = np.zeros_like(wavelengths)
        for ii, wavelength in enumerate(wavelengths):
            rr, _ = multilayer_diel(ns, Ls, wavelength, lam_ref)
            RR[ii] = np.abs(rr) ** 2
    except TypeError:
        rr, _ = multilayer_diel(ns, Ls, wavelengths, lam_ref)
        RR = np.abs(rr) ** 2
    finally:
        return RR
