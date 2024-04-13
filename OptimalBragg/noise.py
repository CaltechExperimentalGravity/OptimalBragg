# Adapted from pygwinc!
import numpy as np
from scipy.constants import k, c, hbar
from scipy.special import exp1, jn, jn_zeros

zeta = jn_zeros(1, 300)
j0m = jn(0, zeta)


def substrate_noise(freq, stack, w_beam, r_mirror=None, d_mirror=None):
    sub_Sbr = substrate_brownian(
        f=freq, stack=stack, w_beam=w_beam, r_mirror=r_mirror, d_mirror=d_mirror
    )
    sub_Str = substrate_thermorefractive(
        f=freq, stack=stack, w_beam=w_beam, d_mirror=d_mirror
    )
    sub_Ste = substrate_thermoelastic(
        f=freq, stack=stack, w_beam=w_beam, r_mirror=r_mirror, d_mirror=d_mirror
    )
    return sub_Sbr, sub_Str, sub_Ste


def coating_noise(
    freq, stack, w_beam, power=None, r_mirror=None, d_mirror=None, m_mirror=None
):
    coat_Sbr = coating_brownian(
        f=freq, stack=stack, w_beam=w_beam, power=power, mass=m_mirror
    )

    coat_Sto, coat_Ste, coat_Str = coating_thermooptic(
        f=freq, stack=stack, w_beam=w_beam, r_mirror=r_mirror, d_mirror=d_mirror
    )

    return coat_Sbr, coat_Ste, coat_Str, coat_Sto


def coating_thermooptic(
    f, stack, w_beam, r_mirror=None, d_mirror=None, **kwargs
):
    """Optical coating thermo-optic displacement noise spectrum

    :f: frequency array in Hz
    :stack: stack dict
    :wBeam: beam radius (at 1 / e**2 power)

    :returns: tuple of:
    SteZ = thermo-elastic component of thermooptic noise
    StrZ = thermo-refractive component of thermooptic noise
    """
    # compute coefficients
    dTO, dTR, dTE, T, junk = getCoatTOPos(stack, w_beam, r_mirror, d_mirror)

    # compute correction factors
    gTO = getCoatThickCorr(f, stack, dTE, dTR)
    gTE = getCoatThickCorr(f, stack, dTE, 0)
    gTR = getCoatThickCorr(f, stack, 0, dTR)

    # compute thermal source spectrum
    SsurfT, junk = getCoatThermal(f, stack, w_beam)

    StoZ = SsurfT * gTO * dTO**2
    SteZ = SsurfT * gTE * dTE**2
    StrZ = SsurfT * gTR * dTR**2

    return StoZ, SteZ, StrZ


def getCoatThickCorr(f, stack, dTE, dTR):
    """Finite coating thickness correction

    :f: frequency array in Hz
    :stack: stack dict
    :wBeam: beam radius (at 1 / e**2 power)

    Uses correction factor from LIGO-T080101, "Thick Coating
    Correction" (Evans).

    See getCoatThermoOptic for example usage.

    """
    ##############################################
    # For comparison in the bTR = 0 limit, the
    # equation from Fejer (PRD70, 2004)
    # modified so that gFC -> 1 as xi -> 0
    #  gTC = (2 ./ (R * xi.**2)) .* (sh - s + R .* (ch - c)) ./ ...
    #    (ch + c + 2 * R * sh + R**2 * (ch - c));
    # which in the limit of xi << 1 becomes
    #  gTC = 1 - xi * (R - 1 / (3 * R));

    # parameter extraction
    pS = stack["sub"]
    Cs = pS.CP * pS.rho
    Ks = pS.Kappa

    # compute coating average parameters
    dc, Cc, Kc, junk = getCoatAvg(stack)

    # R and xi (from T080101, Thick Coating Correction)
    w = 2 * np.pi * f
    R = np.sqrt(Cc * Kc / (Cs * Ks))
    xi = dc * np.sqrt(2 * w * Cc / Kc)

    # trig functions of xi
    s = np.sin(xi)
    c = np.cos(xi)
    sh = np.sinh(xi)
    ch = np.cosh(xi)

    # pR and pE (dTR = -\bar{\beta} lambda, dTE = \Delta \bar{\alpha} d)
    pR = dTR / (dTR + dTE)
    pE = dTE / (dTR + dTE)

    # various parts of gTC
    g0 = 2 * (sh - s) + 2 * R * (ch - c)
    g1 = 8 * np.sin(xi / 2) * (R * np.cosh(xi / 2) + np.sinh(xi / 2))
    g2 = (1 + R**2) * sh + (1 - R**2) * s + 2 * R * ch
    gD = (1 + R**2) * ch + (1 - R**2) * c + 2 * R * sh

    # and finally, the correction factor
    gTC = (pE**2 * g0 + pE * pR * xi * g1 + pR**2 * xi**2 * g2) / (
        R * xi**2 * gD
    )

    return gTC


def getCoatThermal(f, stack, w_beam):
    """Thermal noise spectra for a surface layer

    :f: frequency array in Hz
    :stack: stack dict
    :wBeam: beam radius (at 1 / e**2 power)

    :returns: tuple of:
    SsurfT = power spectra of thermal fluctuations in K**2 / Hz
    rdel = thermal diffusion length at each frequency in m

    """
    pS = stack["sub"]
    C_S = pS.CP * pS.rho
    K_S = pS.Kappa
    kBT2 = k * pS.Temp**2

    # omega
    w = 2 * np.pi * f

    # thermal diffusion length
    rdel = np.sqrt(2 * K_S / (C_S * w))

    # noise equation
    SsurfT = 4 * kBT2 / (np.pi * w * C_S * rdel * w_beam**2)

    return SsurfT, rdel


def getCoatTOPos(stack, w_beam, r_mirror=None, d_mirror=None):
    """Mirror position derivative wrt thermal fluctuations

    :mirror: mirror parameter Struct
    :wavelength: laser wavelength
    :wBeam: beam radius (at 1 / e**2 power)

    :returns: tuple of:
    dTO = total thermo-optic dz/dT
    dTR = thermo-refractive dz/dT
    dTE = thermo-elastic dz/dT
    T = coating power transmission
    R = coating power reflection

    Compute thermal fluctuations with getCoatThermal.

    See LIGO-T080101.

    """
    # parameters
    wavelength = stack["lam_ref"]
    sub = stack["sub"]
    coat = stack["thin_films"]
    ns, Ls = stack["ns"], stack["Ls"]
    nS = sub.Index
    dOpt = stack["Ls"] * stack["ns"][1:-1] / wavelength

    # compute refractive index, effective alpha and beta
    nLayer, aLayer, bLayer, dLayer, sLayer = getCoatLayers(stack)

    # compute coating average parameters
    dc, Cc, Kc, aSub = getCoatAvg(stack)

    # compute reflectivity and parameters
    dphi_dT, dphi_TE, dphi_TR, rCoat = getCoatTOPhase(
        1, nS, nLayer, dOpt, aLayer, bLayer, sLayer
    )
    R = abs(rCoat) ** 2
    T = 1 - R

    # for debugging
    # print(T / 1e-6)

    # convert from phase to meters, subtracting substrate
    dTR = dphi_TR * wavelength / (4 * np.pi)
    dTE = dphi_TE * wavelength / (4 * np.pi) - aSub * dc

    # mirror finite size correction
    if r_mirror is not None and d_mirror is not None:
        Cfsm = getCoatFiniteCorr(stack, w_beam, r_mirror, d_mirror)
    else:
        Cfsm = 1

    dTE = dTE * Cfsm

    # add TE and TR effects (sign is already included)
    dTO = dTE + dTR

    return dTO, dTR, dTE, T, R


def getCoatFiniteCorr(stack, w_beam, r_mirror, d_mirror):
    """Finite mirror size correction

    :stack: mirror parameter dict
    :wBeam: beam radius (at 1 / e**2 power)

    Uses correction factor from PLA 2003 vol 312 pg 244-255
    "Thermodynamical fluctuations in optical mirror coatings"
    by V. B. Braginsky and S. P. Vyatchanin
    http://arxiv.org/abs/cond-mat/0302617

    (see getCoatTOPos for example usage)

    version 1 by Sam Wald, 2008

    """
    # parameter extraction
    wavelength = stack["lam_ref"]
    nC = stack["ns"][1:-1]
    dC = stack["Ls"]
    dOpt = dC * nC / wavelength

    alphaS = stack["sub"].Alpha
    C_S = stack["sub"].CP * stack["sub"].rho
    Y_S = stack["sub"].Y
    sigS = stack["sub"].Sigma

    alphaC = stack["ctes"]
    Cv_C = stack["Cvs"]
    Y_C = stack["Ys"]
    sigmaC = stack["sigmas"]

    # coating sums
    dcsum = np.sum(dC)

    # AVERAGE SPECIFIC HEAT (simple volume average for coating)

    Cf = np.sum(Cv_C * dC) / dcsum
    Cr = Cf / C_S

    # COATING AVERAGE VALUE X = ALPHAF*(1+POISSONf)/(1-POISSONf) avg
    xxC = alphaC * (1 + sigmaC) / (1 - sigmaC)
    Xf = np.sum(xxC * dC) / dcsum
    Xr = Xf / alphaS

    # COATING AVERAGE VALUE Y = ALPHAF* YOUNGSF/(1-POISSONF) avg
    yyC = alphaC * Y_C / (1 - sigmaC)
    Yf = np.sum(yyC * dC) / dcsum
    Yr = Yf / (alphaS * Y_S)

    # COATING AVERAGE VALUE Z = 1/(1-POISSONF) avg
    zzC = 1 / (1 - sigmaC)
    Zf = np.sum(zzC * dC) / dcsum

    #################################### FINITE SIZE CORRECTION CALCULATION

    # beam size parameter used by Braginsky
    r0 = w_beam / np.sqrt(2)

    # between eq 77 and 78
    km = zeta / r_mirror
    Qm = np.exp(-2 * km * d_mirror)
    pm = (
        np.exp(-(km**2) * r0**2 / 4) / j0m
    )  # left out factor of pi * R**2 in denominator

    # eq 88
    Lm = (
        Xr
        - Zf * (1 + sigS)
        + (Yr * (1 - 2 * sigS) + Zf - 2 * Cr)
        * (1 + sigS)
        * (1 - Qm) ** 2
        / ((1 - Qm) ** 2 - 4 * km**2 * d_mirror**2 * Qm)
    )

    # eq 90 and 91
    S1 = (12 * r_mirror**2 / d_mirror**2) * np.sum(pm / zeta**2)
    S2 = np.sum(pm**2 * Lm**2)
    P = (Xr - 2 * sigS * Yr - Cr + S1 * (Cr - Yr * (1 - sigS))) ** 2 + S2

    # eq 60 and 70
    LAMBDA = -Cr + (Xr / (1 + sigS) + Yr * (1 - 2 * sigS)) / 2

    # eq 92
    Cfsm = np.sqrt(
        (r0**2 * P) / (2 * r_mirror**2 * (1 + sigS) ** 2 * LAMBDA**2)
    )

    return Cfsm


def getCoatLayers(stack):
    """Layer vectors for refractive index, effective alpha and beta and geometrical thickness

    :mirror: mirror parameter Struct
    :wavelength: laser wavelength

    :returns: tuple of:
    nLayer = refractive index of each layer, ordered input to output (N x 1)
    aLayer = change in geometrical thickness with temperature
           = the effective thermal expansion coeffient of the coating layer
    bLayer = change in refractive index with temperature
           = dn/dT
    dLayer = geometrical thicness of each layer
    sLayer = Yamamoto thermo-refractive correction
           = alpha * (1 + sigma) / (1 - sigma)

    """
    # coating parameters
    wavelength = stack["lam_ref"]
    nC = stack["ns"][1:-1]
    dOpt = stack["Ls"] * nC / wavelength

    Y_S = stack["sub"].Y
    sigS = stack["sub"].Sigma

    Y_C = stack["Ys"]
    cteC = stack["ctes"]
    sigmaC = stack["sigmas"]
    betaC = stack["betas"]

    Nlayer = len(dOpt)

    # compute effective alpha
    def getExpansionRatio(Y_C, sigC, Y_S, sigS):
        ##############################################
        # Y_C and sigC are for the coating material (can also be substrate)
        # Y_S and sigS are for the substrate material
        #
        ce = ((1 + sigS) / (1 - sigC)) * (
            ((1 + sigC) / (1 + sigS)) + (1 - 2 * sigS) * Y_C / Y_S
        )
        return ce

    aLayer = cteC * getExpansionRatio(Y_C, sigmaC, Y_S, sigS)
    bLayer = betaC
    dLayer = wavelength * np.asarray(dOpt) / nC
    sLayer = cteC * (1 + sigmaC) / (1 - sigmaC)

    return nC, aLayer, bLayer, dLayer, sLayer


def getCoatAvg(stack):
    """Coating average properties

    :stack: dict stack

    :returns: tuple of:
    dc = total thickness (meters)
    Cc = heat capacity
    Kc = thermal diffusivity
    aSub = effective substrate thermal expansion (weighted by heat capacity)

    """
    # coating parameters
    wavelength = stack["lam_ref"]
    nC = stack["ns"][1:-1]
    dOpt = stack["Ls"] * nC / wavelength

    alphaS = stack["sub"].Alpha
    C_S = stack["sub"].CP * stack["sub"].rho
    sigS = stack["sub"].Sigma

    Cv_C = stack["Cvs"]
    Kd_C = stack["thermaldiffs"]
    dLayer = stack["Ls"]

    # heat capacity and thermal diffusivity
    pattern = stack["pattern"]
    dcsum = np.sum(dLayer)
    Cc = np.sum(Cv_C * dLayer) / dcsum
    Kc = dcsum / np.sum(dLayer / Kd_C)

    # effective substrate thermal expansion
    aSub = 2 * alphaS * (1 + sigS) * Cc / C_S

    return dcsum, Cc, Kc, aSub


def getCoatTOPhase(nIn, nOut, nLayer, dOpt, aLayer, bLayer, sLayer):
    """Coating reflection phase derivatives w.r.t. temperature

    :nIn: refractive index of input medium (e.g., vacuum = 1)
    :nOut: refractive index of output medium (e.g., SiO2 = 1.45231 @ 1064nm)
    :nLayer: refractive index of each layer, ordered input to output (N x 1)
    :dOpt: optical thickness / lambda of each layer
           = geometrical thickness * refractive index / lambda
    :aLayer: change in geometrical thickness with temperature
             = the effective thermal expansion coeffient of the coating layer
    :bLayer: change in refractive index with temperature
             = dn/dT
             = dd/dT - n * a

    :returns: tuple of:
    dphi_dT = total thermo-optic phase derivative with respect to temperature
            = dphi_TE + dphi_TR
    dphi_TE = thermo-elastic phase derivative (dphi / dT)
    dphi_TR = thermo-refractive phase derivative (dphi / dT)
    rCoat = amplitude reflectivity of coating (complex)

    Note about aLayer: on a SiO2 substrate,
    a_Ta2O5 ~ 3.5 * alpha_Ta2O5
    a_SiO2 ~ 2.3 * alpha_SiO2

    See :getCoatTOPos: for more information.

    See LIGO-T080101.

    """
    # coating reflectivity calc
    rCoat, dcdp = getCoatRefl2(nIn, nOut, nLayer, dOpt)[:2]

    # geometrical distances
    dGeo = np.asarray(dOpt) / nLayer

    # phase derivatives
    dphi_dd = 4 * np.pi * dcdp

    # thermo-refractive coupling
    dphi_TR = np.sum(dphi_dd * (bLayer + sLayer * nLayer) * dGeo)

    # thermo-elastic
    dphi_TE = 4 * np.pi * np.sum(aLayer * dGeo)

    # total
    dphi_dT = dphi_TR + dphi_TE

    return dphi_dT, dphi_TE, dphi_TR, rCoat


def substrate_thermoelastic(f, stack, w_beam, r_mirror=None, d_mirror=None):
    """Substrate thermal displacement noise spectrum from thermoelastic fluctuations

    :f: frequency array in Hz
    :stack: stack dict
    :wBeam: beam radius (at 1 / e^2 power)

    :returns: displacement noise power spectrum at :f:, in meters

    """
    sub = stack["sub"]
    sigma = sub.Sigma
    rho = sub.rho
    kappa = sub.Kappa  # thermal conductivity
    alpha = sub.Alpha  # thermal expansion
    CM = sub.CP  # heat capacity @ constant mass
    Temp = sub.Temp  # temperature
    kBT = k * Temp

    S = (
        8 * (1 + sigma) ** 2 * kappa * alpha**2 * Temp * kBT
    )  # note kBT has factor Temp
    S /= np.sqrt(2 * np.pi) * (CM * rho) ** 2
    S /= (w_beam / np.sqrt(2)) ** 3  # LT 18 less factor 1/omega^2

    # Corrections for finite test masses:
    if r_mirror is not None and d_mirror is not None:
        S *= substrate_thermoelastic_FiniteCorr(
            stack, w_beam, r_mirror, d_mirror
        )

    return S / (2 * np.pi * f) ** 2


def substrate_thermoelastic_FiniteCorr(stack, w_beam, r_mirror, d_mirror):
    """Substrate thermoelastic noise finite-size test mass correction

    :materials: gwinc optic materials structure
    :wBeam: beam radius (at 1 / e^2 power)

    :returns: correction factor

    (Liu & Thorne gr-qc/0002055 equation 46)

    Equation references to Bondu, et al. Physics Letters A 246 (1998)
    227-236 (hereafter BHV) or Liu and Thorne gr-qc/0002055 (hereafter LT)

    """
    sub = stack["sub"]
    a = r_mirror
    h = d_mirror
    sigma = sub.Sigma

    # LT uses power e-folding
    r0 = w_beam / np.sqrt(2)
    km = zeta / a

    Qm = np.exp(-2 * km * h)  # LT 35a

    pm = np.exp(-((km * r0) ** 2) / 4) / (np.pi * (a * j0m) ** 2)  # LT 37

    c0 = 6 * (a / h) ** 2 * sum(j0m * pm / zeta**2)  # LT 32
    c1 = -2 * c0 / h  # LT 32
    p0 = 1 / (np.pi * a**2)  # LT 28
    c1 += p0 / (2 * h)  # LT 40

    coeff = (1 - Qm) * ((1 - Qm) * (1 + Qm) + 8 * h * km * Qm)
    coeff += 4 * (h * km) ** 2 * Qm * (1 + Qm)
    coeff *= km * (pm * j0m) ** 2 * (1 - Qm)
    coeff /= ((1 - Qm) ** 2 - 4 * (h * km) ** 2 * Qm) ** 2
    coeff = sum(coeff) + h * c1**2 / (1 + sigma) ** 2
    coeff *= (np.sqrt(2 * np.pi) * r0) ** 3 * a**2  # LT 46

    return coeff


def substrate_brownian(f, stack, w_beam, r_mirror=None, d_mirror=None):
    """Thermal displacement noise ASD due to substrate mechanical loss

    :f: frequency array in Hz
    :stack: stack dict
    :wBeam: beam radius (at 1 / e^2 power)

    :returns: displacement noise power spectrum at :f:, in meters

    Args:
        f (array): Fourier frequency in Hz
        stack (dict): Description
        w_beam (TYPE): Description
        r_mirror (None, optional): Description
        d_mirror (None, optional): Description

    Returns:
        TYPE: Description

    """
    wavelength = stack["lam_ref"]
    sub = stack["sub"]
    coat = stack["thin_films"]
    Y = sub.Y
    sigma = sub.Sigma
    c2 = sub.c2
    n = sub.MechanicalLossExponent
    alphas = sub.Alpha
    kBT = k * sub.Temp

    if r_mirror is not None and d_mirror is not None:
        cftm, aftm = substrate_brownian_FiniteCorr(
            stack, w_beam, r_mirror, d_mirror
        )
    else:
        cftm, aftm = 0, 0

    # Bulk substrate contribution
    phibulk = c2 * f**n
    cbulk = 8 * kBT * aftm * phibulk / (2 * np.pi * f)

    # Surface loss contribution
    # csurf = alphas/(Y*pi*wBeam^2)
    csurf = alphas * (1 - 2 * sigma) / ((1 - sigma) * Y * np.pi * w_beam**2)
    csurf *= 8 * kBT / (2 * np.pi * f)

    return csurf + cbulk


def substrate_thermorefractive(f, stack, w_beam, d_mirror, exact=False):
    """Substrate thermal displacement noise spectrum from thermorefractive fluctuations

    :f: frequency array in Hz
    :materials: gwinc optic materials structure
    :wBeam: beam radius (at 1 / e^2 power)
    :exact: whether to use adiabatic approximation or exact calculation (False)

    :returns: displacement noise power spectrum at :f:, in meters

    """
    wavelength = stack["lam_ref"]
    sub = stack["sub"]
    coat = stack["thin_films"]
    H = d_mirror
    kBT = k * sub.Temp
    Temp = sub.Temp
    rho = sub.rho
    beta = sub.Beta
    C = sub.CP
    kappa = sub.Kappa
    r0 = w_beam / np.sqrt(2)
    omega = 2 * np.pi * f

    if exact:
        # arXiv:cond-mat/0402650, Eq. E7
        w = omega * r0**2 * rho * C / (2 * kappa)
        psd = np.abs(
            H
            * beta**2
            * kBT
            * Temp
            / (2 * np.pi * kappa)
            * (np.exp(1j * w) * exp1(1j * w) + np.exp(-1j * w) * exp1(-1j * w))
        )

    else:
        # arXiv:cond-mat/0402650, Eq. 5.3; P1400084, Eq. 18
        psd = (
            4
            * H
            * beta**2
            * kappa
            * kBT
            * Temp
            / (np.pi * r0**4 * omega**2 * (rho * C) ** 2)
        )

    return psd


def coating_brownian(f, stack, w_beam, power=None, mass=None, **kwargs):
    """Coating brownian noise for a given collection of coating layers
    adapted from gwinc.noises which uses Hong et al . PRD 87, 082001 (2013).

    ***Important Note***
    Inside this function phi is used for denoting the phase shift suffered
    by light in one way propagation through a layer. This is in conflict
    with present nomenclature everywhere else where it is used as loss angle.

    The layers are assumed to be alernating low-n high-n layers, with
    low-n first.

    Inputs:
             f = frequency vector in Hz
        mirror = mirror properties Struct
         wBeam = beam radius (at 1 / e**2 power)
         power = laser power falling on the mirror (W)
         mass = mirror mass (kg) when parsing power

    If the power argument is present and is not None, the amplitude noise due
    to coating brownian noise will be calculated and its effect on the phase
    noise will be added (assuming the susceptibility is that of a free mass)

    ***The following parameters are experimental and unsupported as yet***
    The following optional parameters are available in the Materials object
    to provide separate Bulk and Shear loss angles and to include the effect
    of photoelasticity:
    lossBlown = Coating Bulk Loss Angle of Low Refrac.Index layer @ 100Hz
    lossSlown = Coating Shear Loss Angle of Low Refrac. Index layer @ 100Hz
    lossBhighn = Coating Bulk Loss Angle of High Refrac. Index layer @ 100Hz
    lossShighn = Coating Shear Loss Angle of High Refrac. Index layer @ 100Hz
    lossBlown_slope = Coating Bulk Loss Angle Slope of Low Refrac. Index layer
    lossSlown_slope = Coating Shear Loss Angle Slope of Low Refrac. Index layer
    lossBhighn_slope = Coating Bulk Loss Angle Slope of High Refrac. Index layer
    lossShighn_slope = Coating Shear Loss Angle Slope of High Refrac. Index layer
    PETlown = Relevant component of Photoelastic Tensor of High n layer*
    PEThighn = Relevant component of Photoelastic Tensor of Low n layer*

    Returns:
      SbrZ = Brownian noise spectra for one mirror in m**2 / Hz

    *
    Choice of PETlown and PEThighn can be inspired from sec. A.1. of the paper.
    There, values are chosen to get the longitudnal coefficent of
    photoelasticity as -0.5 for Tantala and -0.27 for Silica.
    These values also need to be added in Materials object.
    *
    If the optional arguments are not present, Phihighn and Philown will be
    used as both Bulk and Shear loss angles and PET coefficients will be set
    to 0.

    """
    wavelength = stack["lam_ref"]
    sub, coat = stack["sub"], stack["thin_films"]
    kBT = k * sub.Temp

    # substrate properties
    Ysub = sub.Y
    pratsub = sub.Sigma
    nsub = sub.Index

    # coating properties
    dN = stack["Ls"]
    nN = stack["ns"][1:-1]
    dOpt = dN * nN / wavelength
    yN = stack["Ys"]
    pratN = stack["sigmas"]
    phiBN = stack["phis"]
    phiSN = stack["phis"]

    try:
        PETs = stack["phels"]
    except KeyError:
        PETs = 0.0

    # Loss angles, no fancy interpretation from shear and bulk and leave constant and equal
    lossB, lossS = (
        lambda f: phiBN,
        lambda f: phiSN,
    )

    nol = len(dOpt)  # Number of layers
    # Coefficient of photoelastic effect
    # PRD 87, 082001 Eq (A6)
    CPE = -0.5 * PETs * nN**3

    #  Adapted from the more general calculation in
    #  Yam, W., Gras, S., & Evans, M. Multimaterial coatings with reduced
    #  thermal noise. Physical Review D, 91(4), 042002.  (2015).
    #  http://doi.org/10.1103/PhysRevD.91.042002 )

    # Dispersion of light in each layer
    lambdaN = wavelength / nN

    # Calculate rho and derivatives of rho
    # with respect to both phi_k and r_j
    rho, dLogRho_dPhik, dLogRho_dRk, r = getCoatReflAndDer(nN, nsub, dOpt)

    # Define the function epsilon as per Eq (25)
    # Split epsilon function as:
    # Epsilon = Ep1 - Ep2 * cos(2k0n(z-zjp1)) - Ep3 * sin(2k0n(z-zjp1))

    # Part 1 of epsilon function
    Ep1 = (nN + CPE) * dLogRho_dPhik[:-1]
    # Part 2 of epsilon function (Prefactor of cosine term)
    Ep2 = CPE * (
        dLogRho_dPhik[:-1] * (1 - r[:-1] ** 2) / (2 * r[:-1])
        - (dLogRho_dPhik[1:] * (1 + r[:-1] ** 2) / (2 * r[:-1]))
    )
    # Part 3 of epsilon function (Prefactor of sine term)
    Ep3 = (1 - r[:-1] ** 2) * CPE * dLogRho_dRk[:-1]

    # Define (1 - Im(epsilon)/2)
    Ip1 = 1 - Ep1.imag / 2  # First part of (1 - Im(epsilon)/2)
    Ip2 = Ep2.imag / 2  # Prefactor to cosine in (1 - Im(epsilon)/2)
    Ip3 = Ep3.imag / 2  # Prefactor to sine in (1 - Im(epsilon)/2)

    # Define transfer functions from bulk and shear noise fields to layer
    # thickness and surface height as per Table I in paper
    C_B = np.sqrt(0.5 * (1 + pratN))
    C_SA = np.sqrt(1 - 2 * pratN)
    D_B = (
        (1 - pratsub - 2 * pratsub**2)
        * yN
        / (np.sqrt(2 * (1 + pratN)) * Ysub)
    )
    D_SA = -(
        (1 - pratsub - 2 * pratsub**2)
        * yN
        / (2 * np.sqrt(1 - 2 * pratN) * Ysub)
    )
    D_SB = (
        np.sqrt(3)
        * (1 - pratN)
        * (1 - pratsub - 2 * pratsub**2)
        * yN
        / (2 * np.sqrt(1 - 2 * pratN) * (1 + pratN) * Ysub)
    )

    # Calculating effective beam area on each layer
    # Assuming the beam radius does not change significantly through the
    # depth of the mirror.
    Aeff = np.pi * (w_beam**2)

    # PSD at single layer with thickness equal to wavelength
    # in the medium Eq (96)
    # Bulk & Shear
    S_Bk = np.zeros((nol, len(f)))
    S_Sk = np.zeros((nol, len(f)))
    for layer in range(nol):
        S_Bk[layer, :] = (
            4
            * kBT
            * lambdaN[layer]
            * lossB(f)[layer]
            * (1 - pratN[layer] - 2 * pratN[layer] ** 2)
            / (3 * np.pi * f * yN[layer] * ((1 - pratN[layer]) ** 2) * Aeff)
        )
        S_Sk[layer, :] = (
            4
            * kBT
            * lambdaN[layer]
            * lossS(f)[layer]
            * (1 - pratN[layer] - 2 * pratN[layer] ** 2)
            / (3 * np.pi * f * yN[layer] * ((1 - pratN[layer]) ** 2) * Aeff)
        )

    # Coefficients q_j from Eq (94
    # See https://dcc.ligo.org/T2000552 for derivation
    k0 = 2 * np.pi / wavelength
    q_Bk = (
        +8 * C_B * (D_B + C_B * Ip1) * Ip3
        + 2 * C_B**2 * Ip2 * Ip3
        + 4
        * (
            2 * D_B**2
            + 4 * C_B * D_B * Ip1
            + C_B**2 * (2 * Ip1**2 + Ip2**2 + Ip3**2)
        )
        * k0
        * nN
        * dN
        - 8 * C_B * (D_B + C_B * Ip1) * Ip3 * np.cos(2 * k0 * nN * dN)
        - 2 * C_B**2 * Ip2 * Ip3 * np.cos(4 * k0 * nN * dN)
        + 8 * C_B * (D_B + C_B * Ip1) * Ip2 * np.sin(2 * k0 * nN * dN)
        + C_B**2 * (Ip2 - Ip3) * (Ip2 + Ip3) * np.sin(4 * k0 * nN * dN)
    ) / (8 * k0 * lambdaN * nN)

    q_Sk = (
        D_SB**2 * 8 * k0 * nN * dN
        + 8 * C_SA * (D_SA + C_SA * Ip1) * Ip3
        + 2 * C_SA**2 * Ip2 * Ip3
        + 4
        * (
            2 * D_SA**2
            + 4 * C_SA * D_SA * Ip1
            + C_SA**2 * (2 * Ip1**2 + Ip2**2 + Ip3**2)
        )
        * k0
        * nN
        * dN
        - 8 * C_SA * (D_SA + C_SA * Ip1) * Ip3 * np.cos(2 * k0 * nN * dN)
        - 2 * C_SA**2 * Ip2 * Ip3 * np.cos(4 * k0 * nN * dN)
        + 8 * C_SA * (D_SA + C_SA * Ip1) * Ip2 * np.sin(2 * k0 * nN * dN)
        + C_SA**2 * (Ip2 - Ip3) * (Ip2 + Ip3) * np.sin(4 * k0 * nN * dN)
    ) / (8 * k0 * lambdaN * nN)

    # S_Xi as per Eq(94)
    S_Xi = np.tensordot(q_Bk, S_Bk, axes=1) + np.tensordot(q_Sk, S_Sk, axes=1)

    # From Sec II.E. Eq.(41)
    # Conversion of brownian amplitude noise to displacement noise
    if power is not None:
        # get/calculate optic transmittance
        mTi = stack["T_ref"]

        # Define Re(epsilon)/2
        Rp1 = Ep1.real / 2  # First part of Re(epsilon)/2
        Rp2 = -Ep2.real / 2  # Prefactor to cosine in Re(epsilon)/2
        Rp3 = -Ep3.real / 2  # Prefactor to sine in Re(epsilon)/2
        # Coefficients p_j from Eq (95)
        # See https://dcc.ligo.org/T2000552 for derivation
        p_BkbyC = (
            +8 * Rp1 * Rp3
            + 2 * Rp2 * Rp3
            + 4 * (2 * Rp1**2 + Rp2**2 + Rp3**2) * k0 * nN * dN
            - 8 * Rp1 * Rp3 * np.cos(2 * k0 * nN * dN)
            - 2 * Rp2 * Rp3 * np.cos(4 * k0 * nN * dN)
            + 8 * Rp1 * Rp2 * np.sin(2 * k0 * nN * dN)
            + (Rp2 - Rp3) * (Rp2 + Rp3) * np.sin(4 * k0 * nN * dN)
        ) / (8 * k0 * lambdaN * nN)
        p_Bk = p_BkbyC * C_B**2
        p_Sk = p_BkbyC * C_SA**2

        # S_Zeta as per Eq(95)
        S_Zeta = np.tensordot(p_Bk, S_Bk, axes=1) + np.tensordot(
            p_Sk, S_Sk, axes=1
        )

        AmpToDispConvFac = (32 * power) / (
            mass * wavelength * f**2 * c * 2 * np.pi * np.sqrt(mTi)
        )
        # Adding the two pathways of noise contribution as correlated noise
        SbrZ = (np.sqrt(S_Xi) + AmpToDispConvFac * np.sqrt(S_Zeta)) ** 2
    else:
        SbrZ = S_Xi

    return SbrZ


def getCoatRefl(materials, dOpt):
    """Amplitude reflectivity, with phase, of a coating

    :materials: gwinc optic materials sturcutre
    :dOpt: coating layer thickness array (Nlayer x 1)

    :returns: see return value of :geteCoatRefl2:

    """
    pS = materials.Substrate
    pC = materials.Coating

    nS = pS.RefractiveIndex
    nL = pC.Indexlown
    nH = pC.Indexhighn

    Nlayer = len(dOpt)

    # refractive index of input, coating, and output materials
    nAll = np.zeros(Nlayer + 2)
    nAll[0] = 1  # vacuum input
    nAll[1::2] = nL
    nAll[2::2] = nH
    nAll[-1] = nS  # substrate output

    # backend calculation
    return getCoatRefl2(nAll[0], nAll[-1], nAll[1:-1], dOpt)


def getCoatRefl2(nIn, nOut, nLayer, dOpt):
    """Coating reflection and phase derivatives

    :nIn: refractive index of input medium (e.g., vacuum = 1)
    :nOut: refractive index of output medium (e.g., SiO2 = 1.45231 @ 1064nm)
    :nLayer: refractive index of each layer, ordered input to output (N x 1)
    :dOpt: optical thickness / lambda of each layer,
           geometrical thickness * refractive index / lambda

    :returns: tuple of:
    rCoat = amplitude reflectivity of coating (complex) = rbar(0)
    dcdp = d reflection phase / d round-trip layer phase
    rbar = amplitude reflectivity of coating from this layer down
    r = amplitude reflectivity of this interface (r(1) is nIn to nLayer(1))

    See LIGO-T080101.

    """
    # Z-dir (1 => away from the substrate, -1 => into the substrate)
    zdir = 1

    # vector of all refractive indexs
    nAll = np.concatenate(([nIn], nLayer, [nOut]))

    # reflectivity of each interface
    r = (nAll[:-1] - nAll[1:]) / (nAll[:-1] + nAll[1:])

    # combine reflectivities
    rbar = np.zeros(r.size, dtype=complex)
    ephi = np.zeros(r.size, dtype=complex)

    # round-trip phase in each layer
    ephi[0] = 1
    ephi[1:] = np.exp(4j * zdir * np.pi * np.asarray(dOpt))

    rbar[-1] = ephi[-1] * r[-1]
    for n in range(len(dOpt), 0, -1):
        # accumulate reflectivity
        rbar[n - 1] = (
            ephi[n - 1] * (r[n - 1] + rbar[n]) / (1 + r[n - 1] * rbar[n])
        )

    # reflectivity derivatives
    dr_dphi = ephi[:-1] * (1 - r[:-1] ** 2) / (1 + r[:-1] * rbar[1:]) ** 2
    dr_dphi = (1j * zdir * rbar[1:]) * np.multiply.accumulate(dr_dphi)

    # shift rbar index
    rCoat = rbar[0]
    rbar = rbar[1:]

    # phase derivatives
    dcdp = -np.imag(dr_dphi / rCoat)  ### Where did this minus come from???

    return rCoat, dcdp, rbar, r


def getCoatReflAndDer(nN, nsub, dOpt):
    """
    Helper function for coating_brownian_hong().
    Follows Hong et al . PRD 87, 082001 (2013) Sec V.A.
    This function calculates derivatives of complex reflectivity of Coating
    with respect to phase shifts through each layer and reflectivities of
    each interface
    Input:

      nN = Refractive indices of coatings layers
    nsub = Refractive Index of Substrate
    dOpt = optical thickness / lambda of each layer,
           geometrical thickness * refractive index / lambda

    Returns:
     delLogRho_delPhik = Partial derivative of log of total effective
                         reflectivity of coating with respect to phase shifts
                         in each layer.
    delLogRho_delReflk = Partial derivative of log of total effective
                         reflectivity of coating with respect to reflectivity
                         of each interface.
    """
    nol = len(dOpt)  # Number of layers in coating
    # Reflectivities and transmitivities
    # r[j] is reflectivity from (j-1)th and (j)th layer interface
    # Here r[0] is reflectivity from air and 0th layer
    # and r[-1] is reflectivity between last layer and substrate
    Refl = np.zeros(nol + 1)
    Refl[0] = (1 - nN[0]) / (1 + nN[0])
    Refl[1:-1] = (nN[:-1] - nN[1:]) / (nN[:-1] + nN[1:])
    Refl[-1] = (nN[-1] - nsub) / (nN[-1] + nsub)
    # Note the shift from nomenclature
    # Phi is reserved for denoting one-way phase shift suffered by light
    # during propagation through a layer
    Phi = np.asarray(dOpt) * 2 * np.pi

    # Define rho_k as reflectivity of
    # k layers starting from (N-k-1)th lyer to (N-1)th layer
    # So rhoN[-1] is reflectivity for  no layers but interface from
    #                                              last layer to substrate
    # rhoN[0] is total complex reflectivity of the coating stack.
    rhoN = np.zeros_like(Refl, dtype=np.complex128)

    phiNmkm1 = np.flip(Phi, axis=0)  # phi_{N-k-1}
    rNmkm1 = np.flip(Refl[:-1], axis=0)  # r_{N-k-1}
    exp2iphiNmkm1 = np.exp(2j * phiNmkm1)  # exp(2i phi_{N-k-1})

    # Recursion relation for complex reflectivity
    # See https://dcc.ligo.org/T2000552 for derivation
    rhoN[0] = Refl[-1]
    for k in range(len(Refl) - 1):
        rhoN[k + 1] = (rNmkm1[k] + exp2iphiNmkm1[k] * rhoN[k]) / (
            1 + exp2iphiNmkm1[k] * rNmkm1[k] * rhoN[k]
        )

    denTerm = (1 + exp2iphiNmkm1 * rNmkm1 * rhoN[:-1]) ** 2

    # Derivatives of rho_{k+1} wrt to rho_{k}, r_{N-k-1} and phi_{N-k-1}
    delRhokp1_delRhok = exp2iphiNmkm1 * (1 - rNmkm1**2) / denTerm
    delRhokp1_delRNmkm1 = np.append(
        1, ((1 - (exp2iphiNmkm1 * rhoN[:-1]) ** 2) / denTerm)
    )
    delRhokp1_delPhiNmkm1 = np.append(0, -2j * rhoN[:-1] * delRhokp1_delRhok)

    # Derivative of rho_{N} wrt to rho_{N-j}
    delRhoN_delRhoNmj = np.append(1, np.cumprod(np.flipud(delRhokp1_delRhok)))

    # Derivative of rho_{N} wrt to r_k and phi_k
    delRho_delRk = -delRhoN_delRhoNmj * np.flipud(delRhokp1_delRNmkm1)
    delRho_delPhik = -delRhoN_delRhoNmj * np.flipud(delRhokp1_delPhiNmkm1)
    delLogRho_delReflk = delRho_delRk / rhoN[-1]
    delLogRho_delPhik = delRho_delPhik / rhoN[-1]
    delLogRho_delPhik[-1] = 0  # Define this as per Eq (26)

    return rhoN[-1], delLogRho_delPhik, delLogRho_delReflk, Refl


def substrate_brownian_FiniteCorr(stack, w_beam, r_mirror, d_mirror):
    """Substrate brownian noise finite-size test mass correction

    :stack: stack dict
    :wBeam: beam radius (at 1 / e^2 power)

    :returns: correction factors tuple:
    cftm = finite mirror correction factor
    aftm = amplitude coefficient for thermal noise:
           thermal noise contribution to displacement noise is
           S_x(f) = (8 * kB * T / (2*pi*f)) * Phi(f) * aftm

    Equation references to Bondu, et al. Physics Letters A 246 (1998)
    227-236 (hereafter BHV) and Liu and Thorne gr-qc/0002055 (hereafter LT)

    """
    sub = stack["sub"]
    Y = sub.Y
    sigma = sub.Sigma

    # LT uses e-folding of power
    r0 = w_beam / np.sqrt(2)
    km = zeta / r_mirror

    Qm = np.exp(-2 * km * d_mirror)  # LT eq. 35a

    Um = (1 - Qm) * (1 + Qm) + 4 * d_mirror * km * Qm
    Um = Um / (
        (1 - Qm) ** 2 - 4 * (km * d_mirror) ** 2 * Qm
    )  # LT 53 (BHV eq. btwn 29 & 30)

    x = np.exp(-((zeta * r0 / r_mirror) ** 2) / 4)
    s = sum(x / (zeta**2 * j0m))  # LT 57

    x2 = x * x
    U0 = sum(Um * x2 / (zeta * j0m**2))
    U0 = (
        U0 * (1 - sigma) * (1 + sigma) / (np.pi * r_mirror * Y)
    )  # LT 56 (BHV eq. 3)

    p0 = 1 / (np.pi * r_mirror**2)  # LT 28
    DeltaU = (np.pi * d_mirror**2 * p0) ** 2
    DeltaU = DeltaU + 12 * np.pi * d_mirror**2 * p0 * sigma * s
    DeltaU = DeltaU + 72 * (1 - sigma) * s**2
    DeltaU = DeltaU * r_mirror**2 / (6 * np.pi * d_mirror**3 * Y)  # LT 54

    # LT 58 (eq. following BHV 31)
    aftm = DeltaU + U0

    # amplitude coef for infinite TM, LT 59
    # factored out: (8 * kB * T * Phi) / (2 * pi * f)
    aitm = (1 - sigma**2) / (2 * np.sqrt(2 * np.pi) * Y * r0)

    # finite mirror correction
    cftm = aftm / aitm

    return cftm, aftm


def brownian_proxy(stack):
    """Material properties for speedy evaluation of
    coating Brownian noise.

    Args:
        stack (dict): Dict with coating parameters
    Returns:
        gam (float): Prefactor for use in the proxy
                     function for Brownian noise
                     as per E0900068 pg4.
    """
    pattern = stack["pattern"]
    phiC = stack["phis"]
    nC = stack["ns"][1:-1]
    Y_C = stack["Ys"]
    Y_S = stack["sub"].Y

    phiX, nX, YX, Xs = [], [], [], []
    for j, X in enumerate(pattern):
        if X not in Xs:
            phiX.append(phiC[j])
            nX.append(nC[j])
            YX.append(Y_C[j])
            Xs.append(X)

    # Referenced to lowest index material
    n_low = min(nX)
    phi_low = np.take(phiX, np.argwhere(nX == n_low).flatten())[0]
    Y_low = np.take(YX, np.argwhere(nX == n_low).flatten())[0]
    X_low = np.take(Xs, np.argwhere(nX == n_low).flatten())[0]
    gams = {}
    for j, X in enumerate(Xs):
        if X != X_low:
            gams[X] = (
                (phiX[j] / phi_low)
                * (n_low / nX[j])
                * (YX[j] / Y_S + Y_S / YX[j])
                / (Y_low / Y_S + Y_S / Y_low)
            )
    return gams
