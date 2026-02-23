"""Thermal noise models for dielectric mirror coatings and substrates.

Combines two complementary implementations:

- **JIT thermo-optic** (Numba): fast single-frequency evaluator for
  the optimization hot path (~10 us/call).
- **Pure-numpy coating & substrate noise**: full spectral models for
  analysis and plotting, including Hong et al. Brownian noise.

All public functions take a *stack* dict (built by :func:`OptimalBragg.qw_stack`)
plus beam/mirror geometry parameters.  No gwinc dependency.

References
----------
.. [Evans]  LIGO-T080101, thermo-optic noise model
.. [BV2003] Braginsky & Vyatchanin, PLA 312, 244 (2003)
.. [Hong]   Hong et al., PRD 87, 082001 (2013)
.. [LT]     Liu & Thorne, gr-qc/0002055
.. [BHV]    Bondu et al., Phys. Lett. A 246, 227 (1998)
.. [E0900068] LIGO-E0900068, Brownian noise proxy
"""

import numpy as np
import numba
from scipy.constants import k as kB
from scipy.constants import c as c_light
from scipy.special import exp1, jn, jn_zeros

# Pre-computed Bessel zeros and J0 values (300 terms)
_BESSEL_ZEROS = jn_zeros(1, 300)
_J0M = jn(0, _BESSEL_ZEROS)


# =====================================================================
#  Numba JIT thermo-optic noise (optimization hot path)
# =====================================================================

def extract_stack_params(stack, r_mirror, d_mirror):
    """Extract material parameters from stack dict into a flat tuple.

    Called once per cost-function evaluation (Python overhead, ~1 us).
    The returned tuple is passed directly to :func:`coating_thermooptic_fast`.

    Parameters
    ----------
    stack : dict
        Stack dict from :func:`OptimalBragg.qw_stack`.
    r_mirror : float
        Mirror radius [m].
    d_mirror : float
        Mirror thickness [m].

    Returns
    -------
    tuple
        Flat tuple of scalars for the JIT function.
    """
    sub = stack["sub"]
    tf = stack["thin_films"]
    L_mat, H_mat = tf["L"], tf["H"]
    return (
        # Substrate
        sub.Index, sub.Y, sub.Sigma, sub.Alpha,
        sub.MassCM, sub.MassDensity, sub.MassKappa, sub.Temp,
        # Coating low-n
        L_mat.Index, L_mat.Alpha, L_mat.Beta,
        L_mat.Y, L_mat.Sigma, L_mat.CV, L_mat.ThermalDiffusivity,
        # Coating high-n
        H_mat.Index, H_mat.Alpha, H_mat.Beta,
        H_mat.Y, H_mat.Sigma, H_mat.CV, H_mat.ThermalDiffusivity,
        # Mirror geometry
        r_mirror, d_mirror,
    )


@numba.njit(cache=True)
def _coating_thermooptic_jit(
    f, dOpt, wavelength, wBeam,
    # Substrate
    nS, Y_S, sigS, alphaS, CM_S, rhoS, K_S, Temp,
    # Coating low-n
    nL, alphaL, betaL, Y_L, sigL, C_L, K_L,
    # Coating high-n
    nH, alphaH, betaH, Y_H, sigH, C_H, K_H,
    # Mirror geometry
    MassRadius, MassThickness,
    # Pre-computed Bessel arrays
    bessel_zeros, j0m,
):
    """JIT-compiled thermo-optic noise PSD at frequency *f*.

    Combines getCoatLayers, getCoatAvg, getCoatTOPhase, getCoatRefl2,
    getCoatFiniteCorr, getCoatThickCorr, and getCoatThermal into one
    function to eliminate redundant calculations and Python overhead.
    """
    Nlayer = len(dOpt)
    C_S = CM_S * rhoS
    pi = np.pi

    # getCoatLayers
    nLayer = np.empty(Nlayer)
    aLayer = np.empty(Nlayer)
    bLayer = np.empty(Nlayer)
    sLayer = np.empty(Nlayer)

    ceL = ((1 + sigS) / (1 - sigL)) * (
        (1 + sigL) / (1 + sigS) + (1 - 2 * sigS) * Y_L / Y_S
    )
    ceH = ((1 + sigS) / (1 - sigH)) * (
        (1 + sigH) / (1 + sigS) + (1 - 2 * sigS) * Y_H / Y_S
    )

    for i in range(Nlayer):
        if i % 2 == 0:
            nLayer[i] = nL
            aLayer[i] = alphaL * ceL
            bLayer[i] = betaL
            sLayer[i] = alphaL * (1 + sigL) / (1 - sigL)
        else:
            nLayer[i] = nH
            aLayer[i] = alphaH * ceH
            bLayer[i] = betaH
            sLayer[i] = alphaH * (1 + sigH) / (1 - sigH)

    dLayer = np.empty(Nlayer)
    for i in range(Nlayer):
        dLayer[i] = wavelength * dOpt[i] / nLayer[i]

    # getCoatAvg
    dc = 0.0
    dL_sum = 0.0
    dH_sum = 0.0
    for i in range(Nlayer):
        dc += dLayer[i]
        if i % 2 == 0:
            dL_sum += dLayer[i]
        else:
            dH_sum += dLayer[i]

    Cc = (C_L * dL_sum + C_H * dH_sum) / dc
    Kc = dc / (dL_sum / K_L + dH_sum / K_H)
    aSub = 2 * alphaS * (1 + sigS) * Cc / C_S

    # getCoatRefl2
    N = Nlayer + 2
    nAll = np.empty(N)
    nAll[0] = 1.0
    for i in range(Nlayer):
        nAll[i + 1] = nLayer[i]
    nAll[N - 1] = nS

    r = np.empty(N - 1)
    for i in range(N - 1):
        r[i] = (nAll[i] - nAll[i + 1]) / (nAll[i] + nAll[i + 1])

    ephi = np.empty(N - 1, dtype=np.complex128)
    ephi[0] = 1.0 + 0.0j
    for i in range(Nlayer):
        ephi[i + 1] = np.exp(4.0j * pi * dOpt[i])

    rbar = np.empty(N - 1, dtype=np.complex128)
    rbar[N - 2] = ephi[N - 2] * r[N - 2]
    for n in range(Nlayer, 0, -1):
        rbar[n - 1] = ephi[n - 1] * (r[n - 1] + rbar[n]) / (
            1.0 + r[n - 1] * rbar[n]
        )

    rCoat = rbar[0]

    dr_dphi = np.empty(N - 2, dtype=np.complex128)
    for i in range(N - 2):
        dr_dphi[i] = ephi[i] * (1.0 - r[i] ** 2) / (
            1.0 + r[i] * rbar[i + 1]
        ) ** 2

    acc = dr_dphi[0]
    dr_dphi[0] = 1.0j * rbar[1] * acc
    for i in range(1, N - 2):
        acc *= dr_dphi[i]
        dr_dphi[i] = 1.0j * rbar[i + 1] * acc

    dcdp = np.empty(Nlayer)
    for i in range(Nlayer):
        dcdp[i] = -(dr_dphi[i] / rCoat).imag

    # getCoatTOPhase
    dGeo = np.empty(Nlayer)
    for i in range(Nlayer):
        dGeo[i] = dOpt[i] / nLayer[i]

    dphi_dd = np.empty(Nlayer)
    for i in range(Nlayer):
        dphi_dd[i] = 4.0 * pi * dcdp[i]

    dphi_TR = 0.0
    dphi_TE = 0.0
    for i in range(Nlayer):
        dphi_TR += dphi_dd[i] * (bLayer[i] + sLayer[i] * nLayer[i]) * dGeo[i]
        dphi_TE += aLayer[i] * dGeo[i]
    dphi_TE *= 4.0 * pi

    dTR = dphi_TR * wavelength / (4.0 * pi)
    dTE = dphi_TE * wavelength / (4.0 * pi) - aSub * dc

    # getCoatFiniteCorr
    R_mir = MassRadius
    H_mir = MassThickness

    dL_fc = 0.0
    dH_fc = 0.0
    for i in range(Nlayer):
        if i % 2 == 0:
            dL_fc += dOpt[i]
        else:
            dH_fc += dOpt[i]
    dL_fc *= wavelength / nL
    dH_fc *= wavelength / nH
    dc_fc = dL_fc + dH_fc

    Cf = (C_L * dL_fc + C_H * dH_fc) / dc_fc
    Cr = Cf / C_S

    xxL = alphaL * (1 + sigL) / (1 - sigL)
    xxH = alphaH * (1 + sigH) / (1 - sigH)
    Xf = (xxL * dL_fc + xxH * dH_fc) / dc_fc
    Xr = Xf / alphaS

    yyL = alphaL * Y_L / (1 - sigL)
    yyH = alphaH * Y_H / (1 - sigH)
    Yf = (yyL * dL_fc + yyH * dH_fc) / dc_fc
    Yr = Yf / (alphaS * Y_S)

    zzL = 1.0 / (1 - sigL)
    zzH = 1.0 / (1 - sigH)
    Zf = (zzL * dL_fc + zzH * dH_fc) / dc_fc

    r0 = wBeam / np.sqrt(2.0)

    S1 = 0.0
    S2 = 0.0
    for m in range(len(bessel_zeros)):
        km = bessel_zeros[m] / R_mir
        Qm = np.exp(-2.0 * km * H_mir)
        pm = np.exp(-km ** 2 * r0 ** 2 / 4.0) / j0m[m]
        denom = (1.0 - Qm) ** 2 - 4.0 * km ** 2 * H_mir ** 2 * Qm
        Lm = (
            Xr - Zf * (1 + sigS)
            + (Yr * (1 - 2 * sigS) + Zf - 2 * Cr)
            * (1 + sigS) * (1 - Qm) ** 2 / denom
        )
        S1 += pm / bessel_zeros[m] ** 2
        S2 += pm ** 2 * Lm ** 2

    S1 *= 12.0 * R_mir ** 2 / H_mir ** 2
    P = (Xr - 2 * sigS * Yr - Cr + S1 * (Cr - Yr * (1 - sigS))) ** 2 + S2
    LAMBDA = -Cr + (Xr / (1 + sigS) + Yr * (1 - 2 * sigS)) / 2.0
    Cfsm = np.sqrt(
        (r0 ** 2 * P) / (2.0 * R_mir ** 2 * (1 + sigS) ** 2 * LAMBDA ** 2)
    )

    dTE = dTE * Cfsm
    dTO = dTE + dTR

    # getCoatThickCorr
    w = 2.0 * pi * f
    R_tc = np.sqrt(Cc * Kc / (C_S * K_S))
    xi = dc * np.sqrt(2.0 * w * Cc / Kc)

    s = np.sin(xi)
    c = np.cos(xi)
    sh = np.sinh(xi)
    ch = np.cosh(xi)

    g0 = 2.0 * (sh - s) + 2.0 * R_tc * (ch - c)
    g1 = 8.0 * np.sin(xi / 2.0) * (
        R_tc * np.cosh(xi / 2.0) + np.sinh(xi / 2.0)
    )
    g2 = (1.0 + R_tc ** 2) * sh + (1.0 - R_tc ** 2) * s + 2.0 * R_tc * ch
    gD = (1.0 + R_tc ** 2) * ch + (1.0 - R_tc ** 2) * c + 2.0 * R_tc * sh

    pR_TO = dTR / (dTR + dTE)
    pE_TO = dTE / (dTR + dTE)
    gTO = (
        pE_TO ** 2 * g0 + pE_TO * pR_TO * xi * g1 + pR_TO ** 2 * xi ** 2 * g2
    ) / (R_tc * xi ** 2 * gD)

    # getCoatThermal
    kBT2 = 1.380649e-23 * Temp ** 2
    rdel = np.sqrt(2.0 * K_S / (C_S * w))
    SsurfT = 4.0 * kBT2 / (pi * w * C_S * rdel * wBeam ** 2)

    # Final
    StoZ = SsurfT * gTO * dTO ** 2
    return StoZ


def coating_thermooptic_fast(f, dOpt, wavelength, w_beam, stack_params):
    """Fast JIT thermo-optic noise for the optimization hot path.

    Parameters
    ----------
    f : float
        Frequency [Hz].
    dOpt : array_like
        Optical thicknesses of coating layers.
    wavelength : float
        Laser wavelength [m].
    w_beam : float
        Beam radius at 1/e^2 power [m].
    stack_params : tuple
        Material parameters from :func:`extract_stack_params`.

    Returns
    -------
    float
        Thermo-optic displacement noise PSD [m^2/Hz].
    """
    return _coating_thermooptic_jit(
        f, np.asarray(dOpt, dtype=np.float64), wavelength, w_beam,
        *stack_params, _BESSEL_ZEROS, _J0M,
    )


# =====================================================================
#  Helper functions (pure numpy, from OptimalBragg)
# =====================================================================

def getCoatLayers(stack):
    """Layer vectors for refractive index, effective alpha, beta, thickness.

    Returns (nLayer, aLayer, bLayer, dLayer, sLayer).
    """
    wavelength = stack["lam_ref"]
    nC = stack["ns"][1:-1]
    dOpt = stack["Ls"] * nC / wavelength

    Y_S = stack["sub"].Y
    sigS = stack["sub"].Sigma
    Y_C = stack["Ys"]
    cteC = stack["ctes"]
    sigmaC = stack["sigmas"]
    betaC = stack["betas"]

    def _expansion_ratio(Y_C, sigC, Y_S, sigS):
        return ((1 + sigS) / (1 - sigC)) * (
            (1 + sigC) / (1 + sigS) + (1 - 2 * sigS) * Y_C / Y_S
        )

    aLayer = cteC * _expansion_ratio(Y_C, sigmaC, Y_S, sigS)
    bLayer = betaC
    dLayer = wavelength * np.asarray(dOpt) / nC
    sLayer = cteC * (1 + sigmaC) / (1 - sigmaC)

    return nC, aLayer, bLayer, dLayer, sLayer


def getCoatAvg(stack):
    """Coating average properties.

    Returns (dc, Cc, Kc, aSub).
    """
    wavelength = stack["lam_ref"]
    alphaS = stack["sub"].Alpha
    C_S = stack["sub"].CP * stack["sub"].rho
    sigS = stack["sub"].Sigma

    Cv_C = stack["Cvs"]
    Kd_C = stack["thermaldiffs"]
    dLayer = stack["Ls"]

    dcsum = np.sum(dLayer)
    Cc = np.sum(Cv_C * dLayer) / dcsum
    Kc = dcsum / np.sum(dLayer / Kd_C)
    aSub = 2 * alphaS * (1 + sigS) * Cc / C_S

    return dcsum, Cc, Kc, aSub


def getCoatRefl2(nIn, nOut, nLayer, dOpt):
    """Coating reflection and phase derivatives.

    Returns (rCoat, dcdp, rbar, r).  See LIGO-T080101.
    """
    nAll = np.concatenate(([nIn], nLayer, [nOut]))
    r = (nAll[:-1] - nAll[1:]) / (nAll[:-1] + nAll[1:])

    rbar = np.zeros(r.size, dtype=complex)
    ephi = np.zeros(r.size, dtype=complex)
    ephi[0] = 1
    ephi[1:] = np.exp(4j * np.pi * np.asarray(dOpt))

    rbar[-1] = ephi[-1] * r[-1]
    for n in range(len(dOpt), 0, -1):
        rbar[n - 1] = (
            ephi[n - 1] * (r[n - 1] + rbar[n]) / (1 + r[n - 1] * rbar[n])
        )

    dr_dphi = ephi[:-1] * (1 - r[:-1] ** 2) / (1 + r[:-1] * rbar[1:]) ** 2
    dr_dphi = 1j * rbar[1:] * np.multiply.accumulate(dr_dphi)

    rCoat = rbar[0]
    rbar = rbar[1:]
    dcdp = -np.imag(dr_dphi / rCoat)

    return rCoat, dcdp, rbar, r


def getCoatReflAndDer(nN, nsub, dOpt):
    """Coating reflectivity derivatives for Hong Brownian noise.

    Returns (rho_total, dLogRho_dPhik, dLogRho_dReflk, Refl).
    See Hong et al. PRD 87, 082001 (2013), Sec. V.A.
    """
    nol = len(dOpt)
    Refl = np.zeros(nol + 1)
    Refl[0] = (1 - nN[0]) / (1 + nN[0])
    Refl[1:-1] = (nN[:-1] - nN[1:]) / (nN[:-1] + nN[1:])
    Refl[-1] = (nN[-1] - nsub) / (nN[-1] + nsub)

    Phi = np.asarray(dOpt) * 2 * np.pi

    rhoN = np.zeros_like(Refl, dtype=np.complex128)
    phiNmkm1 = np.flip(Phi, axis=0)
    rNmkm1 = np.flip(Refl[:-1], axis=0)
    exp2iphiNmkm1 = np.exp(2j * phiNmkm1)

    rhoN[0] = Refl[-1]
    for kk in range(len(Refl) - 1):
        rhoN[kk + 1] = (rNmkm1[kk] + exp2iphiNmkm1[kk] * rhoN[kk]) / (
            1 + exp2iphiNmkm1[kk] * rNmkm1[kk] * rhoN[kk]
        )

    denTerm = (1 + exp2iphiNmkm1 * rNmkm1 * rhoN[:-1]) ** 2
    delRhokp1_delRhok = exp2iphiNmkm1 * (1 - rNmkm1 ** 2) / denTerm
    delRhokp1_delRNmkm1 = np.append(
        1, (1 - (exp2iphiNmkm1 * rhoN[:-1]) ** 2) / denTerm
    )
    delRhokp1_delPhiNmkm1 = np.append(
        0, -2j * rhoN[:-1] * delRhokp1_delRhok
    )

    delRhoN_delRhoNmj = np.append(
        1, np.cumprod(np.flipud(delRhokp1_delRhok))
    )

    delRho_delRk = -delRhoN_delRhoNmj * np.flipud(delRhokp1_delRNmkm1)
    delRho_delPhik = -delRhoN_delRhoNmj * np.flipud(delRhokp1_delPhiNmkm1)
    delLogRho_delReflk = delRho_delRk / rhoN[-1]
    delLogRho_delPhik = delRho_delPhik / rhoN[-1]
    delLogRho_delPhik[-1] = 0

    return rhoN[-1], delLogRho_delPhik, delLogRho_delReflk, Refl


def getCoatTOPhase(nIn, nOut, nLayer, dOpt, aLayer, bLayer, sLayer):
    """Coating reflection phase derivatives w.r.t. temperature.

    Returns (dphi_dT, dphi_TE, dphi_TR, rCoat).  See LIGO-T080101.
    """
    rCoat, dcdp = getCoatRefl2(nIn, nOut, nLayer, dOpt)[:2]
    dGeo = np.asarray(dOpt) / nLayer
    dphi_dd = 4 * np.pi * dcdp

    dphi_TR = np.sum(dphi_dd * (bLayer + sLayer * nLayer) * dGeo)
    dphi_TE = 4 * np.pi * np.sum(aLayer * dGeo)
    dphi_dT = dphi_TR + dphi_TE

    return dphi_dT, dphi_TE, dphi_TR, rCoat


def getCoatThermal(f, stack, w_beam):
    """Thermal source spectrum.

    Returns (SsurfT, rdel).
    """
    pS = stack["sub"]
    C_S = pS.CP * pS.rho
    K_S = pS.Kappa
    kBT2 = kB * pS.Temp ** 2
    w = 2 * np.pi * f
    rdel = np.sqrt(2 * K_S / (C_S * w))
    SsurfT = 4 * kBT2 / (np.pi * w * C_S * rdel * w_beam ** 2)
    return SsurfT, rdel


def getCoatThickCorr(f, stack, dTE, dTR):
    """Finite coating thickness correction (LIGO-T080101)."""
    pS = stack["sub"]
    Cs = pS.CP * pS.rho
    Ks = pS.Kappa
    dc, Cc, Kc, _ = getCoatAvg(stack)

    w = 2 * np.pi * f
    R = np.sqrt(Cc * Kc / (Cs * Ks))
    xi = dc * np.sqrt(2 * w * Cc / Kc)

    s = np.sin(xi)
    c = np.cos(xi)
    sh = np.sinh(xi)
    ch = np.cosh(xi)

    pR = dTR / (dTR + dTE)
    pE = dTE / (dTR + dTE)

    g0 = 2 * (sh - s) + 2 * R * (ch - c)
    g1 = 8 * np.sin(xi / 2) * (R * np.cosh(xi / 2) + np.sinh(xi / 2))
    g2 = (1 + R ** 2) * sh + (1 - R ** 2) * s + 2 * R * ch
    gD = (1 + R ** 2) * ch + (1 - R ** 2) * c + 2 * R * sh

    gTC = (pE ** 2 * g0 + pE * pR * xi * g1 + pR ** 2 * xi ** 2 * g2) / (
        R * xi ** 2 * gD
    )
    return gTC


def getCoatFiniteCorr(stack, w_beam, r_mirror, d_mirror):
    """Finite mirror size correction (Braginsky & Vyatchanin, PLA 2003)."""
    wavelength = stack["lam_ref"]
    nC = stack["ns"][1:-1]
    dC = stack["Ls"]

    alphaS = stack["sub"].Alpha
    C_S = stack["sub"].CP * stack["sub"].rho
    Y_S = stack["sub"].Y
    sigS = stack["sub"].Sigma

    alphaC = stack["ctes"]
    Cv_C = stack["Cvs"]
    Y_C = stack["Ys"]
    sigmaC = stack["sigmas"]

    dcsum = np.sum(dC)
    Cf = np.sum(Cv_C * dC) / dcsum
    Cr = Cf / C_S

    xxC = alphaC * (1 + sigmaC) / (1 - sigmaC)
    Xf = np.sum(xxC * dC) / dcsum
    Xr = Xf / alphaS

    yyC = alphaC * Y_C / (1 - sigmaC)
    Yf = np.sum(yyC * dC) / dcsum
    Yr = Yf / (alphaS * Y_S)

    zzC = 1 / (1 - sigmaC)
    Zf = np.sum(zzC * dC) / dcsum

    r0 = w_beam / np.sqrt(2)
    km = _BESSEL_ZEROS / r_mirror
    Qm = np.exp(-2 * km * d_mirror)
    pm = np.exp(-km ** 2 * r0 ** 2 / 4) / _J0M

    Lm = (
        Xr - Zf * (1 + sigS)
        + (Yr * (1 - 2 * sigS) + Zf - 2 * Cr)
        * (1 + sigS) * (1 - Qm) ** 2
        / ((1 - Qm) ** 2 - 4 * km ** 2 * d_mirror ** 2 * Qm)
    )

    S1 = (12 * r_mirror ** 2 / d_mirror ** 2) * np.sum(pm / _BESSEL_ZEROS ** 2)
    S2 = np.sum(pm ** 2 * Lm ** 2)
    P = (Xr - 2 * sigS * Yr - Cr + S1 * (Cr - Yr * (1 - sigS))) ** 2 + S2
    LAMBDA = -Cr + (Xr / (1 + sigS) + Yr * (1 - 2 * sigS)) / 2

    Cfsm = np.sqrt(
        (r0 ** 2 * P) / (2 * r_mirror ** 2 * (1 + sigS) ** 2 * LAMBDA ** 2)
    )
    return Cfsm


def getCoatTOPos(stack, w_beam, r_mirror=None, d_mirror=None):
    """Mirror position derivative w.r.t. thermal fluctuations.

    Returns (dTO, dTR, dTE, T, R).  See LIGO-T080101.
    """
    wavelength = stack["lam_ref"]
    nS = stack["sub"].Index
    dOpt = stack["Ls"] * stack["ns"][1:-1] / wavelength

    nLayer, aLayer, bLayer, dLayer, sLayer = getCoatLayers(stack)
    dc, Cc, Kc, aSub = getCoatAvg(stack)

    dphi_dT, dphi_TE, dphi_TR, rCoat = getCoatTOPhase(
        1, nS, nLayer, dOpt, aLayer, bLayer, sLayer
    )
    R = abs(rCoat) ** 2
    T = 1 - R

    dTR = dphi_TR * wavelength / (4 * np.pi)
    dTE = dphi_TE * wavelength / (4 * np.pi) - aSub * dc

    if r_mirror is not None and d_mirror is not None:
        Cfsm = getCoatFiniteCorr(stack, w_beam, r_mirror, d_mirror)
    else:
        Cfsm = 1

    dTE = dTE * Cfsm
    dTO = dTE + dTR

    return dTO, dTR, dTE, T, R


# =====================================================================
#  Coating thermo-optic noise (pure numpy, for analysis)
# =====================================================================

def coating_thermooptic(f, stack, w_beam, r_mirror=None, d_mirror=None):
    """Coating thermo-optic displacement noise spectrum.

    Parameters
    ----------
    f : float or array_like
        Frequency [Hz].
    stack : dict
        Stack dict.
    w_beam : float
        Beam radius at 1/e^2 power [m].
    r_mirror, d_mirror : float, optional
        Mirror radius and thickness [m] for finite-size correction.

    Returns
    -------
    StoZ : ndarray
        Total thermo-optic noise PSD [m^2/Hz].
    SteZ : ndarray
        Thermo-elastic component.
    StrZ : ndarray
        Thermo-refractive component.
    """
    dTO, dTR, dTE, T, R = getCoatTOPos(stack, w_beam, r_mirror, d_mirror)

    gTO = getCoatThickCorr(f, stack, dTE, dTR)
    gTE = getCoatThickCorr(f, stack, dTE, 0)
    gTR = getCoatThickCorr(f, stack, 0, dTR)

    SsurfT, _ = getCoatThermal(f, stack, w_beam)

    StoZ = SsurfT * gTO * dTO ** 2
    SteZ = SsurfT * gTE * dTE ** 2
    StrZ = SsurfT * gTR * dTR ** 2

    return StoZ, SteZ, StrZ


# =====================================================================
#  Coating Brownian noise — Hong et al. PRD 87, 082001 (2013)
# =====================================================================

def coating_brownian(f, stack, w_beam, power=None, mass=None):
    """Coating Brownian noise using Hong et al. PRD 87, 082001.

    Parameters
    ----------
    f : float or array_like
        Frequency [Hz].
    stack : dict
        Stack dict.
    w_beam : float
        Beam radius at 1/e^2 power [m].
    power : float, optional
        Laser power on mirror [W].  If given, includes amplitude noise.
    mass : float, optional
        Mirror mass [kg].  Required if *power* is given.

    Returns
    -------
    SbrZ : ndarray
        Brownian noise PSD [m^2/Hz].
    """
    f = np.atleast_1d(np.asarray(f, dtype=float))
    wavelength = stack["lam_ref"]
    sub = stack["sub"]
    kBT = kB * sub.Temp

    Ysub = sub.Y
    pratsub = sub.Sigma
    nsub = sub.Index

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

    lossB = lambda freq: phiBN
    lossS = lambda freq: phiSN

    nol = len(dOpt)
    CPE = -0.5 * PETs * nN ** 3

    lambdaN = wavelength / nN

    rho, dLogRho_dPhik, dLogRho_dRk, r = getCoatReflAndDer(nN, nsub, dOpt)

    # Epsilon function (Eq. 25)
    Ep1 = (nN + CPE) * dLogRho_dPhik[:-1]
    # Guard against r[i]=0 (same-index interfaces, e.g. SiO2 cap next to SiO2 layer):
    # physically, zero Fresnel coefficient means no contribution from that interface.
    r_safe = r[:-1].copy()
    zero_mask = r_safe == 0
    r_safe[zero_mask] = 1.0  # dummy to avoid division by zero
    Ep2 = CPE * (
        dLogRho_dPhik[:-1] * (1 - r_safe ** 2) / (2 * r_safe)
        - dLogRho_dPhik[1:] * (1 + r_safe ** 2) / (2 * r_safe)
    )
    Ep2[zero_mask] = 0.0  # zero contribution from same-index interfaces
    Ep3 = (1 - r[:-1] ** 2) * CPE * dLogRho_dRk[:-1]

    Ip1 = 1 - Ep1.imag / 2
    Ip2 = Ep2.imag / 2
    Ip3 = Ep3.imag / 2

    # Transfer functions (Table I)
    C_B = np.sqrt(0.5 * (1 + pratN))
    C_SA = np.sqrt(1 - 2 * pratN)
    D_B = (
        (1 - pratsub - 2 * pratsub ** 2) * yN
        / (np.sqrt(2 * (1 + pratN)) * Ysub)
    )
    D_SA = -(
        (1 - pratsub - 2 * pratsub ** 2) * yN
        / (2 * np.sqrt(1 - 2 * pratN) * Ysub)
    )
    D_SB = (
        np.sqrt(3) * (1 - pratN)
        * (1 - pratsub - 2 * pratsub ** 2) * yN
        / (2 * np.sqrt(1 - 2 * pratN) * (1 + pratN) * Ysub)
    )

    Aeff = np.pi * w_beam ** 2

    # PSD per layer (Eq. 96)
    S_Bk = np.zeros((nol, len(f)))
    S_Sk = np.zeros((nol, len(f)))
    for layer in range(nol):
        common = (
            4 * kBT * lambdaN[layer]
            * (1 - pratN[layer] - 2 * pratN[layer] ** 2)
            / (3 * np.pi * f * yN[layer] * (1 - pratN[layer]) ** 2 * Aeff)
        )
        S_Bk[layer, :] = common * lossB(f)[layer]
        S_Sk[layer, :] = common * lossS(f)[layer]

    # Coefficients q_j (Eq. 94)
    k0 = 2 * np.pi / wavelength
    q_Bk = (
        8 * C_B * (D_B + C_B * Ip1) * Ip3
        + 2 * C_B ** 2 * Ip2 * Ip3
        + 4 * (
            2 * D_B ** 2 + 4 * C_B * D_B * Ip1
            + C_B ** 2 * (2 * Ip1 ** 2 + Ip2 ** 2 + Ip3 ** 2)
        ) * k0 * nN * dN
        - 8 * C_B * (D_B + C_B * Ip1) * Ip3 * np.cos(2 * k0 * nN * dN)
        - 2 * C_B ** 2 * Ip2 * Ip3 * np.cos(4 * k0 * nN * dN)
        + 8 * C_B * (D_B + C_B * Ip1) * Ip2 * np.sin(2 * k0 * nN * dN)
        + C_B ** 2 * (Ip2 - Ip3) * (Ip2 + Ip3) * np.sin(4 * k0 * nN * dN)
    ) / (8 * k0 * lambdaN * nN)

    q_Sk = (
        D_SB ** 2 * 8 * k0 * nN * dN
        + 8 * C_SA * (D_SA + C_SA * Ip1) * Ip3
        + 2 * C_SA ** 2 * Ip2 * Ip3
        + 4 * (
            2 * D_SA ** 2 + 4 * C_SA * D_SA * Ip1
            + C_SA ** 2 * (2 * Ip1 ** 2 + Ip2 ** 2 + Ip3 ** 2)
        ) * k0 * nN * dN
        - 8 * C_SA * (D_SA + C_SA * Ip1) * Ip3 * np.cos(2 * k0 * nN * dN)
        - 2 * C_SA ** 2 * Ip2 * Ip3 * np.cos(4 * k0 * nN * dN)
        + 8 * C_SA * (D_SA + C_SA * Ip1) * Ip2 * np.sin(2 * k0 * nN * dN)
        + C_SA ** 2 * (Ip2 - Ip3) * (Ip2 + Ip3) * np.sin(4 * k0 * nN * dN)
    ) / (8 * k0 * lambdaN * nN)

    S_Xi = np.tensordot(q_Bk, S_Bk, axes=1) + np.tensordot(q_Sk, S_Sk, axes=1)

    if power is not None:
        from OptimalBragg.layers import multilayer_diel
        r_coat, _ = multilayer_diel(
            stack["ns"], stack["Ls"], wavelength
        )
        mTi = 1 - np.abs(r_coat) ** 2

        Rp1 = Ep1.real / 2
        Rp2 = -Ep2.real / 2
        Rp3 = -Ep3.real / 2

        p_BkbyC = (
            8 * Rp1 * Rp3 + 2 * Rp2 * Rp3
            + 4 * (2 * Rp1 ** 2 + Rp2 ** 2 + Rp3 ** 2) * k0 * nN * dN
            - 8 * Rp1 * Rp3 * np.cos(2 * k0 * nN * dN)
            - 2 * Rp2 * Rp3 * np.cos(4 * k0 * nN * dN)
            + 8 * Rp1 * Rp2 * np.sin(2 * k0 * nN * dN)
            + (Rp2 - Rp3) * (Rp2 + Rp3) * np.sin(4 * k0 * nN * dN)
        ) / (8 * k0 * lambdaN * nN)
        p_Bk = p_BkbyC * C_B ** 2
        p_Sk = p_BkbyC * C_SA ** 2

        S_Zeta = (
            np.tensordot(p_Bk, S_Bk, axes=1)
            + np.tensordot(p_Sk, S_Sk, axes=1)
        )

        AmpToDispConvFac = (32 * power) / (
            mass * wavelength * f ** 2 * c_light * 2 * np.pi * np.sqrt(mTi)
        )
        SbrZ = (np.sqrt(S_Xi) + AmpToDispConvFac * np.sqrt(S_Zeta)) ** 2
    else:
        SbrZ = S_Xi

    return SbrZ


# =====================================================================
#  Substrate noise models
# =====================================================================

def substrate_brownian(f, stack, w_beam, r_mirror=None, d_mirror=None):
    """Substrate Brownian displacement noise PSD.

    Parameters
    ----------
    f : float or array_like
        Frequency [Hz].
    stack : dict
        Stack dict.
    w_beam : float
        Beam radius at 1/e^2 power [m].
    r_mirror, d_mirror : float, optional
        Mirror radius and thickness [m] for finite-size correction.

    Returns
    -------
    ndarray
        Displacement noise PSD [m^2/Hz].
    """
    sub = stack["sub"]
    Y = sub.Y
    sigma = sub.Sigma
    c2 = sub.c2
    n_exp = sub.MechanicalLossExponent
    surf_loss = sub.Alphas  # surface loss limit
    kBT = kB * sub.Temp

    if r_mirror is not None and d_mirror is not None:
        cftm, aftm = _substrate_brownian_finite_corr(
            stack, w_beam, r_mirror, d_mirror
        )
    else:
        cftm, aftm = 0, 0

    # Bulk contribution
    phibulk = c2 * f ** n_exp
    cbulk = 8 * kBT * aftm * phibulk / (2 * np.pi * f)

    # Surface loss contribution
    csurf = surf_loss * (1 - 2 * sigma) / ((1 - sigma) * Y * np.pi * w_beam ** 2)
    csurf *= 8 * kBT / (2 * np.pi * f)

    return csurf + cbulk


def _substrate_brownian_finite_corr(stack, w_beam, r_mirror, d_mirror):
    """Substrate Brownian finite-size correction (Liu & Thorne, BHV).

    Returns (cftm, aftm).
    """
    sub = stack["sub"]
    Y = sub.Y
    sigma = sub.Sigma

    r0 = w_beam / np.sqrt(2)
    km = _BESSEL_ZEROS / r_mirror
    Qm = np.exp(-2 * km * d_mirror)

    Um = (1 - Qm) * (1 + Qm) + 4 * d_mirror * km * Qm
    Um = Um / ((1 - Qm) ** 2 - 4 * (km * d_mirror) ** 2 * Qm)

    x = np.exp(-(_BESSEL_ZEROS * r0 / r_mirror) ** 2 / 4)
    s = np.sum(x / (_BESSEL_ZEROS ** 2 * _J0M))

    x2 = x * x
    U0 = np.sum(Um * x2 / (_BESSEL_ZEROS * _J0M ** 2))
    U0 = U0 * (1 - sigma) * (1 + sigma) / (np.pi * r_mirror * Y)

    p0 = 1 / (np.pi * r_mirror ** 2)
    DeltaU = (np.pi * d_mirror ** 2 * p0) ** 2
    DeltaU += 12 * np.pi * d_mirror ** 2 * p0 * sigma * s
    DeltaU += 72 * (1 - sigma) * s ** 2
    DeltaU *= r_mirror ** 2 / (6 * np.pi * d_mirror ** 3 * Y)

    aftm = DeltaU + U0
    aitm = (1 - sigma ** 2) / (2 * np.sqrt(2 * np.pi) * Y * r0)
    cftm = aftm / aitm

    return cftm, aftm


def substrate_thermoelastic(f, stack, w_beam, r_mirror=None, d_mirror=None):
    """Substrate thermoelastic displacement noise PSD.

    Parameters
    ----------
    f : float or array_like
        Frequency [Hz].
    stack : dict
        Stack dict.
    w_beam : float
        Beam radius at 1/e^2 power [m].
    r_mirror, d_mirror : float, optional
        Mirror radius and thickness [m] for finite-size correction.

    Returns
    -------
    ndarray
        Displacement noise PSD [m^2/Hz].
    """
    sub = stack["sub"]
    sigma = sub.Sigma
    rho = sub.rho
    kappa = sub.Kappa
    alpha = sub.Alpha
    CM = sub.CP
    Temp = sub.Temp
    kBT = kB * Temp

    S = 8 * (1 + sigma) ** 2 * kappa * alpha ** 2 * Temp * kBT
    S /= np.sqrt(2 * np.pi) * (CM * rho) ** 2
    S /= (w_beam / np.sqrt(2)) ** 3

    if r_mirror is not None and d_mirror is not None:
        S *= _substrate_thermoelastic_finite_corr(
            stack, w_beam, r_mirror, d_mirror
        )

    return S / (2 * np.pi * f) ** 2


def _substrate_thermoelastic_finite_corr(stack, w_beam, r_mirror, d_mirror):
    """Substrate thermoelastic finite-size correction (Liu & Thorne, Eq. 46)."""
    sub = stack["sub"]
    sigma = sub.Sigma

    r0 = w_beam / np.sqrt(2)
    km = _BESSEL_ZEROS / r_mirror
    Qm = np.exp(-2 * km * d_mirror)

    pm = (
        np.exp(-km ** 2 * r0 ** 2 / 4) / (np.pi * (r_mirror * _J0M) ** 2)
    )

    c0 = 6 * (r_mirror / d_mirror) ** 2 * np.sum(_J0M * pm / _BESSEL_ZEROS ** 2)
    c1 = -2 * c0 / d_mirror
    p0 = 1 / (np.pi * r_mirror ** 2)
    c1 += p0 / (2 * d_mirror)

    coeff = (1 - Qm) * ((1 - Qm) * (1 + Qm) + 8 * d_mirror * km * Qm)
    coeff += 4 * (d_mirror * km) ** 2 * Qm * (1 + Qm)
    coeff *= km * (pm * _J0M) ** 2 * (1 - Qm)
    coeff /= ((1 - Qm) ** 2 - 4 * (d_mirror * km) ** 2 * Qm) ** 2
    coeff = np.sum(coeff) + d_mirror * c1 ** 2 / (1 + sigma) ** 2
    coeff *= (np.sqrt(2 * np.pi) * r0) ** 3 * r_mirror ** 2

    return coeff


def substrate_thermorefractive(f, stack, w_beam, d_mirror):
    """Substrate thermorefractive displacement noise PSD.

    Parameters
    ----------
    f : float or array_like
        Frequency [Hz].
    stack : dict
        Stack dict.
    w_beam : float
        Beam radius at 1/e^2 power [m].
    d_mirror : float
        Mirror thickness [m].

    Returns
    -------
    ndarray
        Displacement noise PSD [m^2/Hz].
    """
    sub = stack["sub"]
    H = d_mirror
    kBT = kB * sub.Temp
    Temp = sub.Temp
    rho = sub.rho
    beta = sub.Beta
    C = sub.CP
    kappa = sub.Kappa
    r0 = w_beam / np.sqrt(2)
    omega = 2 * np.pi * f

    # Adiabatic approximation (cond-mat/0402650, Eq. 5.3)
    psd = (
        4 * H * beta ** 2 * kappa * kBT * Temp
        / (np.pi * r0 ** 4 * omega ** 2 * (rho * C) ** 2)
    )
    return psd


# =====================================================================
#  Brownian noise proxy (fast estimate for optimization)
# =====================================================================

def brownian_proxy(stack):
    """Brownian noise proxy prefactors (LIGO-E0900068).

    Parameters
    ----------
    stack : dict
        Stack dict.

    Returns
    -------
    dict
        Mapping from pattern character to gamma ratio.
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

    n_low = min(nX)
    idx_low = nX.index(n_low)
    phi_low = phiX[idx_low]
    Y_low = YX[idx_low]
    X_low = Xs[idx_low]

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


# =====================================================================
#  Convenience wrappers
# =====================================================================

def coating_noise(f, stack, w_beam, power=None, r_mirror=None,
                  d_mirror=None, m_mirror=None):
    """All coating noise components.

    Returns (coat_Sbr, coat_Ste, coat_Str, coat_Sto).
    """
    coat_Sbr = coating_brownian(
        f=f, stack=stack, w_beam=w_beam, power=power, mass=m_mirror,
    )
    coat_Sto, coat_Ste, coat_Str = coating_thermooptic(
        f=f, stack=stack, w_beam=w_beam, r_mirror=r_mirror, d_mirror=d_mirror,
    )
    return coat_Sbr, coat_Ste, coat_Str, coat_Sto


def substrate_noise(f, stack, w_beam, r_mirror=None, d_mirror=None):
    """All substrate noise components.

    Returns (sub_Sbr, sub_Str, sub_Ste).
    """
    sub_Sbr = substrate_brownian(
        f=f, stack=stack, w_beam=w_beam, r_mirror=r_mirror, d_mirror=d_mirror,
    )
    sub_Str = substrate_thermorefractive(
        f=f, stack=stack, w_beam=w_beam, d_mirror=d_mirror,
    )
    sub_Ste = substrate_thermoelastic(
        f=f, stack=stack, w_beam=w_beam, r_mirror=r_mirror, d_mirror=d_mirror,
    )
    return sub_Sbr, sub_Str, sub_Ste
