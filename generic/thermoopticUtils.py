"""JIT-compiled thermo-optic noise calculation.

Replaces gwinc.noise.coatingthermal.coating_thermooptic for the
optimization hot path. All physics is identical; the speedup comes
from Numba JIT compilation and eliminating redundant sub-calls.

References:
    LIGO-T080101 (Evans, thermo-optic noise model)
    Braginsky & Vyatchanin, PLA 312, 244 (2003) (finite mirror correction)
"""
import numpy as np
import numba
from scipy.special import jn_zeros, jn

# Pre-computed Bessel zeros and J0 values (same as gwinc.const)
_BESSEL_ZEROS = jn_zeros(1, 300)
_J0M = jn(0, _BESSEL_ZEROS)
_kB = 1.380649e-23  # Boltzmann constant, J/K


def extract_ifo_params(ifo):
    """Extract all material parameters from ifo struct into a flat tuple.

    Called once per getMirrorCost evaluation (Python overhead, ~1 us).
    Returns a tuple of scalars that can be passed to the JIT function.
    """
    pS = ifo.Optics.ETM.Substrate
    pC = ifo.Optics.ETM.Coating
    return (
        # Substrate properties
        pS.RefractiveIndex,      # nS
        pS.MirrorY,             # Y_S
        pS.MirrorSigma,         # sigS
        pS.MassAlpha,           # alphaS
        pS.MassCM,              # CM_S
        pS.MassDensity,         # rhoS
        pS.MassKappa,           # K_S
        pS.Temp,                # Temp
        # Coating low-n properties
        pC.Indexlown,           # nL
        pC.Alphalown,           # alphaL
        pC.Betalown,            # betaL
        pC.Ylown,               # Y_L
        pC.Sigmalown,           # sigL
        pC.CVlown,              # C_L
        pC.ThermalDiffusivitylown,  # K_L
        # Coating high-n properties
        pC.Indexhighn,          # nH
        pC.Alphahighn,          # alphaH
        pC.Betahighn,           # betaH
        pC.Yhighn,              # Y_H
        pC.Sigmahighn,          # sigH
        pC.CVhighn,             # C_H
        pC.ThermalDiffusivityhighn,  # K_H
        # Mirror geometry
        ifo.Optics.ETM.MassRadius,      # MassRadius
        ifo.Optics.ETM.MassThickness,   # MassThickness
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
    """JIT-compiled thermo-optic noise PSD at frequency f.

    Combines getCoatLayers, getCoatAvg, getCoatTOPhase, getCoatRefl2,
    getCoatFiniteCorr, getCoatThickCorr, and getCoatThermal into one
    function to eliminate redundant calculations and Python overhead.
    """
    Nlayer = len(dOpt)
    C_S = CM_S * rhoS
    pi = np.pi

    # ================================================================
    # getCoatLayers — build layer property vectors (computed ONCE)
    # ================================================================
    nLayer = np.empty(Nlayer)
    aLayer = np.empty(Nlayer)
    bLayer = np.empty(Nlayer)
    sLayer = np.empty(Nlayer)

    # Expansion ratio helper (inlined)
    ceL = ((1 + sigS) / (1 - sigL)) * ((1 + sigL) / (1 + sigS) + (1 - 2 * sigS) * Y_L / Y_S)
    ceH = ((1 + sigS) / (1 - sigH)) * ((1 + sigH) / (1 + sigS) + (1 - 2 * sigS) * Y_H / Y_S)

    for i in range(Nlayer):
        if i % 2 == 0:  # low-n layer
            nLayer[i] = nL
            aLayer[i] = alphaL * ceL
            bLayer[i] = betaL
            sLayer[i] = alphaL * (1 + sigL) / (1 - sigL)
        else:  # high-n layer
            nLayer[i] = nH
            aLayer[i] = alphaH * ceH
            bLayer[i] = betaH
            sLayer[i] = alphaH * (1 + sigH) / (1 - sigH)

    # Geometrical thickness
    dLayer = np.empty(Nlayer)
    for i in range(Nlayer):
        dLayer[i] = wavelength * dOpt[i] / nLayer[i]

    # ================================================================
    # getCoatAvg — coating average properties (computed ONCE)
    # ================================================================
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

    # ================================================================
    # getCoatRefl2 — reflectivity and phase derivatives
    # ================================================================
    # Vector of all refractive indices: [nIn, nLayer..., nOut]
    N = Nlayer + 2
    nAll = np.empty(N)
    nAll[0] = 1.0  # vacuum
    for i in range(Nlayer):
        nAll[i + 1] = nLayer[i]
    nAll[N - 1] = nS

    # Interface reflectivities
    r = np.empty(N - 1)
    for i in range(N - 1):
        r[i] = (nAll[i] - nAll[i + 1]) / (nAll[i] + nAll[i + 1])

    # Round-trip phase factors
    ephi = np.empty(N - 1, dtype=np.complex128)
    ephi[0] = 1.0 + 0.0j
    for i in range(Nlayer):
        ephi[i + 1] = np.exp(4.0j * pi * dOpt[i])

    # Recursive reflectivity (backward)
    rbar = np.empty(N - 1, dtype=np.complex128)
    rbar[N - 2] = ephi[N - 2] * r[N - 2]
    for n in range(Nlayer, 0, -1):
        rbar[n - 1] = ephi[n - 1] * (r[n - 1] + rbar[n]) / (1.0 + r[n - 1] * rbar[n])

    rCoat = rbar[0]

    # Phase derivatives
    dr_dphi = np.empty(N - 2, dtype=np.complex128)
    for i in range(N - 2):
        dr_dphi[i] = ephi[i] * (1.0 - r[i]**2) / (1.0 + r[i] * rbar[i + 1])**2

    # Cumulative product of dr_dphi
    acc = dr_dphi[0]
    dr_dphi[0] = 1.0j * rbar[1] * acc
    for i in range(1, N - 2):
        acc *= dr_dphi[i]
        dr_dphi[i] = 1.0j * rbar[i + 1] * acc

    # dcdp = -imag(dr_dphi / rCoat)
    dcdp = np.empty(Nlayer)
    for i in range(Nlayer):
        dcdp[i] = -(dr_dphi[i] / rCoat).imag

    # ================================================================
    # getCoatTOPhase — phase derivatives w.r.t. temperature
    # ================================================================
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

    R_coat = abs(rCoat)**2
    T_coat = 1.0 - R_coat

    # Convert phase to meters
    dTR = dphi_TR * wavelength / (4.0 * pi)
    dTE = dphi_TE * wavelength / (4.0 * pi) - aSub * dc

    # ================================================================
    # getCoatFiniteCorr — finite mirror size correction
    # ================================================================
    R_mir = MassRadius
    H_mir = MassThickness

    # Coating sums for finite correction (reuses dOpt directly)
    dL_fc = wavelength * 0.0
    dH_fc = wavelength * 0.0
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

    # Bessel sum (300 terms)
    S1 = 0.0
    S2 = 0.0
    for m in range(len(bessel_zeros)):
        km = bessel_zeros[m] / R_mir
        Qm = np.exp(-2.0 * km * H_mir)
        pm = np.exp(-km**2 * r0**2 / 4.0) / j0m[m]
        denom = (1.0 - Qm)**2 - 4.0 * km**2 * H_mir**2 * Qm
        Lm = Xr - Zf * (1 + sigS) + (Yr * (1 - 2 * sigS) + Zf - 2 * Cr) * \
             (1 + sigS) * (1 - Qm)**2 / denom
        S1 += pm / bessel_zeros[m]**2
        S2 += pm**2 * Lm**2

    S1 *= 12.0 * R_mir**2 / H_mir**2
    P = (Xr - 2 * sigS * Yr - Cr + S1 * (Cr - Yr * (1 - sigS)))**2 + S2
    LAMBDA = -Cr + (Xr / (1 + sigS) + Yr * (1 - 2 * sigS)) / 2.0
    Cfsm = np.sqrt((r0**2 * P) / (2.0 * R_mir**2 * (1 + sigS)**2 * LAMBDA**2))

    # Apply finite size correction
    dTE = dTE * Cfsm
    dTO = dTE + dTR

    # ================================================================
    # getCoatThickCorr — finite coating thickness correction (x3)
    # ================================================================
    w = 2.0 * pi * f
    R_tc = np.sqrt(Cc * Kc / (C_S * K_S))
    xi = dc * np.sqrt(2.0 * w * Cc / Kc)

    s = np.sin(xi)
    c = np.cos(xi)
    sh = np.sinh(xi)
    ch = np.cosh(xi)

    g0 = 2.0 * (sh - s) + 2.0 * R_tc * (ch - c)
    g1 = 8.0 * np.sin(xi / 2.0) * (R_tc * np.cosh(xi / 2.0) + np.sinh(xi / 2.0))
    g2 = (1.0 + R_tc**2) * sh + (1.0 - R_tc**2) * s + 2.0 * R_tc * ch
    gD = (1.0 + R_tc**2) * ch + (1.0 - R_tc**2) * c + 2.0 * R_tc * sh

    # gTO: both dTE and dTR
    pR_TO = dTR / (dTR + dTE)
    pE_TO = dTE / (dTR + dTE)
    gTO = (pE_TO**2 * g0 + pE_TO * pR_TO * xi * g1 + pR_TO**2 * xi**2 * g2) / (R_tc * xi**2 * gD)

    # gTE: dTR = 0
    # pR = 0, pE = 1
    gTE = g0 / (R_tc * xi**2 * gD)

    # gTR: dTE = 0
    # pR = 1, pE = 0
    gTR = g2 / gD

    # ================================================================
    # getCoatThermal — thermal source spectrum
    # ================================================================
    kBT2 = _kB * Temp**2
    rdel = np.sqrt(2.0 * K_S / (C_S * w))
    SsurfT = 4.0 * kBT2 / (pi * w * C_S * rdel * wBeam**2)

    # ================================================================
    # Final: combine everything
    # ================================================================
    StoZ = SsurfT * gTO * dTO**2

    return StoZ


def coating_thermooptic_fast(f, dOpt, wavelength, wBeam, ifo_params):
    """Fast thermo-optic noise calculation.

    Parameters
    ----------
    f : float
        Frequency in Hz.
    dOpt : array_like
        Optical thicknesses of the coating layers.
    wavelength : float
        Laser wavelength in meters.
    wBeam : float
        Beam radius at 1/e^2 power in meters.
    ifo_params : tuple
        Material parameters from extract_ifo_params().

    Returns
    -------
    StoZ : float
        Thermo-optic displacement noise PSD at frequency f.
    """
    return _coating_thermooptic_jit(
        f, np.asarray(dOpt, dtype=np.float64), wavelength, wBeam,
        *ifo_params,
        _BESSEL_ZEROS, _J0M,
    )
