# Set of functions used to evaluate a cost function for optimization

import copy
import numpy as np
from generic.coatingUtils import *
from generic.thermoopticUtils import coating_thermooptic_fast, extract_ifo_params
import os
import gwinc.noise

def transmissionCost(target, n, L, lamb=1, theta=0, pol='te'):
    """Evaluate the transmission cost for a dielectric coating.

    Computes the power transmission of the stack and returns a
    normalized squared-error cost relative to the target.

    Parameters
    ----------
    target : float
        Target power transmission.
    n : array_like
        Refractive indices, including incident and transmitted media.
    L : array_like
        Optical thicknesses of the dielectric stack.
    lamb : float or array_like, optional
        Wavelength(s) normalized to design wavelength. Default is 1.
    theta : float, optional
        Angle of incidence in degrees. Default is 0.
    pol : {'te', 'tm'}, optional
        Polarization. Default is ``'te'``.

    Returns
    -------
    cost : float
        Normalized squared error: ``|target - T|^2 / target^2``.
    T : float or ndarray
        Power transmission of the coating.
    """
    r, _ = multidiel1(n, L, lamb, theta, pol)
    T = 1 - np.abs(r)**2
    return (np.abs((target - T)/target)**2)[0], T


def sensitivityCost(target, n, L, lamb=1, theta=0, pol='te'):
    """Evaluate the sensitivity of transmission to a 1% thickness error.

    Perturbs all layer thicknesses by +1% and computes the transmission
    cost at the perturbed thicknesses.

    Parameters
    ----------
    target : float
        Target power transmission.
    n : array_like
        Refractive indices, including incident and transmitted media.
    L : array_like
        Optical thicknesses of the dielectric stack.
    lamb : float or array_like, optional
        Wavelength(s) normalized to design wavelength. Default is 1.
    theta : float, optional
        Angle of incidence in degrees. Default is 0.
    pol : {'te', 'tm'}, optional
        Polarization. Default is ``'te'``.

    Returns
    -------
    cost : float
        Transmission cost evaluated at ``1.01 * L``.
    """
    return transmissionCost(target, n, 1.01*L, lamb, theta, pol)[0]


def surfEfieldCost(target, n, L, lamb=1, theta=0, pol='te'):
    """Evaluate the surface electric field cost.

    Uses an ``arcsinh`` mapping to penalize high surface E-fields
    while remaining smooth near zero.

    Parameters
    ----------
    target : float
        Not used directly; the cost is ``50 * arcsinh(|1 + r|^2)``.
    n : array_like
        Refractive indices, including incident and transmitted media.
    L : array_like
        Optical thicknesses of the dielectric stack.
    lamb : float or array_like, optional
        Wavelength(s) normalized to design wavelength. Default is 1.
    theta : float, optional
        Angle of incidence in degrees. Default is 0.
    pol : {'te', 'tm'}, optional
        Polarization. Default is ``'te'``.

    Returns
    -------
    cost : float
        Surface E-field cost: ``50 * arcsinh(|1 + r|^2)``.
    """
    r, _ = multidiel1(n, L, lamb, theta, pol)
    return (50 * np.arcsinh(np.abs(1 + r)**2))[0]


def stdevLCost(target, L,):
    """Evaluate layer thickness uniformity cost.

    Penalizes deviations of the mean-to-standard-deviation ratio from
    a target value. Uniform stacks (zero std) get ``relative_stdev = 0``.

    Parameters
    ----------
    target : float
        Target value of ``mean(L) / std(L)``.
    L : array_like
        Optical thicknesses of the dielectric stack.

    Returns
    -------
    cost : float
        Normalized squared error: ``|target - relative_stdev|^2 / target^2``.
    """
    if np.std(np.array(L)):
        relative_stdev = np.mean(np.array(L)) / np.std(np.array(L))
    else:
        relative_stdev = 0.0
    return np.abs((target - relative_stdev) / target)**2


def brownianProxy(ifo):
    """Pre-compute the Brownian noise proxy factor gamma.

    Evaluates material properties from the ``ifo`` structure to build a
    fast proxy for coating Brownian thermal noise, per LIGO-E0900068 p. 4.

    Parameters
    ----------
    ifo : gwinc.Struct
        Interferometer model containing material properties.

    Returns
    -------
    gam : float
        Brownian noise proxy pre-factor.
    """
    phi_high = ifo.Materials.Coating.Phihighn
    n_high = ifo.Materials.Coating.Indexhighn
    Y_high = ifo.Materials.Coating.Yhighn
    phi_low = ifo.Materials.Coating.Philown
    n_low = ifo.Materials.Coating.Indexlown
    Y_low = ifo.Materials.Coating.Ylown
    Y_sub = ifo.Materials.Substrate.MirrorY
    a = phi_high / phi_low
    b = n_low / n_high
    c = Y_high / Y_sub + Y_sub / Y_high
    d = Y_low / Y_sub + Y_sub / Y_low
    gam = a*b*c/d
    return(gam)


def brownianCost(target, L, gam,):
    """Compute Brownian noise proxy cost for a coating.

    Uses the formula from LIGO-E0900068 p. 4, normalized by a
    reference value (e.g. quarter-wave stack proxy) so the cost is O(1).

    Parameters
    ----------
    target : float
        Reference Brownian proxy value for normalization (e.g. QW stack).
    L : array_like
        Optical thicknesses of the dielectric stack.
    gam : float
        Pre-computed Brownian proxy factor from :func:`brownianProxy`.

    Returns
    -------
    cost : float
        Normalized Brownian noise proxy: ``SBrZ / target``.
    """
    zLow = np.sum(L[::2])    # Sum of thicknesses of low index layers
    zHigh = np.sum(L[1::2])   # Sum of thicknesses of high index layers
    SBrZ = zLow + gam*zHigh  # Proxy brownian noise
    return SBrZ / target


def thermoopticCost(target, fTarget, L, ifo, ifo_params=None):
    """Compute thermo-optic noise cost.

    Uses a Numba JIT-compiled implementation of the thermo-optic noise
    calculation (extracted from gwinc), normalized by a reference value
    (e.g. quarter-wave stack S_TO) so the cost is O(1).

    Parameters
    ----------
    target : float
        Reference S_TO value for normalization (e.g. QW stack at fTarget).
    fTarget : float
        Frequency at which to evaluate TO noise, in Hz.
    L : array_like
        Optical thicknesses of the dielectric stack.
    ifo : gwinc.Struct
        Interferometer model containing material and optic properties.
    ifo_params : tuple, optional
        Pre-computed material parameters from ``extract_ifo_params(ifo)``.
        If None, extracted on each call.

    Returns
    -------
    cost : float
        Normalized thermo-optic noise cost: ``S_TO(fTarget) / target``.
    """
    if ifo_params is None:
        ifo_params = extract_ifo_params(ifo)
    StoZ = coating_thermooptic_fast(
        fTarget, L, ifo.Laser.Wavelength,
        ifo.Optics.ETM.BeamRadius, ifo_params)
    return StoZ / target


def precompute_misc(costs, ifo, misc):
    """Pre-compute cached values for getMirrorCost hot loop.

    Call once before ``differential_evolution`` to populate ``misc``
    with values that are constant across all evaluations, avoiding
    redundant work on every cost-function call.

    Parameters
    ----------
    costs : dict
        Cost specifications (same as for :func:`getMirrorCost`).
    ifo : gwinc.Struct
        Interferometer model.
    misc : dict
        Will be updated in-place with cached keys:
        ``_ifo_n``, ``_active_costs``, ``_wl_array``, ``_wl_map``,
        ``_sens_wl``, ``_ifo_params``, ``_lambdaAUX``.
    """
    Npairs = misc.get('Npairs', 0)
    Ncopies = misc.get('Ncopies', 0)
    Nfixed = misc.get('Nfixed', 0)

    # Total number of layers after copies/fixed
    nLayers = 2 * Npairs + 1
    if Ncopies > 0:
        nLayers += (2 * Npairs) * Ncopies
    if Nfixed > 0:
        nLayers += 2 * Nfixed

    # Pre-build the refractive index array
    n_low = ifo.Materials.Coating.Indexlown
    n_high = ifo.Materials.Coating.Indexhighn
    n_sub = ifo.Materials.Substrate.RefractiveIndex
    doublet = np.tile(np.array([n_low, n_high]), nLayers // 2)
    if len(doublet) != nLayers:
        doublet = np.append(doublet, n_low)
    n = np.empty(len(doublet) + 2)
    n[0] = 1.0
    n[1:-1] = doublet
    n[-1] = n_sub
    misc['_ifo_n'] = n

    # Active costs (nonzero weight)
    active = frozenset(c for c, s in costs.items() if s['weight'])
    misc['_active_costs'] = active

    # AUX wavelength ratio
    lambdaAUX = misc.get('lambdaAUX', 1550 / 2050)
    misc['_lambdaAUX'] = lambdaAUX

    # Wavelength array and map for consolidated multidiel1 call
    wl_list, wl_map = [], {}
    if 'Trans1064' in active or 'Esurf' in active:
        wl_map['PSL'] = len(wl_list)
        wl_list.append(1.0)
    if 'Trans532' in active:
        wl_map['AUX'] = len(wl_list)
        wl_list.append(lambdaAUX)
    if 'TransOPLEV' in active:
        wl_map['OPL'] = len(wl_list)
        wl_list.append(0.297)
    misc['_wl_array'] = np.array(wl_list) if wl_list else None
    misc['_wl_map'] = wl_map

    # Sensitivity wavelengths
    if 'Lsens' in active:
        misc['_sens_wl'] = np.array([1.0, lambdaAUX])

    # Pre-extract ifo params for thermooptic JIT
    if 'Thermooptic' in active:
        misc['_ifo_params'] = extract_ifo_params(ifo)

    return misc


def getMirrorCost(L, costs, ifo, gam, verbose=False, misc={}):
    """Master cost function for coating optimization.

    Evaluates a multiplicative product of sub-costs:
    ``C = prod( 1 + w_i * c_i )`` where each ``c_i`` is O(1).
    Every factor is >= 1, so no single term can zero the product.
    This prevents the optimizer from sacrificing one objective for
    another, unlike an additive weighted sum.

    Uses consolidated ``multidiel1`` calls: one call with all active
    wavelengths for transmission/E-field costs, and one call with
    perturbed thicknesses for sensitivity cost.

    Parameters
    ----------
    L : array_like
        Optical thicknesses of the candidate dielectric stack.
    costs : dict
        Cost specifications. Each key maps to
        ``{'target': float, 'weight': float}``. Supported keys:
        ``Trans1064``, ``Trans532``, ``TransOPLEV``, ``Brownian``,
        ``Thermooptic``, ``Lsens``, ``Esurf``, ``Lstdev``, ``Absorption``.
    ifo : gwinc.Struct
        Interferometer model containing material properties.
    gam : float
        Pre-computed Brownian proxy factor from :func:`brownianProxy`.
    verbose : bool, optional
        If True, return ``(scalar_cost, output_dict)`` with detailed
        results. Default is False (return scalar only).
    misc : dict, optional
        Additional parameters: ``aoi``, ``pol``, ``fTO``, ``Npairs``,
        ``Ncopies``, ``Nfixed``, ``lambdaAUX``.
        Pre-computed cache keys (from :func:`precompute_misc`):
        ``_ifo_n``, ``_active_costs``, ``_wl_array``, ``_wl_map``.

    Returns
    -------
    scalar_cost : float
        Weighted sum of all active cost terms.
    output : dict
        Only returned when ``verbose=True``. Contains keys ``'n'``,
        ``'L'``, ``'scalarCost'``, ``'vectorCost'``, ``'T1064'``,
        ``'R1064'``, and additional transmission values if enabled.
    """

    # Add Ncopies of single variable stack
    if misc.get('Ncopies', 0) > 0:
        copiedLayers = np.tile(L[:-2].copy(), misc['Ncopies'])
        L = np.append(L, copiedLayers)

    # Add fixed layers at end? Default to 0
    if misc.get('Nfixed', 0) > 0:
        fixedLayers = np.tile(L[-2:].copy(), misc['Nfixed'])
        L = np.append(L, fixedLayers)

    # Use pre-computed n array if available, else build it
    n = misc.get('_ifo_n')
    if n is None:
        doublet = np.tile(np.array([ifo.Materials.Coating.Indexlown,
                                    ifo.Materials.Coating.Indexhighn]),
                              int(np.floor(len(L)/2)))
        if len(doublet) != len(L):
            doublet = np.append(doublet, doublet[0])
        n = np.append(1, doublet)
        n = np.append(n, ifo.Materials.Substrate.RefractiveIndex)

    # Use pre-computed active set or build it
    active = misc.get('_active_costs')
    if active is None:
        active = {c for c, s in costs.items() if s['weight']}

    # Use pre-computed lambdaAUX or read from misc
    lambdaAUX = misc.get('_lambdaAUX')
    if lambdaAUX is None:
        lambdaAUX = misc.get('lambdaAUX', 1550/2050)

    vector_cost, output = {}, {}
    aoi, pol = misc['aoi'], misc['pol']

    # --- Consolidated multidiel1 call for main wavelengths ---
    # Use pre-computed wavelength array/map if available
    wl_arr = misc.get('_wl_array')
    wl_map = misc.get('_wl_map')
    if wl_map is None:
        wl_list, wl_map = [], {}
        if 'Trans1064' in active or 'Esurf' in active:
            wl_map['PSL'] = len(wl_list)
            wl_list.append(1.0)
        if 'Trans532' in active:
            wl_map['AUX'] = len(wl_list)
            wl_list.append(lambdaAUX)
        if 'TransOPLEV' in active:
            wl_map['OPL'] = len(wl_list)
            wl_list.append(0.297)
        wl_arr = np.array(wl_list) if wl_list else None

    if wl_arr is not None and len(wl_arr) > 0:
        r_main, _ = multidiel1(n, L, wl_arr, aoi, pol)
        T_main = 1 - np.abs(r_main)**2

    # Extract per-wavelength results
    if 'Trans1064' in active:
        idx = wl_map['PSL']
        T1064 = T_main[idx]
        vector_cost['Trans1064'] = np.abs((costs['Trans1064']['target'] - T1064)
                                         / costs['Trans1064']['target'])**2
        if verbose:
            output['T1064'] = T1064
            output['R1064'] = 1 - T1064

    if 'Trans532' in active:
        idx = wl_map['AUX']
        T532 = T_main[idx]
        vector_cost['Trans532'] = np.abs((costs['Trans532']['target'] - T532)
                                         / costs['Trans532']['target'])**2
        if verbose:
            output['T532'] = T532
            output['R532'] = 1 - T532

    if 'TransOPLEV' in active:
        idx = wl_map['OPL']
        TOPL = T_main[idx]
        vector_cost['TransOPLEV'] = np.abs((costs['TransOPLEV']['target'] - TOPL)
                                           / costs['TransOPLEV']['target'])**2
        if verbose:
            output['TOPL'] = TOPL
            output['ROPL'] = 1 - TOPL

    if 'Esurf' in active:
        idx = wl_map['PSL']
        vector_cost['Esurf'] = 50 * np.arcsinh(np.abs(1 + r_main[idx])**2)

    # --- Brownian (no multidiel1 needed, inlined for speed) ---
    if 'Brownian' in active:
        vector_cost['Brownian'] = (
            L[::2].sum() + gam * L[1::2].sum()) / costs['Brownian']['target']

    # --- Thermooptic (JIT-compiled, no gwinc in hot path) ---
    if 'Thermooptic' in active:
        ifo_params = misc.get('_ifo_params')
        if ifo_params is None:
            ifo_params = extract_ifo_params(ifo)
        vector_cost['Thermooptic'] = thermoopticCost(
            costs['Thermooptic']['target'], misc['fTO'], L, ifo,
            ifo_params=ifo_params)

    # --- Absorption (E-field profile integration, slow — only runs if weight > 0) ---
    if 'Absorption' in active:
        wavelength = ifo.Laser.Wavelength
        L_phys = wavelength * op2phys(L, n[1:-1])
        nPts = 4  # coarse grid (matches OptimalBragg default)
        _, Esq = fieldDepth(L_phys, n, lam=wavelength, theta=aoi,
                            pol='s' if pol == 'te' else 'p', nPts=nPts)
        alpha_low = getattr(ifo.Optics.ETM.Coating, 'Absorptionlown', 0)
        alpha_high = getattr(ifo.Optics.ETM.Coating, 'Absorptionhighn', 0)
        absorp = calcAbsorption(Esq, L_phys, nPts, alpha_low, alpha_high)
        target_abs = costs['Absorption']['target']
        vector_cost['Absorption'] = np.abs((target_abs - absorp) / target_abs)

    # --- Layer thickness std dev (no multidiel1 needed) ---
    if 'Lstdev' in active:
        vector_cost['Lstdev'] = stdevLCost(costs['Lstdev']['target'], L)

    # --- Consolidated multidiel1 call for sensitivity (perturbed L) ---
    if 'Lsens' in active:
        L_pert = 1.01 * L
        sens_wl = misc.get('_sens_wl')
        if sens_wl is None:
            sens_wl = np.array([1.0, lambdaAUX])
        r_pert, _ = multidiel1(n, L_pert, sens_wl, aoi, pol)
        T_pert = 1 - np.abs(r_pert)**2
        target_sens = costs['Lsens']['target']
        sensPSL = np.abs((target_sens - T_pert[0]) / target_sens)**2
        sensAUX = np.abs((target_sens - T_pert[1]) / target_sens)**2
        vector_cost['Lsens'] = np.sqrt(sensPSL**2 + sensAUX**2)

    # Multiplicative product: C = prod( 1 + w*c_i )
    # Each factor >= 1, so no single zero can collapse the product.
    scalar_cost = 1.0
    for cost_name, cost_val in vector_cost.items():
        w = costs[cost_name]['weight']
        scalar_cost *= (1 + w * cost_val)

    if verbose:
        for cost_name in vector_cost:
            print(cost_name + f' cost = {vector_cost[cost_name]:.4f}')
        output['n'] = n
        output['L'] = L
        output['scalarCost'] = scalar_cost
        output['vectorCost'] = vector_cost
        return scalar_cost, output
    else:
        return scalar_cost
