# Set of functions used to evaluate a cost function for optimization

import copy
import numpy as np
from generic.coatingUtils import *
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

    Uses the formula from LIGO-E0900068 p. 4:
    ``cost = target * (sum_low + gamma * sum_high)``.

    Parameters
    ----------
    target : float
        Scaling factor for the Brownian noise cost.
    L : array_like
        Optical thicknesses of the dielectric stack.
    gam : float
        Pre-computed Brownian proxy factor from :func:`brownianProxy`.

    Returns
    -------
    cost : float
        Brownian noise proxy cost.
    """
    zLow = np.sum(L[::2])    # Sum of thicknesses of low index layers
    zHigh = np.sum(L[1::2])   # Sum of thicknesses of high index layers
    SBrZ = zLow + gam*zHigh  # Proxy brownian noise
    return target * SBrZ


def thermoopticCost(target, fTarget, L, ifo):
    """Compute thermo-optic noise cost via pygwinc.

    Evaluates the thermo-optic noise power spectral density at a single
    frequency using ``gwinc.noise.coatingthermal.coating_thermooptic``.

    Parameters
    ----------
    target : float
        Scaling factor for the thermo-optic cost.
    fTarget : float
        Frequency at which to evaluate TO noise, in Hz.
    L : array_like
        Optical thicknesses of the dielectric stack.
    ifo : gwinc.Struct
        Interferometer model containing material and optic properties.

    Returns
    -------
    cost : float
        Thermo-optic noise cost: ``target * S_TO(fTarget)``.
    """
    # Get the TO noise PSD
    # Build up a "mirror" structure as required by pygwinc
    # Use shallow copies to avoid mutating the shared ifo object
    mir = copy.copy(ifo.Optics.ETM)
    mir.Coating = copy.copy(ifo.Optics.ETM.Coating)
    mir.Coating.dOpt = L
    StoZ, _, _, _ = gwinc.noise.coatingthermal.coating_thermooptic(
        fTarget, mir, ifo.Laser.Wavelength, ifo.Optics.ETM.BeamRadius)
    return target * StoZ


def getMirrorCost(L, costs, ifo, gam, verbose=False, misc={}):
    """Master cost function for coating optimization.

    Evaluates a weighted sum of sub-costs (transmission, Brownian noise,
    thermo-optic noise, sensitivity, surface E-field, thickness
    uniformity) for a candidate coating design. Called by
    ``scipy.optimize.differential_evolution`` on every candidate.

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
        ``TransPSL``, ``TransAUX``, ``TransOPLEV``, ``Brownian``,
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

    Returns
    -------
    scalar_cost : float
        Weighted sum of all active cost terms.
    output : dict
        Only returned when ``verbose=True``. Contains keys ``'n'``,
        ``'L'``, ``'scalarCost'``, ``'vectorCost'``, ``'TPSL'``,
        ``'RPSL'``, and additional transmission values if enabled.
    """

    # Add Ncopies of single variable stack
    if misc.get('Ncopies', 0) > 0:
        copiedLayers = np.tile(L[:-2].copy(), misc['Ncopies'])
        L = np.append(L, copiedLayers)

    # Add fixed layers at end? Default to 0
    if misc.get('Nfixed', 0) > 0:
        fixedLayers = np.tile(L[-2:].copy(), misc['Nfixed'])
        L = np.append(L, fixedLayers)

    # Build up the array of refractive indices
    doublet = np.tile(np.array([ifo.Materials.Coating.Indexlown,
                                ifo.Materials.Coating.Indexhighn]),
                          int(np.floor(len(L)/2)))

    if len(doublet) != len(L):
        # Add another low index layer at the bottom of the stack
        doublet = np.append(doublet, doublet[0])

    # Add air, fixed extension, and substrate
    n = np.append(1, doublet)
    n = np.append(n, ifo.Materials.Substrate.RefractiveIndex)

    # AUX wavelength ratio (configurable per project)
    lambdaAUX = misc.get('lambdaAUX', 1550/2050)

    # Determine which costs are active (nonzero weight)
    active = {c for c, s in costs.items() if s['weight']}

    vector_cost, output = {}, {}
    scalar_cost = 0.0
    aoi, pol = misc['aoi'], misc['pol']

    # --- Consolidated multidiel1 call for main wavelengths ---
    # Build wavelength array and index map for all transmission/Esurf costs
    wl_list, wl_map = [], {}
    if 'TransPSL' in active or 'Esurf' in active:
        wl_map['PSL'] = len(wl_list)
        wl_list.append(1.0)
    if 'TransAUX' in active:
        wl_map['AUX'] = len(wl_list)
        wl_list.append(lambdaAUX)
    if 'TransOPLEV' in active:
        wl_map['OPL'] = len(wl_list)
        wl_list.append(0.297)

    if wl_list:
        r_main, _ = multidiel1(n, L, np.array(wl_list), aoi, pol)
        T_main = 1 - np.abs(r_main)**2

    # Extract per-wavelength results
    if 'TransPSL' in active:
        idx = wl_map['PSL']
        TPSL = T_main[idx]
        vector_cost['TransPSL'] = np.abs((costs['TransPSL']['target'] - TPSL)
                                         / costs['TransPSL']['target'])**2
        if verbose:
            output['TPSL'] = TPSL
            output['RPSL'] = 1 - TPSL

    if 'TransAUX' in active:
        idx = wl_map['AUX']
        TAUX = T_main[idx]
        vector_cost['TransAUX'] = np.abs((costs['TransAUX']['target'] - TAUX)
                                         / costs['TransAUX']['target'])**2
        if verbose:
            output['TAUX'] = TAUX
            output['RAUX'] = 1 - TAUX

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

    # --- Brownian (no multidiel1 needed) ---
    if 'Brownian' in active:
        vector_cost['Brownian'] = brownianCost(costs['Brownian']['target'], L, gam)

    # --- Thermooptic (uses gwinc, not multidiel1) ---
    if 'Thermooptic' in active:
        vector_cost['Thermooptic'] = thermoopticCost(
            costs['Thermooptic']['target'], misc['fTO'], L, ifo)

    # --- Absorption (not implemented) ---

    # --- Layer thickness std dev (no multidiel1 needed) ---
    if 'Lstdev' in active:
        vector_cost['Lstdev'] = stdevLCost(costs['Lstdev']['target'], L)

    # --- Consolidated multidiel1 call for sensitivity (perturbed L) ---
    if 'Lsens' in active:
        L_pert = 1.01 * L
        sens_wl = [1.0, lambdaAUX]
        r_pert, _ = multidiel1(n, L_pert, np.array(sens_wl), aoi, pol)
        T_pert = 1 - np.abs(r_pert)**2
        target_sens = costs['Lsens']['target']
        sensPSL = np.abs((target_sens - T_pert[0]) / target_sens)**2
        sensAUX = np.abs((target_sens - T_pert[1]) / target_sens)**2
        vector_cost['Lsens'] = np.sqrt(sensPSL**2 + sensAUX**2)

    # Weighted sum
    for cost_name, cost_val in vector_cost.items():
        scalar_cost += costs[cost_name]['weight'] * cost_val

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
