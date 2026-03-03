"""Cost functions for coating optimization.

Uses a multiplicative scalar cost ``C = prod(1 + w_i * c_i)`` where
each factor ``c_i`` is O(1).  Every factor is >= 1, so no single term
can zero the product — the optimizer cannot sacrifice one objective
for another, unlike an additive weighted sum.

Two APIs:

- **Individual cost functions** (``transmissionCost``, ``brownianCost``,
  etc.) — each returns a single scalar cost, useful for analysis.
- **Master evaluator** :func:`getMirrorCost` — the hot-path function
  called by ``differential_evolution``.  Uses consolidated ``multidiel1``
  calls and pre-computed caches for speed.

All functions work with a stack dict (from :func:`OptimalBragg.qw_stack`)
and optical thicknesses.  No gwinc dependency.
"""

import numpy as np

from OptimalBragg.layers import multidiel1, field_zmag, calc_abs, op2phys
from OptimalBragg.noise import (
    coating_thermooptic_fast,
    extract_stack_params,
    brownian_proxy,
)


# ── Normalization decorator ──────────────────────────────────────────

def norm(norm_arg):
    """Decorator to apply a normalization to a cost function's return value.

    Supported norms: ``"l1"`` (identity), ``"l2"`` (square),
    ``"arcsinh"`` (smooth compression of large values).
    """
    def decorator(func):
        def wrapper(*args, **kwargs):
            ci = func(*args, **kwargs)
            if norm_arg == "l2":
                return ci ** 2
            elif norm_arg == "arcsinh":
                return np.arcsinh(ci)
            return ci  # l1
        wrapper.__name__ = func.__name__
        wrapper.__doc__ = func.__doc__
        return wrapper
    return decorator


# ── Individual cost functions ────────────────────────────────────────

def transmissionCost(target, n, L, lamb=1, theta=0, pol='te'):
    """Transmission cost: ``|target - T|^2 / target^2``.

    Parameters
    ----------
    target : float
        Target power transmission.
    n : array_like
        Refractive indices [incident, layers..., substrate].
    L : array_like
        Optical thicknesses (fraction of design wavelength).
    lamb : float or array_like
        Normalized wavelength(s). Default 1.
    theta : float
        Angle of incidence [degrees]. Default 0.
    pol : str
        Polarization. Default ``'te'``.

    Returns
    -------
    cost : float
        Normalized squared error.
    T : float or ndarray
        Power transmission.
    """
    r, _ = multidiel1(n, L, lamb, theta, pol)
    T = 1 - np.abs(r) ** 2
    return (np.abs((target - T) / target) ** 2)[0], T


def sensitivityCost(target, n, L, lamb=1, theta=0, pol='te'):
    """Layer thickness sensitivity: transmission cost at 1.01 * L."""
    return transmissionCost(target, n, 1.01 * L, lamb, theta, pol)[0]


def surfEfieldCost(target, n, L, lamb=1, theta=0, pol='te'):
    """Surface E-field threshold penalty.

    Returns 0 when ``|1+r|^2 <= target``, quadratic penalty above.
    For high reflectivity coatings with HL pattern (high-n on top),
    ``|1+r|^2`` is naturally near zero.
    """
    r, _ = multidiel1(n, L, lamb, theta, pol)
    Esurf_sq = np.abs(1 + r[0]) ** 2
    return max(0.0, Esurf_sq / target - 1.0) ** 2


def stdevLCost(target, L):
    """Layer thickness uniformity: penalizes mean/std deviation from target."""
    if np.std(np.array(L)):
        relative_stdev = np.mean(np.array(L)) / np.std(np.array(L))
    else:
        relative_stdev = 0.0
    return np.abs((target - relative_stdev) / target) ** 2


def brownianCost(target, L, gam):
    """Brownian noise threshold penalty (LIGO-E0900068).

    Returns 0 when noise proxy is at or below *target*, and a quadratic
    penalty ``(SBrZ/target - 1)^2`` when above.

    Parameters
    ----------
    target : float
        Brownian proxy budget (e.g. QW stack reference).
    L : array_like
        Optical thicknesses.
    gam : float or dict
        Brownian proxy factor.  If dict, uses the first value
        (for binary stacks with a single high-n material).

    Returns
    -------
    float
        Threshold penalty: ``max(0, SBrZ/target - 1)^2``.
    """
    if isinstance(gam, dict):
        gam_val = next(iter(gam.values()))
    else:
        gam_val = gam
    zLow = np.sum(L[::2])
    zHigh = np.sum(L[1::2])
    SBrZ = zLow + gam_val * zHigh
    return max(0.0, SBrZ / target - 1.0) ** 2


def thermoopticCost(target, fTarget, L, stack, stack_params=None,
                    w_beam=None):
    """Thermo-optic noise threshold penalty (JIT-compiled).

    Returns 0 when S_TO is at or below *target*, quadratic above.

    Parameters
    ----------
    target : float
        S_TO budget for normalization.
    fTarget : float
        Frequency [Hz].
    L : array_like
        Optical thicknesses.
    stack : dict
        Stack dict.
    stack_params : tuple, optional
        Pre-computed from :func:`extract_stack_params`.
    w_beam : float, optional
        Beam radius [m].  Required if *stack_params* is None.

    Returns
    -------
    float
        Normalized TO noise: ``S_TO / target``.
    """
    wavelength = stack["lam_ref"]
    if stack_params is None:
        raise ValueError("stack_params required (from extract_stack_params)")
    if w_beam is None:
        raise ValueError("w_beam required")
    StoZ = coating_thermooptic_fast(fTarget, L, wavelength, w_beam,
                                    stack_params)
    return max(0.0, StoZ / target - 1.0) ** 2


# ── Pre-computation cache ────────────────────────────────────────────

def _build_wavelength_map(active, costs, lambda2, lambda3):
    """Build consolidated wavelength array and index map for multidiel1.

    Returns ``(wl_array, wl_map)`` where *wl_array* is a 1-D float array
    (or ``None`` if no wavelengths are needed) and *wl_map* maps label
    strings to indices or slices into *wl_array*.
    """
    wl_list, wl_map = [], {}
    if 'Trans1' in active or 'Esurf' in active:
        wl_map['PSL'] = len(wl_list)
        wl_list.append(1.0)
    if 'Trans2' in active:
        bw = costs['Trans2'].get('bandwidth', 0)
        if bw > 0:
            n_bw = 7
            lams = np.linspace(lambda2 * (1 - bw), lambda2 * (1 + bw), n_bw)
            wl_map['AUX'] = slice(len(wl_list), len(wl_list) + n_bw)
            wl_map['AUX_center'] = len(wl_list) + n_bw // 2
            wl_list.extend(lams.tolist())
        else:
            wl_map['AUX'] = len(wl_list)
            wl_list.append(lambda2)
    if 'Trans3' in active:
        bw = costs['Trans3'].get('bandwidth', 0)
        if bw > 0:
            n_bw = 7
            lams = np.linspace(lambda3 * (1 - bw), lambda3 * (1 + bw), n_bw)
            wl_map['OPL'] = slice(len(wl_list), len(wl_list) + n_bw)
            wl_map['OPL_center'] = len(wl_list) + n_bw // 2
            wl_list.extend(lams.tolist())
        else:
            wl_map['OPL'] = len(wl_list)
            wl_list.append(lambda3)
    wl_array = np.array(wl_list) if wl_list else None
    return wl_array, wl_map


def precompute_misc(costs, stack, misc):
    """Pre-compute cached values for the getMirrorCost hot loop.

    Call once before ``differential_evolution`` to populate *misc*
    with values that are constant across all evaluations.

    Parameters
    ----------
    costs : dict
        Cost specifications (same as for :func:`getMirrorCost`).
    stack : dict
        Stack dict.
    misc : dict
        Updated in-place with cache keys.

    Returns
    -------
    dict
        The updated *misc* dict.
    """
    Ncopies = misc.get('Ncopies', 0)
    Nfixed = misc.get('Nfixed', 0)

    # Build refractive index array from stack, extending for copies/fixed
    n_base = stack["ns"]  # [superstrate, layer1, ..., layerN, substrate]
    n_layers_only = n_base[1:-1]  # just the coating layers

    if Ncopies > 0 or Nfixed > 0:
        n_low = stack["thin_films"]["L"].Index
        n_high = stack["thin_films"]["H"].Index
        extra = []
        if Ncopies > 0:
            extra.append(np.tile(n_layers_only[:-2], Ncopies))
        if Nfixed > 0:
            extra.append(np.tile(n_layers_only[-2:], Nfixed))
        n_layers_ext = np.concatenate([n_layers_only] + extra)
        n = np.empty(len(n_layers_ext) + 2)
        n[0] = n_base[0]
        n[1:-1] = n_layers_ext
        n[-1] = n_base[-1]
    else:
        n = n_base

    misc['_ifo_n'] = np.asarray(n, dtype=float)

    # Active costs (nonzero weight)
    active = frozenset(c for c, s in costs.items() if s['weight'])
    misc['_active_costs'] = active

    # Wavelength ratios (generic names)
    lambda2 = misc.get('lambda2', 1550 / 2050)
    misc['_lambda2'] = lambda2

    lambda3 = misc.get('lambda3', 0.297)
    misc['_lambda3'] = lambda3

    # Wavelength array and map for consolidated multidiel1 call
    misc['_wl_array'], misc['_wl_map'] = _build_wavelength_map(
        active, costs, lambda2, lambda3
    )

    # Sensitivity wavelengths
    if 'Lsens' in active:
        misc['_sens_wl'] = np.array([1.0, lambda2])

    # Pre-extract stack params for thermooptic JIT
    if 'Thermooptic' in active:
        r_mirror = misc.get('r_mirror', 0.17)
        d_mirror = misc.get('d_mirror', 0.20)
        misc['_stack_params'] = extract_stack_params(
            stack, r_mirror, d_mirror
        )

    return misc


# ── Master cost evaluator ────────────────────────────────────────────

def getMirrorCost(L, costs, stack, gam, verbose=False, misc=None):
    """Master cost function for coating optimization.

    Evaluates a multiplicative product of sub-costs:
    ``C = prod(1 + w_i * c_i)`` where each ``c_i`` is O(1).

    Uses consolidated ``multidiel1`` calls: one call with all active
    wavelengths for transmission/E-field costs, and one call with
    perturbed thicknesses for sensitivity cost.

    Parameters
    ----------
    L : array_like
        Optical thicknesses of the candidate dielectric stack.
    costs : dict
        Cost specifications. Each key maps to
        ``{'target': float, 'weight': float}``.  Supported keys:
        ``Trans1``, ``Trans2``, ``Trans3``, ``Brownian``,
        ``Thermooptic``, ``Lsens``, ``Sensitivity``, ``Esurf``,
        ``Lstdev``, ``Absorption``.
    stack : dict
        Stack dict from :func:`OptimalBragg.qw_stack`.
    gam : float or dict
        Brownian noise proxy factor from :func:`brownian_proxy`.
    verbose : bool, optional
        If True, return ``(scalar_cost, output_dict)``.
    misc : dict, optional
        Additional parameters and pre-computed cache.

    Returns
    -------
    scalar_cost : float
        Multiplicative product of all active cost terms.
    output : dict
        Only returned when ``verbose=True``.
    """
    if misc is None:
        misc = {}

    # Extend L with copies/fixed layers
    if misc.get('Ncopies', 0) > 0:
        copiedLayers = np.tile(L[:-2].copy(), misc['Ncopies'])
        L = np.append(L, copiedLayers)

    if misc.get('Nfixed', 0) > 0:
        fixedLayers = np.tile(L[-2:].copy(), misc['Nfixed'])
        L = np.append(L, fixedLayers)

    # Use pre-computed n array or build it
    n = misc.get('_ifo_n')
    if n is None:
        n_low = stack["thin_films"]["L"].Index
        n_high = stack["thin_films"]["H"].Index
        n_sub = stack["sub"].Index
        doublet = np.tile(np.array([n_low, n_high]),
                          int(np.floor(len(L) / 2)))
        if len(doublet) != len(L):
            doublet = np.append(doublet, doublet[0])
        n = np.append(1, doublet)
        n = np.append(n, n_sub)

    active = misc.get('_active_costs')
    if active is None:
        active = {c for c, s in costs.items() if s['weight']}

    lambda2 = misc.get('_lambda2')
    if lambda2 is None:
        lambda2 = misc.get('lambda2', 1550 / 2050)

    lambda3 = misc.get('_lambda3')
    if lambda3 is None:
        lambda3 = misc.get('lambda3', 0.297)

    vector_cost, output = {}, {}
    aoi, pol = misc.get('aoi', 0), misc.get('pol', 'te')

    # --- Consolidated multidiel1 call for main wavelengths ---
    wl_arr = misc.get('_wl_array')
    wl_map = misc.get('_wl_map')
    if wl_map is None:
        wl_arr, wl_map = _build_wavelength_map(
            active, costs, lambda2, lambda3
        )

    if wl_arr is not None and len(wl_arr) > 0:
        r_main, _ = multidiel1(n, L, wl_arr, aoi, pol)
        T_main = 1 - np.abs(r_main) ** 2

    # Extract per-wavelength results
    if 'Trans1' in active:
        idx = wl_map['PSL']
        T1 = T_main[idx]
        vector_cost['Trans1'] = np.abs(
            (costs['Trans1']['target'] - T1)
            / costs['Trans1']['target']
        ) ** 2
        if verbose:
            output['T1'] = T1
            output['R1'] = 1 - T1

    if 'Trans2' in active:
        idx = wl_map['AUX']
        T2_band = T_main[idx]
        cost_band = np.abs(
            (costs['Trans2']['target'] - T2_band)
            / costs['Trans2']['target']
        ) ** 2
        vector_cost['Trans2'] = float(np.max(cost_band))
        if verbose:
            if isinstance(idx, slice):
                output['T2'] = float(T_main[wl_map['AUX_center']])
            else:
                output['T2'] = float(T2_band)
            output['R2'] = 1 - output['T2']

    if 'Trans3' in active:
        idx = wl_map['OPL']
        T3_band = T_main[idx]
        cost_band = np.abs(
            (costs['Trans3']['target'] - T3_band)
            / costs['Trans3']['target']
        ) ** 2
        vector_cost['Trans3'] = float(np.max(cost_band))
        if verbose:
            if isinstance(idx, slice):
                output['T3'] = float(T_main[wl_map['OPL_center']])
            else:
                output['T3'] = float(T3_band)
            output['R3'] = 1 - output['T3']

    if 'Esurf' in active:
        idx = wl_map['PSL']
        Esurf_sq = np.abs(1 + r_main[idx]) ** 2
        target_Esurf = costs['Esurf']['target']
        excess = Esurf_sq / target_Esurf - 1.0
        vector_cost['Esurf'] = max(0.0, excess) ** 2
        if verbose:
            output['Esurf'] = Esurf_sq

    # --- Brownian (no multidiel1 needed) ---
    # Threshold penalty: zero at/below target, quadratic above.
    if 'Brownian' in active:
        gam_val = next(iter(gam.values())) if isinstance(gam, dict) else gam
        SBrZ = (L[::2].sum() + gam_val * L[1::2].sum())
        excess = SBrZ / costs['Brownian']['target'] - 1.0
        vector_cost['Brownian'] = max(0.0, excess) ** 2

    # --- Thermooptic (JIT-compiled) ---
    # Threshold penalty: zero at/below target, quadratic above.
    if 'Thermooptic' in active:
        stack_params = misc.get('_stack_params')
        if stack_params is None:
            r_mirror = misc.get('r_mirror', 0.17)
            d_mirror = misc.get('d_mirror', 0.20)
            stack_params = extract_stack_params(stack, r_mirror, d_mirror)
        w_beam = misc.get('w_beam', 0.062)
        StoZ = coating_thermooptic_fast(
            misc.get('fTO', 100.0), L, stack["lam_ref"], w_beam,
            stack_params,
        )
        excess = StoZ / costs['Thermooptic']['target'] - 1.0
        vector_cost['Thermooptic'] = max(0.0, excess) ** 2

    # --- Absorption (E-field profile, slow) ---
    if 'Absorption' in active:
        wavelength = stack["lam_ref"]
        L_phys = wavelength * op2phys(L, n[1:-1])
        nPts = 4
        _, Esq = field_zmag(
            n, L_phys, lam=wavelength,
            aoi=np.radians(aoi) if aoi else 0,
            pol='s' if pol == 'te' else 'p', n_pts=nPts,
        )
        alphas = stack["alphas"]
        absorp = calc_abs(Esq, L_phys, alphas, n_pts=nPts)
        target_abs = costs['Absorption']['target']
        vector_cost['Absorption'] = np.abs(
            (target_abs - absorp) / target_abs
        )

    # --- Layer thickness std dev ---
    if 'Lstdev' in active:
        vector_cost['Lstdev'] = stdevLCost(costs['Lstdev']['target'], L)

    # --- Sensitivity (perturbed L) ---
    if 'Lsens' in active:
        L_pert = 1.01 * L
        sens_wl = misc.get('_sens_wl')
        if sens_wl is None:
            sens_wl = np.array([1.0, lambda2])
        r_pert, _ = multidiel1(n, L_pert, sens_wl, aoi, pol)
        T_pert = 1 - np.abs(r_pert) ** 2
        target_sens = costs['Lsens']['target']
        sensPSL = np.abs((target_sens - T_pert[0]) / target_sens) ** 2
        sensAUX = np.abs((target_sens - T_pert[1]) / target_sens) ** 2
        vector_cost['Lsens'] = np.sqrt(sensPSL ** 2 + sensAUX ** 2)

    # --- First-derivative penalties (3 independent multiplicative terms) ---
    # Each manufacturing DOF gets its own term in the product:
    #   dTdnH  — dT/dn_H  (high-index material error)
    #   dTdnL  — dT/dn_L  (low-index material error)
    #   dTdd   — dT/d_thickness (thickness calibration error)
    # Computed via central finite differences.  log1p compression.
    deriv_active = {'dTdnH', 'dTdnL', 'dTdd'} & active
    if deriv_active:
        # Shared delta from any active derivative cost
        for dname in ('dTdnH', 'dTdnL', 'dTdd'):
            if dname in active:
                delta = costs[dname].get('delta', 1e-2)
                break

        # Collect active wavelengths for derivative evaluation.
        # For broadband costs (bandwidth > 0), use ALL band wavelengths
        # so derivatives capture curvature across the full band.
        sens_wls = []
        if 'Trans1' in active and costs['Trans1']['target'] >= 0.5:
            sens_wls.append(1.0)
        if 'Trans2' in active:
            idx2 = wl_map.get('AUX')
            if isinstance(idx2, slice) and wl_arr is not None:
                sens_wls.extend(wl_arr[idx2].tolist())
            else:
                sens_wls.append(lambda2)
        if 'Trans3' in active:
            idx3 = wl_map.get('OPL')
            if isinstance(idx3, slice) and wl_arr is not None:
                sens_wls.extend(wl_arr[idx3].tolist())
            else:
                sens_wls.append(lambda3)

        if sens_wls:
            wl_sens = np.array(sens_wls)

            # Build high-n / low-n masks
            n_layers_only = n[1:-1]
            n_thresh = np.sqrt(n_layers_only.min() * n_layers_only.max())
            low_mask = np.zeros(len(n), dtype=bool)
            high_mask = np.zeros(len(n), dtype=bool)
            for ii in range(1, len(n) - 1):
                if n[ii] < n_thresh:
                    low_mask[ii] = True
                else:
                    high_mask[ii] = True

            if 'dTdnH' in active:
                n_hup = n.copy(); n_hup[high_mask] *= (1 + delta)
                n_hdn = n.copy(); n_hdn[high_mask] *= (1 - delta)
                r_hup, _ = multidiel1(n_hup, L, wl_sens, aoi, pol)
                r_hdn, _ = multidiel1(n_hdn, L, wl_sens, aoi, pol)
                dT = ((1 - np.abs(r_hup)**2) - (1 - np.abs(r_hdn)**2)) / (2 * delta)
                vector_cost['dTdnH'] = np.mean(np.log1p(np.abs(dT)))

            if 'dTdnL' in active:
                n_lup = n.copy(); n_lup[low_mask] *= (1 + delta)
                n_ldn = n.copy(); n_ldn[low_mask] *= (1 - delta)
                r_lup, _ = multidiel1(n_lup, L, wl_sens, aoi, pol)
                r_ldn, _ = multidiel1(n_ldn, L, wl_sens, aoi, pol)
                dT = ((1 - np.abs(r_lup)**2) - (1 - np.abs(r_ldn)**2)) / (2 * delta)
                vector_cost['dTdnL'] = np.mean(np.log1p(np.abs(dT)))

            if 'dTdd' in active:
                r_tup, _ = multidiel1(n, L * (1 + delta), wl_sens, aoi, pol)
                r_tdn, _ = multidiel1(n, L * (1 - delta), wl_sens, aoi, pol)
                dT = ((1 - np.abs(r_tup)**2) - (1 - np.abs(r_tdn)**2)) / (2 * delta)
                vector_cost['dTdd'] = np.mean(np.log1p(np.abs(dT)))
        else:
            for dname in deriv_active:
                vector_cost[dname] = 0.0

    # Multiplicative product: C = prod(1 + w * c_i)
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
