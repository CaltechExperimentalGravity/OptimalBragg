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
    wl_list, wl_map = [], {}
    if 'Trans1' in active or 'Esurf' in active:
        wl_map['PSL'] = len(wl_list)
        wl_list.append(1.0)
    if 'Trans2' in active:
        wl_map['AUX'] = len(wl_list)
        wl_list.append(lambda2)
    if 'Trans3' in active:
        wl_map['OPL'] = len(wl_list)
        wl_list.append(lambda3)
    misc['_wl_array'] = np.array(wl_list) if wl_list else None
    misc['_wl_map'] = wl_map

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
        wl_list, wl_map = [], {}
        if 'Trans1' in active or 'Esurf' in active:
            wl_map['PSL'] = len(wl_list)
            wl_list.append(1.0)
        if 'Trans2' in active:
            wl_map['AUX'] = len(wl_list)
            wl_list.append(lambda2)
        if 'Trans3' in active:
            wl_map['OPL'] = len(wl_list)
            wl_list.append(lambda3)
        wl_arr = np.array(wl_list) if wl_list else None

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
        T2 = T_main[idx]
        vector_cost['Trans2'] = np.abs(
            (costs['Trans2']['target'] - T2)
            / costs['Trans2']['target']
        ) ** 2
        if verbose:
            output['T2'] = T2
            output['R2'] = 1 - T2

    if 'Trans3' in active:
        idx = wl_map['OPL']
        T3 = T_main[idx]
        vector_cost['Trans3'] = np.abs(
            (costs['Trans3']['target'] - T3)
            / costs['Trans3']['target']
        ) ** 2
        if verbose:
            output['T3'] = T3

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

    # --- Sensitivity: per-layer curvature (d²T/dL²) penalty ---
    # Fragile designs sit at narrow T peaks where dT/dL ≈ 0 but
    # d²T/dL² is large and negative — small perturbations cause T
    # to fall off a cliff.  We penalize the mean squared curvature
    # across all layers and active transmission targets.
    if 'Sensitivity' in active:
        eps = 5e-3  # perturbation in optical thickness units
        # Collect active transmission targets and their nominal T values
        sens_targets = []  # (cost_name, wavelength_ratio)
        sens_T_nom = []    # nominal T at each wavelength
        if 'Trans1' in active:
            sens_targets.append(('Trans1', 1.0))
            sens_T_nom.append(T1)
        if 'Trans2' in active:
            sens_targets.append(('Trans2', lambda2))
            sens_T_nom.append(T2)
        if 'Trans3' in active:
            sens_targets.append(('Trans3', lambda3))
            sens_T_nom.append(T3)

        if sens_targets:
            wl_sens = np.array([s[1] for s in sens_targets])
            T_nom_arr = np.array(sens_T_nom)
            N_layers = len(L)
            curv_sq_sum = 0.0

            for j in range(N_layers):
                Lp = L.copy()
                Lp[j] += eps
                Lm = L.copy()
                Lm[j] = max(L[j] - eps, 0.01)
                eps_m = L[j] - Lm[j]
                rp, _ = multidiel1(n, Lp, wl_sens, aoi, pol)
                rm, _ = multidiel1(n, Lm, wl_sens, aoi, pol)
                Tp = 1 - np.abs(rp) ** 2
                Tm = 1 - np.abs(rm) ** 2
                d2TdL2 = (Tp + Tm - 2 * T_nom_arr) / (eps * eps_m)
                curv_sq_sum += np.sum(d2TdL2 ** 2)

            n_terms = N_layers * len(sens_targets)
            mean_curv_sq = curv_sq_sum / n_terms
            # Threshold penalty: zero at or below target, quadratic above.
            # Target is the curvature budget (mean squared d²T/dL²).
            target_curv = costs['Sensitivity']['target']
            excess = mean_curv_sq / target_curv - 1.0
            vector_cost['Sensitivity'] = max(0.0, excess) ** 2
        else:
            vector_cost['Sensitivity'] = 0.0

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
