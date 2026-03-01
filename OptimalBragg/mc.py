"""Monte Carlo sensitivity analysis for coating designs.

Perturbs refractive indices and layer thicknesses to evaluate
sensitivity of transmission, thermal noise, and surface E-field
to manufacturing tolerances.

Perturbation parameters (3D):
  - High-n index (0.5% Gaussian)
  - Low-n index (0.5% Gaussian)
  - Layer thickness (0.5% Gaussian)

Usage::

    from OptimalBragg.mc import run_mc
    result = run_mc('Data/ETM/ETM_Layers_260221.hdf5', n_samples=5000)

Or from command line::

    python -m OptimalBragg.mc Data/ETM/ETM_Layers.hdf5 output_MC.hdf5 5000
"""

import numpy as np
import emcee
import h5py
import tqdm

from OptimalBragg.layers import multidiel1, op2phys


def _extract_wavelength_ratios(layers_hdf5):
    """Read lambda2/lambda3 from the params YAML stored in the HDF5.

    Returns
    -------
    lambda2 : float
        Second wavelength ratio (default 0.5 if not found).
    lambda3 : float or None
        Third wavelength ratio (None if not present).
    """
    from OptimalBragg.io import yamlread
    from pathlib import Path

    with h5py.File(layers_hdf5, 'r') as f:
        if 'params_file' not in f:
            return 0.5, None
        raw = np.array(f['params_file'])
        pf_str = raw.item().decode('utf-8') if isinstance(raw.item(), bytes) else str(raw.item())

    params_path = Path(pf_str)
    if not params_path.exists():
        return 0.5, None

    params = yamlread(str(params_path))
    misc = params.get('misc', {})
    return misc.get('lambda2', 0.5), misc.get('lambda3', None)


def _generate_perturbations(n_samples, n_dim=3, n_walkers=128,
                            width=0.005):
    """Generate Gaussian perturbation samples using emcee.

    Parameters
    ----------
    n_samples : int
        Number of perturbation samples to generate.
    n_dim : int
        Number of perturbation dimensions (3: high-n, low-n, thickness).
    n_walkers : int
        Number of emcee walkers.
    width : float
        Standard deviation of Gaussian perturbations (default 0.5%).

    Returns
    -------
    perturbs : ndarray
        Shape ``(n_samples, n_dim)`` — fractional perturbations
        centered at zero.
    """
    means = np.zeros(n_dim)
    cov = np.diag(width * np.ones(n_dim))
    cov = np.dot(cov, cov)
    icov = np.linalg.inv(cov)

    def lnprob(x, mu, icov):
        diff = x - mu
        return -0.5 * np.dot(diff, np.dot(icov, diff))

    p0 = np.random.rand(n_dim * n_walkers).reshape((n_walkers, n_dim))
    sampler = emcee.EnsembleSampler(n_walkers, n_dim, lnprob,
                                    args=[means, icov])

    # Burn-in
    pos, _, _ = sampler.run_mcmc(p0, 2000, progress=False)
    sampler.reset()

    # Production
    n_steps = max(n_samples // n_walkers + 1, 200)
    sampler.run_mcmc(pos, n_steps, progress=False)

    return sampler.get_chain(flat=True)[:n_samples, :]


def run_mc(layers_hdf5, n_samples=5000, lambda2=None, lambda3=None,
           w_beam=0.062, r_mirror=0.17, d_mirror=0.20,
           freq=None):
    """Run Monte Carlo sensitivity analysis on an optimized coating.

    Parameters
    ----------
    layers_hdf5 : str
        Path to optimizer output HDF5 (must contain ``diffevo_output/n``
        and ``diffevo_output/L``).
    n_samples : int
        Number of MC samples.
    lambda2 : float, optional
        Second wavelength ratio (lambda_2 / lambda_design).
        Auto-detected from params YAML if None.
    lambda3 : float, optional
        Third wavelength ratio. Auto-detected from params YAML if None.
    w_beam : float
        Beam radius [m].
    r_mirror, d_mirror : float
        Mirror radius and thickness [m].
    freq : array_like, optional
        Frequency array [Hz]. Default ``logspace(1, 3, 50)``.

    Returns
    -------
    result : dict
        Keys: ``MCout`` (5 or 6 x n_samples), ``TOnoise``, ``Brnoise``,
        ``freq``, ``perturbs``, ``has_lambda3``.
    """
    from OptimalBragg.noise import coating_thermooptic, coating_brownian

    if freq is None:
        freq = np.logspace(1, 3, 50)

    # Auto-detect wavelength ratios from params YAML in HDF5
    _lambda2, _lambda3 = _extract_wavelength_ratios(layers_hdf5)
    if lambda2 is None:
        lambda2 = _lambda2
    if lambda3 is None:
        lambda3 = _lambda3

    # Read optimizer output
    with h5py.File(layers_hdf5, 'r') as f:
        n_out = np.array(f['diffevo_output/n'])
        L_opt = np.array(f['diffevo_output/L'])

    L_phys_frac = op2phys(L_opt, n_out[1:-1])  # physical thickness as fraction of lambda

    # Build reference stack (for material properties and wavelength)
    _ref_stack = _build_mc_stack(n_out, L_phys_frac, layers_hdf5)
    _wavelength = _ref_stack['lam_ref']

    # Convert to SI metres for Brownian noise
    L_phys_m = L_phys_frac * _wavelength

    from OptimalBragg.noise import coating_thermooptic_fast, extract_stack_params
    _stack_params = extract_stack_params(_ref_stack, r_mirror, d_mirror)

    # Generate perturbation samples (3D: high-n, low-n, thickness)
    perturbs = _generate_perturbations(n_samples, n_dim=3)
    perturb_factors = 1 + perturbs  # (n_samples, 3)

    # Pre-build perturbed arrays
    # Identify low-n and high-n layers from the refractive indices.
    # The threshold is the geometric mean of the two film indices.
    n_layers_only = n_out[1:-1]
    n_thresh = np.sqrt(n_layers_only.min() * n_layers_only.max())
    low_mask = np.zeros(len(n_out), dtype=bool)
    high_mask = np.zeros(len(n_out), dtype=bool)
    for ii in range(1, len(n_out) - 1):
        if n_out[ii] < n_thresh:
            low_mask[ii] = True
        else:
            high_mask[ii] = True

    n_all = np.tile(n_out, (n_samples, 1))
    n_all[:, low_mask] *= perturb_factors[:, 1:2]   # low-n perturbation
    n_all[:, high_mask] *= perturb_factors[:, 0:1]   # high-n perturbation

    L_m_all = L_phys_m[None, :] * perturb_factors[:, 2:3]  # metres

    # Allocate output arrays
    Tp_IR = np.empty(n_samples)
    Tp_AUX = np.empty(n_samples)
    Tp_lam3 = np.empty(n_samples) if lambda3 is not None else None
    surfField = np.empty(n_samples)
    TOnoise = np.empty((n_samples, len(freq)))
    Brnoise = np.empty((n_samples, len(freq)))

    for jj in tqdm.tqdm(range(n_samples), desc='MC samples'):
        n_j = n_all[jj]
        L_m_j = L_m_all[jj]                       # physical thickness [m]
        L_opt_j = L_m_j * n_j[1:-1] / _wavelength  # optical thickness [frac of lambda]

        # Transmission at PSL wavelength
        r_psl, _ = multidiel1(n_j, L_opt_j, 1.0)
        Tp_IR[jj] = 1 - np.abs(r_psl) ** 2

        # Surface E-field
        surfField[jj] = np.abs(1 + r_psl[0]) ** 2

        # Transmission at AUX wavelength
        r_aux, _ = multidiel1(n_j, L_opt_j, lambda2)
        Tp_AUX[jj] = 1 - np.abs(r_aux) ** 2

        # Transmission at third wavelength (if present)
        if lambda3 is not None:
            r_lam3, _ = multidiel1(n_j, L_opt_j, lambda3)
            Tp_lam3[jj] = 1 - np.abs(r_lam3) ** 2

        # Thermo-optic noise (JIT)
        for fi, f_val in enumerate(freq):
            TOnoise[jj, fi] = coating_thermooptic_fast(
                f_val, L_opt_j, _wavelength, w_beam, _stack_params)

        # Brownian noise — use perturbed stack
        _perturbed_stack = dict(_ref_stack)
        _perturbed_stack['Ls'] = L_m_j.copy()      # metres
        _perturbed_stack['ns'] = n_j.copy()
        SbrZ = coating_brownian(freq, _perturbed_stack, w_beam)
        Brnoise[jj] = SbrZ

    # Package output: T [ppm] for HR wavelength, R [ppm] for AR wavelengths
    idx100 = np.argmin(np.abs(freq - 100))
    rows = [
        1e6 * Tp_IR,               # T_PSL [ppm]
        1e6 * (1 - Tp_AUX),        # R_AUX [ppm]
    ]
    if lambda3 is not None:
        rows.append(1e6 * (1 - Tp_lam3))  # R_lam3 [ppm]
    rows.extend([
        np.sqrt(TOnoise[:, idx100]) * 1e21,  # S_TO [×1e-21]
        np.sqrt(Brnoise[:, idx100]) * 1e21,  # S_Br [×1e-21]
        surfField,                  # E_surface
    ])
    MCout = np.vstack(rows)

    return {
        'MCout': MCout,
        'TOnoise': np.sqrt(TOnoise),
        'Brnoise': np.sqrt(Brnoise),
        'freq': freq,
        'perturbs': perturbs,
        'has_lambda3': lambda3 is not None,
    }


def _build_mc_stack(n_out, L_phys, layers_hdf5):
    """Build a minimal stack dict for MC noise evaluation.

    Reads the params_file or materials_file reference from the HDF5
    to reconstruct the stack with correct material properties.
    """
    from OptimalBragg import qw_stack, Material, load_materials_yaml
    from OptimalBragg.io import yamlread
    from pathlib import Path
    import warnings

    # 1. Try reading params_file stored in the HDF5 by the optimizer
    params_path = None
    with h5py.File(layers_hdf5, 'r') as f:
        if 'params_file' in f:
            raw = np.array(f['params_file'])
            pf_str = raw.item().decode('utf-8') if isinstance(raw.item(), bytes) else str(raw.item())
            candidate = Path(pf_str)
            if candidate.exists():
                params_path = candidate

    # 2. Fall back to directory walking
    if params_path is None:
        hdf5_dir = Path(layers_hdf5).parent
        for candidate in [
            hdf5_dir.parent.parent / 'ETM_params.yml',
            hdf5_dir.parent.parent / 'ITM_params.yml',
        ]:
            if candidate.exists():
                params_path = candidate
                break

    if params_path is not None:
        params = yamlread(str(params_path))
        mat_file = params_path.parent / params['misc'].get(
            'materials_file', 'materials.yml')
        if mat_file.exists():
            materials = load_materials_yaml(str(mat_file))
            from OptimalBragg.materials import air
            hwcap = params['misc'].get('hwcap', '')
            n_layers = len(L_phys)
            Npairs = (n_layers - len(hwcap)) // 2
            stack = qw_stack(
                lam_ref=materials['laser']['wavelength'],
                substrate=materials['substrate'],
                superstrate=Material(air),
                thin_films=materials['thin_films'],
                pattern='LH' * Npairs,
                hwcap=hwcap,
            )
            # Override with the actual optimized thicknesses
            stack['Ls'] = L_phys * stack['lam_ref']   # fraction → meters
            stack['Ls_opt'] = L_phys * n_out[1:-1]    # physical fraction × n = optical fraction
            stack['ns'] = n_out.copy()
            return stack

    # Fallback: build from the n/L arrays with default aLIGO materials
    warnings.warn(
        f"No params YAML found for {layers_hdf5}; falling back to "
        "aLIGO defaults (1064 nm, SiO2/TiTa2O5). This is WRONG for Voyager.",
        stacklevel=2,
    )
    from OptimalBragg.materials import SiO2, TiTa2O5, FusedSilica, air
    n_layers = len(L_phys)
    Npairs = n_layers // 2
    stack = qw_stack(
        lam_ref=1064e-9,
        substrate=Material(FusedSilica),
        superstrate=Material(air),
        thin_films={'L': Material(SiO2), 'H': Material(TiTa2O5)},
        pattern='LH' * Npairs,
    )
    stack['Ls'] = L_phys * stack['lam_ref']   # fraction → meters
    stack['Ls_opt'] = L_phys * n_out[1:-1]    # physical fraction × n = optical fraction
    stack['ns'] = n_out.copy()
    return stack


def save_mc(result, output_path):
    """Save MC results to HDF5.

    Parameters
    ----------
    result : dict
        Output of :func:`run_mc`.
    output_path : str
        HDF5 output path.
    """
    with h5py.File(output_path, 'w') as f:
        f['MCout'] = result['MCout']
        f['TOnoise'] = result['TOnoise']
        f['Brnoise'] = result['Brnoise']
        f['freq'] = result['freq']
        f['perturbs'] = result['perturbs']


if __name__ == '__main__':
    import sys

    if len(sys.argv) < 4:
        print("Usage: python -m OptimalBragg.mc <layers.hdf5> "
              "<output_MC.hdf5> <n_samples>")
        sys.exit(1)

    layers_hdf5 = sys.argv[1]
    output_path = sys.argv[2]
    n_samples = int(sys.argv[3])

    result = run_mc(layers_hdf5, n_samples=n_samples)
    save_mc(result, output_path)
    print(f"Saved {n_samples} MC samples to {output_path}")
