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


def _generate_perturbations(n_samples, n_dim=3, n_walkers=64,
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
    pos, _, _ = sampler.run_mcmc(p0, 1000, progress=False)
    sampler.reset()

    # Production
    n_steps = max(n_samples // n_walkers + 1, 100)
    sampler.run_mcmc(pos, n_steps, progress=False)

    return sampler.get_chain(flat=True)[:n_samples, :]


def run_mc(layers_hdf5, n_samples=5000, lambda_aux=0.5,
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
    lambda_aux : float
        AUX wavelength ratio (lambda_AUX / lambda_PSL).
    w_beam : float
        Beam radius [m].
    r_mirror, d_mirror : float
        Mirror radius and thickness [m].
    freq : array_like, optional
        Frequency array [Hz]. Default ``logspace(1, 3, 50)``.

    Returns
    -------
    result : dict
        Keys: ``MCout`` (5 x n_samples), ``TOnoise``, ``Brnoise``,
        ``freq``, ``perturbs``.
    """
    from OptimalBragg.noise import coating_thermooptic, coating_brownian

    if freq is None:
        freq = np.logspace(1, 3, 50)

    # Read optimizer output
    with h5py.File(layers_hdf5, 'r') as f:
        n_out = np.array(f['diffevo_output/n'])
        L_opt = np.array(f['diffevo_output/L'])

    L_phys = op2phys(L_opt, n_out[1:-1])

    # Generate perturbation samples (3D: high-n, low-n, thickness)
    perturbs = _generate_perturbations(n_samples, n_dim=3)
    perturb_factors = 1 + perturbs  # (n_samples, 3)

    # Pre-build perturbed arrays
    # n_out layout: [superstrate, L, H, L, H, ..., substrate]
    # Even indices (1,3,5,...) = low-n, Odd indices (2,4,6,...) = high-n
    n_all = np.tile(n_out, (n_samples, 1))
    n_all[:, 1:-1:2] *= perturb_factors[:, 1:2]   # low-n perturbation
    n_all[:, 2::2] *= perturb_factors[:, 0:1]      # high-n perturbation

    L_all = L_phys[None, :] * perturb_factors[:, 2:3]  # thickness

    # Allocate output arrays
    Tp_IR = np.empty(n_samples)
    Tp_AUX = np.empty(n_samples)
    surfField = np.empty(n_samples)
    TOnoise = np.empty((n_samples, len(freq)))
    Brnoise = np.empty((n_samples, len(freq)))

    for jj in tqdm.tqdm(range(n_samples), desc='MC samples'):
        n_j = n_all[jj]
        L_j = L_all[jj]
        L_opt_j = L_j * n_j[1:-1]  # physical → optical

        # Transmission at PSL wavelength
        r_psl, _ = multidiel1(n_j, L_opt_j, 1.0)
        Tp_IR[jj] = 1 - np.abs(r_psl) ** 2

        # Surface E-field
        surfField[jj] = np.abs(1 + r_psl[0]) ** 2

        # Transmission at AUX wavelength
        r_aux, _ = multidiel1(n_j, L_opt_j, lambda_aux)
        Tp_AUX[jj] = 1 - np.abs(r_aux) ** 2

        # Thermal noise — build a temporary stack dict
        # We need a minimal stack for the noise functions
        # For speed, use the JIT thermo-optic directly
        from OptimalBragg.noise import (
            coating_thermooptic_fast, extract_stack_params,
        )
        # Build minimal stack for this perturbed sample
        # (Only needed for extract_stack_params on first iteration;
        #  after that we can reuse the params structure)
        if jj == 0:
            # Build a reference stack for extract_stack_params
            # We need substrate/material properties — read from the
            # original stack. For MC we only perturb n and L, not
            # material thermal properties.
            from OptimalBragg.io import yamlread
            _ref_stack = _build_mc_stack(n_out, L_phys, layers_hdf5)
            _stack_params = extract_stack_params(
                _ref_stack, r_mirror, d_mirror)
            _wavelength = _ref_stack['lam_ref']

        # Thermo-optic noise (JIT)
        for fi, f_val in enumerate(freq):
            TOnoise[jj, fi] = coating_thermooptic_fast(
                f_val, L_opt_j, _wavelength, w_beam, _stack_params)

        # Brownian noise
        SbrZ = coating_brownian(freq, _ref_stack, w_beam)
        Brnoise[jj] = SbrZ

    # Package output in same format as legacy doMC.py
    idx100 = np.argmin(np.abs(freq - 100))
    MCout = np.vstack((
        1e6 * Tp_IR,               # T_1064 [ppm]
        1e2 * Tp_AUX,              # T_AUX [%]
        np.sqrt(TOnoise[:, idx100]) * 1e21,  # S_TO [×1e-21]
        np.sqrt(Brnoise[:, idx100]) * 1e21,  # S_Br [×1e-21]
        surfField,                  # E_surface
    ))

    return {
        'MCout': MCout,
        'TOnoise': np.sqrt(TOnoise),
        'Brnoise': np.sqrt(Brnoise),
        'freq': freq,
        'perturbs': perturbs,
    }


def _build_mc_stack(n_out, L_phys, layers_hdf5):
    """Build a minimal stack dict for MC noise evaluation.

    Reads the params_file or materials_file reference from the HDF5
    to reconstruct the stack with correct material properties.
    """
    from OptimalBragg import qw_stack, Material, load_materials_yaml
    from OptimalBragg.io import yamlread
    from pathlib import Path
    import os

    # Try to find materials config from the HDF5 metadata
    hdf5_dir = Path(layers_hdf5).parent

    # Walk up to find a params YAML
    params_path = None
    for candidate in [
        hdf5_dir.parent.parent / 'ETM_params.yml',
        hdf5_dir.parent.parent / 'ITM_params.yml',
        hdf5_dir / '../../ETM_params.yml',
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
            n_layers = len(L_phys)
            Npairs = n_layers // 2
            stack = qw_stack(
                lam_ref=materials['laser']['wavelength'],
                substrate=materials['substrate'],
                superstrate=Material(air),
                thin_films=materials['thin_films'],
                pattern='LH' * Npairs,
            )
            # Override with the actual optimized thicknesses
            stack['Ls'] = L_phys.copy()
            stack['Ls_opt'] = L_phys * n_out[1:-1] / stack['lam_ref']
            stack['ns'] = n_out.copy()
            return stack

    # Fallback: build from the n/L arrays with default materials
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
    stack['Ls'] = L_phys.copy()
    stack['Ls_opt'] = L_phys * n_out[1:-1] / stack['lam_ref']
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
