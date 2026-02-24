"""Integration tests: optimization, HDF5 I/O, plotting, MC, and reports."""
import os
import tempfile
import numpy as np
import pytest
from pathlib import Path
from scipy.optimize import differential_evolution as devo

from OptimalBragg import Material, qw_stack
from OptimalBragg.materials import SiO2, TiTa2O5, FusedSilica, air
from OptimalBragg.costs import getMirrorCost, precompute_misc
from OptimalBragg.noise import brownian_proxy


@pytest.mark.slow
def test_short_optimization():
    """Run a tiny differential_evolution to verify end-to-end pipeline."""
    Npairs = 3
    stack = qw_stack(
        lam_ref=1064e-9,
        substrate=Material(FusedSilica),
        superstrate=Material(air),
        thin_films={"L": Material(SiO2), "H": Material(TiTa2O5)},
        pattern="LH" * Npairs,
    )
    gam = brownian_proxy(stack)

    costs = {
        'Trans1': {'target': 5e-6, 'weight': 5},
        'Brownian': {'target': 20.0, 'weight': 2},
    }
    misc = {
        'pol': 'te', 'aoi': 0,
        'Npairs': Npairs, 'Nfixed': 0, 'Ncopies': 0,
        'lambda2': 0.5,
    }

    precompute_misc(costs, stack, misc)

    n_vars = 2 * Npairs
    bounds = ((0.05, 0.45),) * n_vars

    res = devo(
        func=getMirrorCost,
        bounds=bounds,
        updating='deferred',
        strategy='best1bin',
        popsize=5,
        maxiter=10,
        workers=1,
        args=(costs, stack, gam, False, misc),
        polish=False,
        disp=False,
    )

    assert np.isfinite(res.fun)
    assert len(res.x) == n_vars
    assert all(np.isfinite(res.x))


@pytest.mark.slow
def test_full_pipeline_hdf5_plot_mc_report():
    """End-to-end: optimize → HDF5 → read back → plot → MC → corner → RST.

    Verifies:
    - HDF5 has T1, vectorCost/* as readable datasets (not just attrs)
    - All 4 figure files created on disk
    - MC output is readable
    - RST has all required sections
    """
    import h5py
    import matplotlib
    matplotlib.use('Agg')

    from OptimalBragg.io import h5write, yamlread
    from OptimalBragg.layers import op2phys
    from OptimalBragg.plot import (apply_style, plot_layers, plot_spectral,
                                   plot_noise, plot_starfish)

    # ── Setup ──
    Npairs = 4
    wavelength = 1064e-9
    stack = qw_stack(
        lam_ref=wavelength,
        substrate=Material(FusedSilica),
        superstrate=Material(air),
        thin_films={"L": Material(SiO2), "H": Material(TiTa2O5)},
        pattern="LH" * Npairs,
    )
    gam = brownian_proxy(stack)

    costs = {
        'Trans1': {'target': 5e-6, 'weight': 5},
        'Brownian': {'target': 20.0, 'weight': 2},
    }
    misc = {
        'pol': 'te', 'aoi': 0,
        'Npairs': Npairs, 'Nfixed': 0, 'Ncopies': 0,
        'lambda2': 0.5,
    }
    precompute_misc(costs, stack, misc)

    n_vars = 2 * Npairs
    bounds = ((0.05, 0.45),) * n_vars

    # ── 1. Optimize ──
    res = devo(
        func=getMirrorCost,
        bounds=bounds,
        updating='deferred',
        strategy='best1bin',
        popsize=20,
        maxiter=15,
        workers=1,
        args=(costs, stack, gam, False, misc),
        polish=False,
        disp=False,
    )
    assert np.isfinite(res.fun)

    # Final verbose evaluation to get output dict
    scalar_cost, output = getMirrorCost(
        res.x, costs=costs, stack=stack, gam=gam,
        verbose=True, misc=misc,
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        ts = '260222_120000'
        optic = 'ETM'

        # ── 2. Write HDF5 ──
        vector_cost = dict(output.pop('vectorCost'))
        hdf5_path = os.path.join(tmpdir, f'{optic}_Layers_{ts}.hdf5')

        h5_dict = {
            'trajectory': np.array([scalar_cost]),
            'vec_evol': np.array([]),
            'diffevo_output': {
                'vectorCost': {k: np.float64(v)
                               for k, v in vector_cost.items()},
                'L': res.x,
                'n': stack['ns'],
                **{k: (np.float64(v) if np.isscalar(v) else v)
                   for k, v in output.items()},
            },
        }
        h5write(hdf5_path, h5_dict)

        # ── 3. Read back HDF5 and verify datasets ──
        with h5py.File(hdf5_path, 'r') as f:
            # T1 must be a readable dataset, not just an attr
            T1 = float(np.array(f['diffevo_output/T1']))
            assert np.isfinite(T1)

            # vectorCost entries must be datasets
            for cost_name in vector_cost:
                val = float(np.array(
                    f[f'diffevo_output/vectorCost/{cost_name}']))
                assert np.isfinite(val)

            # Scalar cost
            sc = float(np.array(f['diffevo_output/scalarCost']))
            assert np.isfinite(sc)

            L = np.array(f['diffevo_output/L'])
            n = np.array(f['diffevo_output/n'])

        assert len(L) == n_vars
        assert len(n) == n_vars + 2  # incident + substrate

        # ── 4. Plot all 4 figures ──
        apply_style()
        fig_dir = os.path.join(tmpdir, 'Figures', optic)
        os.makedirs(fig_dir, exist_ok=True)

        layers_path = os.path.join(fig_dir, f'{optic}_Layers_{ts}.svg')
        plot_layers(n, L, wavelength, save_path=layers_path)
        assert os.path.isfile(layers_path)

        # Starfish
        stats = {}
        for cname, cval in vector_cost.items():
            if cval > 0:
                stat = 1 / (1e-11 + cval)
                if cval > 1e3 or cval < 1e-5:
                    stat = 1e-2
                stats[cname] = np.abs(np.log(np.abs(stat)))

        sf_path = os.path.join(fig_dir, f'{optic}_SF{ts}.svg')
        if stats:
            plot_starfish(stats, scale=10, save_path=sf_path,
                          title=f'{optic}: {len(L)} layers')
        assert os.path.isfile(sf_path), "Starfish SVG not created"

        spectral_path = os.path.join(fig_dir, f'{optic}_R{ts}.svg')
        plot_spectral(n, L, wavelength, T1=T1,
                      lambda2=0.5, save_path=spectral_path)
        assert os.path.isfile(spectral_path)

        tn_path = os.path.join(fig_dir, f'{optic}_TN_{ts}.svg')
        stack['Ls_opt'] = L.copy()
        stack['Ls'] = wavelength * op2phys(L, n[1:-1])
        stack['ns'] = n.copy()
        plot_noise(stack, 0.062, save_path=tn_path)
        assert os.path.isfile(tn_path)

        import matplotlib.pyplot as plt
        plt.close('all')

        # ── 5. MC (tiny sample) ──
        from OptimalBragg.mc import run_mc, save_mc
        mc_output = os.path.join(tmpdir, f'{optic}_MC.hdf5')
        mc_result = run_mc(hdf5_path, n_samples=50, lambda2=0.5)
        save_mc(mc_result, mc_output)

        with h5py.File(mc_output, 'r') as f:
            mc_samples = np.array(f['MCout'][:])
        assert mc_samples.shape[1] == 50

        # ── 6. Corner plot ──
        from OptimalBragg.plot import plot_corner
        corner_path = os.path.join(fig_dir,
                                   f'{optic}_nominal_cornerPlt.svg')
        plot_corner(mc_samples, mirror_type=optic, save_path=corner_path)
        assert os.path.isfile(corner_path)
        plt.close('all')

        # ── 7. Generate RST report ──
        import OptimalBragg.report as _report_mod
        from OptimalBragg.report import generate_run_rst

        # Redirect report output to tmpdir so we don't pollute real docs/
        orig_root = _report_mod._PROJECT_ROOT
        _report_mod._PROJECT_ROOT = tmpdir

        try:
            rel_fig = f'Figures/{optic}'
            fig_paths = {
                'layers': f'{rel_fig}/{optic}_Layers_{ts}.svg',
                'spectral': f'{rel_fig}/{optic}_R{ts}.svg',
                'starfish': f'{rel_fig}/{optic}_SF{ts}.svg',
                'thermal_noise': f'{rel_fig}/{optic}_TN_{ts}.svg',
                'corner': f'{rel_fig}/{optic}_nominal_cornerPlt.svg',
            }

            params = {
                'costs': {
                    'Trans1': {'target': 5e-6, 'weight': 5},
                    'Brownian': {'target': 20.0, 'weight': 2},
                },
                '_wavelength': wavelength,
                'misc': misc,
            }

            rst_content, rst_path = generate_run_rst(
                hdf5_path=hdf5_path,
                mirror_type=optic,
                fig_paths=fig_paths,
                params=params,
                mc_hdf5_path=mc_output,
            )
        finally:
            _report_mod._PROJECT_ROOT = orig_root

        # Verify RST content
        assert 'Design Summary' in rst_content
        assert 'T @ 1064 nm' in rst_content
        assert 'Layer Structure' in rst_content
        assert 'Spectral Reflectivity' in rst_content
        assert 'Thermal Noise' in rst_content
        assert 'Monte Carlo' in rst_content
        # Should have image directives for all 5 figures
        assert rst_content.count('.. image::') == 5
