#!/usr/bin/env python
"""Full pipeline demo: optimize + plot + MC + corner + Sphinx report.

Runs ETM and ITM for aLIGO, generates all outputs, builds HTML docs.

Usage:
    source ~/.zshrc && conda activate coatingDev
    python run_full_pipeline.py
"""
import os
import sys
import glob
import time
from pathlib import Path

import matplotlib
matplotlib.use('Agg')

PROJECT = 'projects/aLIGO'


def find_latest_hdf5(data_dir):
    """Find the most recently created HDF5 file."""
    candidates = sorted(glob.glob(str(data_dir / '*Layers*.hdf5')),
                        key=os.path.getctime)
    if not candidates:
        raise FileNotFoundError(f"No HDF5 files in {data_dir}")
    return candidates[-1]


def run_pipeline(optic):
    """Run full pipeline for one optic (ETM or ITM)."""
    params_yml = f'{PROJECT}/{optic}_params.yml'
    print(f"\n{'='*70}")
    print(f"  {optic}: FULL PIPELINE")
    print(f"{'='*70}")

    # --- 1. Optimize ---
    print(f"\n[1/5] Optimizing {optic}...")
    t0 = time.perf_counter()

    from OptimalBragg.optimizer import run_optimization
    result = run_optimization(params_yml)

    dt_opt = time.perf_counter() - t0
    print(f"  Optimization complete in {dt_opt:.1f}s")
    print(f"  Final cost: {result['scalar_cost']:.6f}")
    print(f"  T_1: {result['output']['T1']*1e6:.2f} ppm")
    if 'T2' in result['output']:
        print(f"  T_2:  {result['output']['T2']*100:.2f}%")

    # Find the HDF5 file just created
    data_dir = Path(PROJECT) / 'Data' / optic
    hdf5_path = find_latest_hdf5(data_dir)
    print(f"  Saved: {hdf5_path}")

    # --- 2. Plot ---
    print(f"\n[2/5] Generating plots for {optic}...")
    t0 = time.perf_counter()

    sys.argv = ['optimalbragg', 'plot', params_yml, '--hdf5', hdf5_path]
    from OptimalBragg.cli import cmd_plot
    import argparse
    args = argparse.Namespace(params=params_yml, hdf5=hdf5_path)
    cmd_plot(args)

    dt_plot = time.perf_counter() - t0
    print(f"  Plots complete in {dt_plot:.1f}s")

    # --- 3. Monte Carlo ---
    print(f"\n[3/5] Running MC for {optic} (2000 samples)...")
    t0 = time.perf_counter()

    mc_output = str(data_dir / f'{optic}_MC.hdf5')

    from OptimalBragg.io import yamlread
    opt_params = yamlread(params_yml)
    misc = opt_params['misc']
    lambda2 = misc.get('lambda2', 0.5)
    lambda3 = misc.get('lambda3', None)

    from OptimalBragg.mc import run_mc, save_mc
    mc_result = run_mc(hdf5_path, n_samples=2000, lambda2=lambda2,
                       lambda3=lambda3)
    save_mc(mc_result, mc_output)

    dt_mc = time.perf_counter() - t0
    print(f"  MC complete in {dt_mc:.1f}s")
    print(f"  Saved: {mc_output}")

    # --- 4. Corner plot ---
    print(f"\n[4/5] Generating corner plot for {optic}...")
    t0 = time.perf_counter()

    fig_dir = Path(PROJECT) / 'Figures' / optic
    corner_path = str(fig_dir / f'{optic}_nominal_cornerPlt.svg')

    from OptimalBragg.plot import plot_corner
    import h5py
    import numpy as np
    with h5py.File(mc_output, 'r') as f:
        mc_samples = np.array(f['MCout'][:])
    from OptimalBragg.io import load_materials_yaml
    mat_file = str(Path(PROJECT) / misc.get('materials_file', 'materials.yml'))
    materials = load_materials_yaml(mat_file)
    wl_info = {
        'wavelength': materials['laser']['wavelength'],
        'lambda2': lambda2,
        'lambda3': lambda3,
    }
    plot_corner(mc_samples, mirror_type=optic, save_path=corner_path,
                wavelength_info=wl_info)

    dt_corner = time.perf_counter() - t0
    print(f"  Corner plot complete in {dt_corner:.1f}s")
    print(f"  Saved: {corner_path}")

    # --- 5. Sphinx report ---
    print(f"\n[5/5] Generating Sphinx report for {optic}...")
    t0 = time.perf_counter()

    ts = Path(hdf5_path).stem[-13:]  # YYMMDD_HHMMSS
    fig_paths = {
        'layers': f'{PROJECT}/Figures/{optic}/{optic}_Layers_{ts}.svg',
        'spectral': f'{PROJECT}/Figures/{optic}/{optic}_R{ts}.svg',
        'starfish': f'{PROJECT}/Figures/{optic}/{optic}_SF{ts}.svg',
        'thermal_noise': f'{PROJECT}/Figures/{optic}/{optic}_TN_{ts}.svg',
        'corner': corner_path,
    }

    from OptimalBragg.report import generate_run_rst
    rst_content, rst_path = generate_run_rst(
        hdf5_path=hdf5_path,
        mirror_type=optic,
        fig_paths=fig_paths,
        params=opt_params,
        project_dir=PROJECT,
        mc_hdf5_path=mc_output,
    )

    dt_report = time.perf_counter() - t0
    print(f"  Report complete in {dt_report:.1f}s")
    print(f"  RST: {rst_path}")

    return {
        'hdf5': hdf5_path,
        'mc_hdf5': mc_output,
        'corner': corner_path,
        'rst': rst_path,
        'dt_opt': dt_opt,
        'dt_plot': dt_plot,
        'dt_mc': dt_mc,
        'dt_corner': dt_corner,
        'dt_report': dt_report,
    }


def main():
    t_total = time.perf_counter()

    # Run ETM pipeline
    etm = run_pipeline('ETM')

    # Run ITM pipeline
    itm = run_pipeline('ITM')

    dt_total = time.perf_counter() - t_total

    # Summary
    print(f"\n{'='*70}")
    print(f"  PIPELINE COMPLETE — {dt_total:.0f}s total")
    print(f"{'='*70}")
    print(f"\n  ETM: {etm['dt_opt']:.0f}s opt + {etm['dt_mc']:.0f}s MC")
    print(f"  ITM: {itm['dt_opt']:.0f}s opt + {itm['dt_mc']:.0f}s MC")
    print(f"\n  Outputs:")
    print(f"    {etm['hdf5']}")
    print(f"    {itm['hdf5']}")
    print(f"    {etm['mc_hdf5']}")
    print(f"    {itm['mc_hdf5']}")
    print(f"    {etm['corner']}")
    print(f"    {itm['corner']}")
    print(f"    {etm['rst']}")
    print(f"    {itm['rst']}")

    # Check for Sphinx HTML
    html_dir = Path('docs/_build/html')
    if html_dir.exists():
        print(f"\n  Sphinx HTML: {html_dir.resolve()}")
        etm_ts = Path(etm['hdf5']).stem[-13:]
        itm_ts = Path(itm['hdf5']).stem[-13:]
        for ts in [etm_ts, itm_ts]:
            html_file = html_dir / 'runs' / f'{ts}.html'
            if html_file.exists():
                print(f"    {html_file}")
    else:
        print(f"\n  (Sphinx HTML build may have failed - check above)")


if __name__ == '__main__':
    main()
