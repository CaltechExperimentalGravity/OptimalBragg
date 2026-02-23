"""Command-line interface for OptimalBragg.

Usage::

    optimalbragg optimize projects/aLIGO/ETM_params.yml
    optimalbragg plot projects/aLIGO/ETM_params.yml
    optimalbragg mc Data/ETM/ETM_Layers_260221.hdf5 5000
    optimalbragg corner Data/ETM/ETM_MC.hdf5
"""

import argparse
import sys


def cmd_optimize(args):
    """Run coating optimization."""
    from OptimalBragg.optimizer import run_optimization
    run_optimization(
        args.params,
        save=not args.no_save,
        optic=args.optic,
    )


def cmd_plot(args):
    """Generate plots from optimization output."""
    import glob
    import os
    from pathlib import Path

    import numpy as np
    import h5py
    import matplotlib
    matplotlib.use('Agg')

    from OptimalBragg.io import yamlread, load_materials_yaml
    from OptimalBragg.layers import op2phys
    from OptimalBragg.plot import apply_style, plot_layers, plot_spectral
    from OptimalBragg.plot import plot_noise, plot_starfish

    params_path = Path(args.params)
    params_dir = params_path.parent
    opt_params = yamlread(str(params_path))
    misc = opt_params['misc']

    # Infer optic name
    stem = params_path.stem.upper()
    optic = 'ETM' if 'ETM' in stem else ('ITM' if 'ITM' in stem else 'output')

    # Find the HDF5 file
    data_dir = params_dir / 'Data' / optic
    if args.hdf5:
        hdf5_path = args.hdf5
    else:
        candidates = sorted(glob.glob(str(data_dir / '*Layers*.hdf5')),
                            key=os.path.getctime)
        if not candidates:
            print(f"No HDF5 files found in {data_dir}")
            sys.exit(1)
        hdf5_path = candidates[-1]
        print(f"Using: {hdf5_path}")

    # Load materials
    mat_file = params_dir / misc.get('materials_file', 'materials.yml')
    materials = load_materials_yaml(str(mat_file))

    # Read optimizer output
    with h5py.File(hdf5_path, 'r') as f:
        n = np.array(f['diffevo_output/n'])
        L_opt = np.array(f['diffevo_output/L'])
        T1064 = float(np.array(f['diffevo_output/T1064'])) \
            if 'diffevo_output/T1064' in f else None

    wavelength = materials['laser']['wavelength']
    fig_dir = params_dir / 'Figures' / optic
    os.makedirs(fig_dir, exist_ok=True)

    ts = Path(hdf5_path).stem[-13:]  # YYMMDD_HHMMSS

    apply_style()

    # Layer structure
    name_high = list(materials['thin_films'].keys())[-1]
    name_low = list(materials['thin_films'].keys())[0]
    plot_layers(n, L_opt, wavelength,
                name_high=name_high, name_low=name_low,
                save_path=str(fig_dir / f'{optic}_Layers_{ts}.svg'))
    print(f"  Saved: {optic}_Layers_{ts}.svg")

    # Starfish
    stats = {}
    with h5py.File(hdf5_path, 'r') as f:
        for cost in opt_params['costs']:
            if opt_params['costs'][cost]['weight']:
                try:
                    cv = float(np.array(
                        f[f'diffevo_output/vectorCost/{cost}']))
                    stat = 1 / (1e-11 + cv)
                    if cv > 1e3 or cv < 1e-5:
                        stat = 1e-2
                    stats[cost] = np.abs(np.log(np.abs(stat)))
                except KeyError:
                    pass

    sf_path = str(fig_dir / f'{optic}_SF{ts}.svg')
    if stats:
        plot_starfish(stats, scale=10, save_path=sf_path,
                      title=f'{optic}: {len(L_opt)} layers')
        print(f"  Saved: {optic}_SF{ts}.svg")

    # Spectral reflectivity
    lambda_aux = opt_params['misc'].get('lambdaAUX', 0.5)
    plot_spectral(n, L_opt, wavelength, T1064=T1064,
                  lambda_aux=lambda_aux,
                  save_path=str(fig_dir / f'{optic}_R{ts}.svg'))
    print(f"  Saved: {optic}_R{ts}.svg")

    # Thermal noise
    from OptimalBragg import Material, qw_stack
    from OptimalBragg.materials import air
    from OptimalBragg.noise import brownian_proxy

    hwcap = misc.get('hwcap', '')
    Npairs = misc.get('Npairs', (len(L_opt) - len(hwcap)) // 2)
    stack = qw_stack(
        lam_ref=wavelength,
        substrate=materials['substrate'],
        superstrate=Material(air),
        thin_films=materials['thin_films'],
        pattern='LH' * Npairs,
        hwcap=hwcap,
    )
    # Override with optimized thicknesses
    stack['Ls_opt'] = L_opt.copy()
    stack['Ls'] = wavelength * op2phys(L_opt, n[1:-1])
    stack['ns'] = n.copy()

    sub_overrides = {}
    raw_mat = yamlread(str(mat_file))
    if 'substrate' in raw_mat and 'overrides' in raw_mat['substrate']:
        sub_overrides = raw_mat['substrate']['overrides']
    r_mirror = sub_overrides.get('MassRadius', 0.17)
    d_mirror = sub_overrides.get('MassThickness', 0.20)
    optic_cfg = materials.get('optics', {}).get(optic, {})
    w_beam = optic_cfg.get('beam_radius', 0.062)

    plot_noise(stack, w_beam, r_mirror=r_mirror, d_mirror=d_mirror,
               save_path=str(fig_dir / f'{optic}_TN_{ts}.svg'))
    print(f"  Saved: {optic}_TN_{ts}.svg")

    import matplotlib.pyplot as plt
    plt.close('all')
    print("Done.")


def cmd_mc(args):
    """Run Monte Carlo sensitivity analysis."""
    from OptimalBragg.mc import run_mc, save_mc

    result = run_mc(args.hdf5, n_samples=args.n_samples)
    save_mc(result, args.output)
    print(f"Saved {args.n_samples} MC samples to {args.output}")


def cmd_sweep(args):
    """Sweep Npairs to find minimum that hits all targets."""
    from OptimalBragg.optimizer import sweep_nlayers

    n_range = None
    if args.min_pairs is not None and args.max_pairs is not None:
        n_range = (args.min_pairs, args.max_pairs)

    sweep_nlayers(
        args.params,
        n_range=n_range,
        save=not args.no_save,
        optic=args.optic,
    )


def cmd_corner(args):
    """Generate corner plot from MC output."""
    import matplotlib
    matplotlib.use('Agg')
    import h5py
    import numpy as np
    from OptimalBragg.plot import plot_corner

    with h5py.File(args.hdf5, 'r') as f:
        mc_samples = np.array(f['MCout'][:])

    mirror_type = args.mirror_type
    if mirror_type is None:
        mirror_type = 'ITM' if 'ITM' in args.hdf5 else 'ETM'

    output = args.output
    if output is None:
        import os
        fig_dir = f'Figures/{mirror_type}'
        os.makedirs(fig_dir, exist_ok=True)
        output = f'{fig_dir}/{mirror_type}_nominal_cornerPlt.svg'

    plot_corner(mc_samples, mirror_type=mirror_type, save_path=output)
    print(f"Saved corner plot: {output}")


def main():
    parser = argparse.ArgumentParser(
        prog='optimalbragg',
        description='OptimalBragg — optical coating design and optimization',
    )
    subparsers = parser.add_subparsers(dest='command', required=True)

    # optimize
    p_opt = subparsers.add_parser('optimize', help='Run coating optimization')
    p_opt.add_argument('params', help='Path to cost/optimizer YAML')
    p_opt.add_argument('--optic', help='Optic name (ETM or ITM)')
    p_opt.add_argument('--no-save', action='store_true',
                       help='Skip HDF5 output')
    p_opt.set_defaults(func=cmd_optimize)

    # plot
    p_plot = subparsers.add_parser('plot', help='Generate plots')
    p_plot.add_argument('params', help='Path to cost/optimizer YAML')
    p_plot.add_argument('--hdf5', help='Specific HDF5 file (default: latest)')
    p_plot.set_defaults(func=cmd_plot)

    # mc
    p_mc = subparsers.add_parser('mc', help='Monte Carlo sensitivity')
    p_mc.add_argument('hdf5', help='Path to optimizer output HDF5')
    p_mc.add_argument('output', help='Output HDF5 path')
    p_mc.add_argument('n_samples', type=int, help='Number of MC samples')
    p_mc.set_defaults(func=cmd_mc)

    # sweep
    p_sweep = subparsers.add_parser('sweep',
                                     help='Sweep Npairs to find optimal layer count')
    p_sweep.add_argument('params', help='Path to cost/optimizer YAML')
    p_sweep.add_argument('--optic', help='Optic name (ETM or ITM)')
    p_sweep.add_argument('--min-pairs', type=int, help='Min bilayer count')
    p_sweep.add_argument('--max-pairs', type=int, help='Max bilayer count')
    p_sweep.add_argument('--no-save', action='store_true',
                         help='Skip HDF5 output')
    p_sweep.set_defaults(func=cmd_sweep)

    # corner
    p_corner = subparsers.add_parser('corner', help='Corner plot from MC')
    p_corner.add_argument('hdf5', help='Path to MC output HDF5')
    p_corner.add_argument('--output', '-o', help='Output figure path')
    p_corner.add_argument('--mirror-type', help='ETM or ITM')
    p_corner.set_defaults(func=cmd_corner)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
