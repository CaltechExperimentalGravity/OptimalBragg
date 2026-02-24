"""Trade study: green finesse vs 1064 nm thermal noise for aLIGO ETM.

Reoptimizes the ETM coating for T_532 targets of 2%, 0.5%, 0.1%.
Compares thermal noise spectra, absorption, and T_1064 achieved.

Usage:
    python greenFinesse.py
"""
import copy
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.optimize import differential_evolution as devo

from OptimalBragg import Material, qw_stack, load_materials_yaml
from OptimalBragg.io import yamlread
from OptimalBragg.materials import air
from OptimalBragg.layers import op2phys, fieldDepth
from OptimalBragg.costs import getMirrorCost, precompute_misc
from OptimalBragg.noise import brownian_proxy, coating_thermooptic, coating_brownian

targets_532 = [0.02, 0.005, 0.001]  # 2%, 0.5%, 0.1%


def run_study():
    opt_params = yamlread('ETM_params.yml')
    materials = load_materials_yaml(opt_params['misc']['materials_file'])
    lambdaPSL = materials['laser']['wavelength']
    Npairs = opt_params['misc']['Npairs']
    stack = qw_stack(
        lam_ref=lambdaPSL,
        substrate=materials['substrate'],
        superstrate=Material(air),
        thin_films=materials['thin_films'],
        pattern='LH' * Npairs,
    )
    gam = brownian_proxy(stack)
    w_beam = materials.get('optics', {}).get('ETM', {}).get('beam_radius', 0.062)
    r_mirror = 0.17
    d_mirror = 0.20

    results = {}

    for t532 in targets_532:
        print(f"\n{'='*60}")
        print(f"  Optimizing for T_532 = {t532*100:.1f}%")
        print(f"{'='*60}")

        params = copy.deepcopy(opt_params)
        params['costs']['Trans2']['target'] = t532
        precompute_misc(params['costs'], stack, params['misc'])

        nLayers = 2 * Npairs
        minThick = 10e-9 * stack['ns'][1] / lambdaPSL
        bounds = ((minThick, 0.48),) + ((0.05, 0.48),) * (nLayers - 1)

        res = devo(func=getMirrorCost, bounds=bounds,
                   updating='deferred', strategy='best1bin',
                   mutation=(0.05, 1.5),
                   popsize=params['misc']['Nparticles'],
                   init=params['misc']['init_method'],
                   workers=1, maxiter=2000,
                   atol=params['misc']['atol'],
                   tol=params['misc']['tol'],
                   args=(params['costs'], stack, gam, False, params['misc']),
                   polish=True, disp=True)

        # Get detailed output (no Ncopies/Nfixed for final eval)
        final_misc = copy.deepcopy(params['misc'])
        final_misc.update({'Ncopies': 0, 'Nfixed': 0})
        precompute_misc(params['costs'], stack, final_misc)
        _, output = getMirrorCost(res.x, params['costs'], stack, gam, True, final_misc)

        # Thermal noise — update stack with optimized thicknesses
        opt_stack = copy.deepcopy(stack)
        opt_stack['Ls_opt'] = res.x.copy()
        n = output['n']
        opt_stack['Ls'] = lambdaPSL * op2phys(res.x, n[1:-1])
        opt_stack['ns'] = n.copy()

        ff = np.logspace(0, 4, 500)
        StoZ, SteZ, StrZ = coating_thermooptic(
            ff, opt_stack, w_beam, r_mirror, d_mirror)
        SbrZ = coating_brownian(ff, opt_stack, w_beam)

        # Absorption
        L_phys = opt_stack['Ls']
        _, field = fieldDepth(L_phys, n, pol='p', nPts=300, lam=lambdaPSL)
        alpha_low, alpha_high = 0.5, 5.0
        alphas = np.where(np.arange(len(L_phys)) % 2 == 0, alpha_low, alpha_high)
        from OptimalBragg.layers import calc_abs
        absorption = calc_abs(field, L_phys, alphas)

        # CTN at 100 Hz
        idx100 = np.argmin(np.abs(ff - 100.0))
        CTN_100 = np.sqrt(StoZ[idx100] + SbrZ[idx100])
        Br_100 = np.sqrt(SbrZ[idx100])

        results[t532] = {
            'L_opt': res.x, 'n': n, 'cost': res.fun,
            'T1064': output['T1'], 'T532': output['T2'],
            'ff': ff, 'StoZ': StoZ, 'SbrZ': SbrZ,
            'SteZ': SteZ, 'StrZ': StrZ,
            'absorption': absorption,
            'CTN_100': CTN_100, 'Br_100': Br_100,
        }

    return results


def print_summary(results):
    print(f"\n{'='*72}")
    print(f"  Green Finesse Trade Study — Summary")
    print(f"{'='*72}")
    print(f"{'T_532 target':>14s}  {'T_1064 [ppm]':>13s}  {'T_532 achieved':>14s}  "
          f"{'CTN@100Hz':>12s}  {'Brownian':>12s}  {'Abs [ppm]':>10s}")
    print(f"{'-'*72}")
    for t532 in targets_532:
        r = results[t532]
        print(f"{t532*100:>12.1f}%  {r['T1']*1e6:>13.2f}  "
              f"{r['T2']*100:>12.2f}%  "
              f"{r['CTN_100']/1e-22:>10.2f}e-22  "
              f"{r['Br_100']/1e-22:>10.2f}e-22  "
              f"{r['absorption']:>10.1f}")
    print()


def make_plots(results):
    plt.rcParams.update({
        'text.usetex': False,
        'lines.linewidth': 3,
        'font.size': 18,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'axes.grid': True,
        'grid.alpha': 0.2,
        'legend.fontsize': 'small',
        'figure.figsize': (16, 7),
        'savefig.dpi': 140,
        'savefig.bbox': 'tight',
        'pdf.compression': 9,
    })

    fig, axes = plt.subplots(1, 3, figsize=(20, 7))

    # Panel 1: Thermal noise spectra
    colors = ['xkcd:blue', 'xkcd:orange', 'xkcd:red']
    for (t532, r), c in zip(results.items(), colors):
        CTN = np.sqrt(r['StoZ'] + r['SbrZ'])
        axes[0].loglog(r['ff'], CTN, c=c, lw=2.5,
                       label=f'T_532={t532*100:.1f}%')
        axes[0].loglog(r['ff'], np.sqrt(r['SbrZ']), c=c, lw=1, ls='--', alpha=0.5)
    axes[0].set_xlabel('Frequency [Hz]')
    axes[0].set_ylabel(R'Displacement noise [m/$\sqrt{\mathrm{Hz}}$]')
    axes[0].set_xlim([10, 1e4])
    axes[0].set_ylim([8e-24, 2e-20])
    axes[0].legend()
    axes[0].set_title('Coating Thermal Noise')

    # Panel 2: Bar chart — T_1064 and absorption
    labels = [f'{t*100:.1f}%' for t in targets_532]
    x = np.arange(len(labels))
    width = 0.35

    T1064 = [results[t]['T1'] * 1e6 for t in targets_532]
    absorp = [results[t]['absorption'] for t in targets_532]

    ax2a = axes[1]
    bars1 = ax2a.bar(x - width/2, T1064, width, label='T_1064 [ppm]',
                     color='xkcd:cerulean', alpha=0.8)
    ax2a.set_ylabel('T_1064 [ppm]')
    ax2a.set_xticks(x)
    ax2a.set_xticklabels(labels)
    ax2a.set_xlabel('T_532 target')

    ax2b = ax2a.twinx()
    bars2 = ax2b.bar(x + width/2, absorp, width, label='Absorption [ppm]',
                     color='xkcd:tomato', alpha=0.8)
    ax2b.set_ylabel('Absorption [ppm]')

    ax2a.legend(handles=[bars1, bars2], loc='upper left')
    axes[1].set_title('Transmission & Absorption')

    # Panel 3: CTN at 100 Hz
    CTN_vals = [results[t]['CTN_100'] / 1e-22 for t in targets_532]
    Br_vals = [results[t]['Br_100'] / 1e-22 for t in targets_532]
    axes[2].bar(x - width/2, CTN_vals, width, label='Total CTN',
                color='xkcd:dark grey', alpha=0.8)
    axes[2].bar(x + width/2, Br_vals, width, label='Brownian only',
                color='xkcd:tomato', alpha=0.8)
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(labels)
    axes[2].set_xlabel('T_532 target')
    axes[2].set_ylabel(R'Noise @ 100 Hz [$\times 10^{-22}$ m/$\sqrt{\mathrm{Hz}}$]')
    axes[2].legend()
    axes[2].set_title('CTN @ 100 Hz')

    fig.suptitle('Green Finesse Trade Study: aLIGO ETM', fontsize=20, y=1.02)
    fig.tight_layout()

    os.makedirs('Figures/ETM/', exist_ok=True)
    plt.savefig('Figures/ETM/greenFinesse_comparison.pdf')
    print("Saved: Figures/ETM/greenFinesse_comparison.pdf")


def save_hdf5(results):
    os.makedirs('Data/ETM/', exist_ok=True)
    fname = 'Data/ETM/greenFinesse_results.hdf5'
    with h5py.File(fname, 'w') as f:
        for t532, r in results.items():
            g = f.create_group(f'T2_{t532}')
            g['L_opt'] = r['L_opt']
            g['T1'] = r['T1']
            g['T2'] = r['T2']
            g['absorption'] = r['absorption']
            g['StoZ'] = r['StoZ']
            g['SbrZ'] = r['SbrZ']
            g['SteZ'] = r['SteZ']
            g['StrZ'] = r['StrZ']
            g['ff'] = r['ff']
            g['CTN_100'] = r['CTN_100']
            g['Br_100'] = r['Br_100']
    print(f"Saved: {fname}")


if __name__ == '__main__':
    results = run_study()
    print_summary(results)
    make_plots(results)
    save_hdf5(results)
