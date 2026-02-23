"""Corner plot from Monte Carlo sensitivity analysis using ArviZ.

Usage::

    python cornerPlt.py Data/ETM/ETM_MC.hdf5
    python cornerPlt.py Data/ITM/ITM_MC.hdf5

Reads ``MCout`` array (5 × nSamples) from HDF5, applies 3σ outlier
clipping, and produces a KDE-based pair plot with marginals.
"""

import sys
import os
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import arviz as az

# Observable labels matching MCout row order from doMC.py
LABELS = [
    r'$T_{1064}$ [ppm]',
    r'$T_{532}$ [%]',
    r'$S_\mathrm{TO}$ [$\times 10^{-21}$ m/$\sqrt{\mathrm{Hz}}$]',
    r'$S_\mathrm{Br}$ [$\times 10^{-21}$ m/$\sqrt{\mathrm{Hz}}$]',
    r'$E_\mathrm{surface}$ [V/m]',
]

VAR_NAMES = ['T1064_ppm', 'T532_pct', 'S_TO', 'S_Br', 'E_surf']


def sigma_clip(samples, nsigma=3):
    """Remove samples beyond nsigma from the median in any dimension."""
    mask = np.ones(samples.shape[1], dtype=bool)
    for i in range(samples.shape[0]):
        med = np.median(samples[i])
        sig = np.std(samples[i])
        if sig > 0:
            mask &= np.abs(samples[i] - med) < nsigma * sig
    return samples[:, mask], mask


def make_corner(hdf5_path, output_path=None, mirror_type=None):
    """Generate a corner plot from MC HDF5 output.

    Parameters
    ----------
    hdf5_path : str
        Path to MC output HDF5 (must contain ``MCout`` dataset).
    output_path : str, optional
        Where to save the figure. If None, auto-determined from
        hdf5_path and mirror_type.
    mirror_type : str, optional
        'ETM' or 'ITM'. If None, inferred from hdf5_path.

    Returns
    -------
    output_path : str
        Path to the saved figure.
    """
    if mirror_type is None:
        if 'ITM' in hdf5_path:
            mirror_type = 'ITM'
        else:
            mirror_type = 'ETM'

    if output_path is None:
        fdir = f'Figures/{mirror_type}/'
        os.makedirs(fdir, exist_ok=True)
        output_path = f'{fdir}{mirror_type}_nominal_cornerPlt.svg'

    # Load MC samples
    with h5py.File(hdf5_path, 'r') as f:
        samples = np.array(f['MCout'][:])  # (5, nSamples)
    n_obs, n_samples = samples.shape
    print(f"Loaded {n_samples} MC samples ({n_obs} observables)")

    # 3σ outlier clipping
    samples_clean, mask = sigma_clip(samples, nsigma=3)
    n_clipped = (~mask).sum()
    if n_clipped > 0:
        print(f"Clipped {n_clipped} outliers ({100*n_clipped/n_samples:.1f}%)")

    # Build ArviZ InferenceData
    names = VAR_NAMES[:n_obs]
    posterior = {name: samples_clean[i] for i, name in enumerate(names)}
    idata = az.from_dict(posterior=posterior)

    # Style
    if 'gvELOG' in plt.style.available:
        plt.style.use('gvELOG')
    else:
        plt.style.use('bmh')

    # Make the pair plot
    ax = az.plot_pair(
        idata,
        kind='kde',
        marginals=True,
        point_estimate='median',
        figsize=(18, 18),
        textsize=14,
        kde_kwargs={'contourf_kwargs': {'cmap': 'Reds'},
                    'contour_kwargs': {'colors': 'firebrick'}},
        marginal_kwargs={'color': 'firebrick'},
        point_estimate_kwargs={'color': 'k'},
        point_estimate_marker_kwargs={'marker': '+', 's': 144, 'linewidths': 2},
    )

    # Fix axis labels to use LaTeX versions
    fig = plt.gcf()
    labels = LABELS[:n_obs]
    axes = fig.get_axes()
    n = n_obs

    # Bottom row: x-labels
    for j in range(n):
        idx = (n - 1) * n + j  # bottom-left corner of grid
        if idx < len(axes):
            axes[idx].set_xlabel(labels[j], fontsize=14, fontweight='bold')
    # Left column: y-labels
    for i in range(n):
        idx = i * n
        if idx < len(axes):
            axes[idx].set_ylabel(labels[i], fontsize=14, fontweight='bold')

    fig.suptitle(f'{mirror_type} Monte Carlo Sensitivity ({samples_clean.shape[1]} samples)',
                 fontsize=18, y=0.99)
    plt.savefig(output_path, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f"Saved corner plot: {output_path}")
    return output_path


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python cornerPlt.py <MC_output.hdf5>")
        sys.exit(1)
    make_corner(sys.argv[1])
