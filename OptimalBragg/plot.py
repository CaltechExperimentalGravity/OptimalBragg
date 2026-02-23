"""Visualization functions for coating designs.

Provides plotting of layer structure / E-field, spectral reflectivity,
thermal noise budget, starfish cost chart, and ArviZ corner plots.

All functions accept a stack dict and/or HDF5 data — no gwinc dependency.
"""

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter


# ── Style setup ──────────────────────────────────────────────────────

def apply_style():
    """Apply the project plotting style."""
    plt.style.use('bmh')

    plt.rcParams.update({
        'text.usetex': False,
        'lines.linewidth': 3,
        'font.size': 22,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.labelsize': 'medium',
        'ytick.labelsize': 'medium',
        'axes.labelsize': 'small',
        'axes.titlesize': 'medium',
        'axes.grid.axis': 'both',
        'axes.grid.which': 'both',
        'axes.grid': True,
        'grid.color': 'xkcd:Cerulean',
        'grid.alpha': 0.2,
        'lines.markersize': 12,
        'legend.borderpad': 0.2,
        'legend.fancybox': True,
        'legend.fontsize': 'small',
        'legend.framealpha': 0.8,
        'legend.handletextpad': 0.5,
        'legend.labelspacing': 0.33,
        'legend.loc': 'best',
        'figure.figsize': (12, 8),
        'savefig.dpi': 140,
        'savefig.bbox': 'tight',
        'pdf.compression': 9,
    })


# ── Layer structure & E-field ────────────────────────────────────────

def plot_layers(n, L_opt, wavelength, name_high='H', name_low='L',
                save_path=None):
    """Plot layer structure and E-field profile.

    Parameters
    ----------
    n : array_like
        Refractive indices [superstrate, layers..., substrate].
    L_opt : array_like
        Optical thicknesses (fraction of lambda).
    wavelength : float
        Design wavelength [m].
    name_high, name_low : str
        Display names for high-n and low-n materials.
    save_path : str, optional
        If given, save figure to this path.

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    from OptimalBragg.layers import op2phys, field_zmag, multidiel1, calc_abs

    L_phys = wavelength * op2phys(L_opt, n[1:-1])

    # E-field profile
    Z, field = field_zmag(n, L_phys, lam=wavelength, pol='p', n_pts=300)

    # Absorption estimate
    alpha_low, alpha_high = 0.5, 5.0
    alphas = np.where(np.arange(len(L_phys)) % 2 == 0, alpha_low, alpha_high)
    intAbs = calc_abs(field[:len(L_phys) * 300 + 1],
                      L_phys, alphas) if hasattr(calc_abs, '__call__') else 0

    layers = np.cumsum(1e6 * L_phys)
    layers = np.append(0, layers)

    fig, ax = plt.subplots(2, 1, sharex=True)

    # E-field
    ax[0].plot(Z * 1e6, field, color='xkcd:electric purple', alpha=0.97)
    absStr = (Rf"$|\vec E_{{\mathrm{{surface}}}}| = {1e6*field[0]:.0f}$ ppm "
              Rf"of $|\vec E_{{\mathrm{{inc}}}}|$")
    ax[0].text(0.5, 0.7, absStr, transform=ax[0].transAxes, fontsize=14)

    # Vertical lines at layer boundaries
    ax[0].vlines(np.cumsum(L_phys)[1:-1:2] * 1e6, 1e-5, 0.55,
                 color='xkcd:bright teal', linewidth=0.6, linestyle='--',
                 alpha=0.75)
    ax[0].vlines(np.cumsum(L_phys)[::2] * 1e6, 1e-5, 0.55,
                 color='xkcd:deep purple', linewidth=0.6, linestyle='--',
                 alpha=0.75)

    # Layer thicknesses
    ax[1].bar(layers[:-1:2], 1e9 * L_phys[::2], width=1e6 * L_phys[::2],
              align='edge', color='xkcd:bright teal', alpha=0.4,
              label=name_low)
    ax[1].bar(layers[1:-1:2], 1e9 * L_phys[1::2], width=1e6 * L_phys[1::2],
              align='edge', color='xkcd:deep purple', alpha=0.4,
              label=name_high)
    ax[1].legend()
    ax[1].yaxis.set_major_formatter(FormatStrFormatter('%3d'))
    ax[0].set_ylabel(R"Normalized $|E(z)|^2$")
    ax[1].set_ylabel(R"Physical layer thickness [nm]")
    ax[1].set_xlabel(R"Distance from air interface, $[\mu \mathrm{m}]$")

    fig.subplots_adjust(hspace=0.01, left=0.09, right=0.95, top=0.92)
    fig.suptitle(f'{name_high}:{name_low} coating electric field')

    if save_path:
        plt.savefig(save_path)
    return fig


# ── Spectral reflectivity ────────────────────────────────────────────

def plot_spectral(n, L_opt, wavelength, T1064=None, T532=None,
                  lambda_aux=None, save_path=None):
    """Plot spectral transmission/reflectivity.

    Parameters
    ----------
    n : array_like
        Refractive indices.
    L_opt : array_like
        Optical thicknesses.
    wavelength : float
        Design wavelength [m].
    T1064, T532 : float, optional
        Transmission values to annotate.
    lambda_aux : float, optional
        AUX wavelength ratio (e.g. 0.5 for 532 nm / 1064 nm).
        If given, the plot range extends to cover the AUX wavelength.
    save_path : str, optional
        If given, save figure.

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    from OptimalBragg.layers import multidiel1

    # Extend range to cover AUX wavelength if provided
    lo = 0.75
    if lambda_aux and lambda_aux < lo:
        lo = lambda_aux * 0.9
    wavelengths = np.linspace(lo, 1.25, 512)
    rr, _ = multidiel1(n, L_opt, wavelengths)
    RR = np.abs(rr) ** 2
    TT = 1 - RR

    fig, ax = plt.subplots(1, 1)
    wl_um = 1e6 * wavelengths * wavelength
    ax.semilogy(wl_um, TT, lw=1.5,
                label='Transmissivity', c='xkcd:Red')
    ax.semilogy(wl_um, RR, lw=1.5,
                label='Reflectivity', c='xkcd:electric blue', alpha=0.7)

    if T1064 and T1064 > 0:
        ax.axvline(1e6 * wavelength, ls='--', color='blue', alpha=0.5,
                   label=f'T={T1064*1e6:.2f} ppm @ {1e6*wavelength:.3f} um')

    if lambda_aux:
        ax.axvline(1e6 * lambda_aux * wavelength, ls=':', color='green',
                   alpha=0.6,
                   label=f'AUX @ {1e6*lambda_aux*wavelength:.3f} um')

    ax.set_xlabel(R"Wavelength [$\mu \mathrm{m}$]")
    ax.set_ylabel(R"T or R")
    ax.set_ylim((5e-8, 1.0))
    ax.legend(loc='lower left')

    if save_path:
        plt.savefig(save_path)
    return fig


# ── Thermal noise budget ─────────────────────────────────────────────

def plot_noise(stack, w_beam, r_mirror=0.17, d_mirror=0.20,
               freq=None, save_path=None):
    """Plot thermal noise budget using OptimalBragg noise models.

    Parameters
    ----------
    stack : dict
        Stack dict (with optimized ``Ls_opt``).
    w_beam : float
        Beam radius [m].
    r_mirror, d_mirror : float
        Mirror radius and thickness [m].
    freq : array_like, optional
        Frequency array [Hz]. Default ``logspace(0, 4, 500)``.
    save_path : str, optional
        If given, save figure.

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    from OptimalBragg.noise import (
        coating_thermooptic, coating_brownian,
        substrate_brownian, substrate_thermoelastic,
    )

    if freq is None:
        freq = np.logspace(0, 4, 500)

    StoZ, SteZ, StrZ = coating_thermooptic(
        freq, stack, w_beam, r_mirror, d_mirror)
    SbrZ = coating_brownian(freq, stack, w_beam)
    subBrown = substrate_brownian(freq, stack, w_beam, r_mirror, d_mirror)
    subTE = substrate_thermoelastic(freq, stack, w_beam, r_mirror, d_mirror)

    CTNtot = np.sqrt(StoZ + SbrZ)

    fig, ax = plt.subplots(1, 1)
    ax.loglog(freq, np.sqrt(np.abs(StoZ)), label='Thermo-Optic',
              c='xkcd:Purplish Blue')
    ax.loglog(freq, np.sqrt(SteZ), label='Thermo-Elastic',
              c='xkcd:Golden')
    ax.loglog(freq, np.sqrt(StrZ), label='Thermo-Refractive',
              c='xkcd:Puke')

    idx100 = np.argmin(np.abs(freq - 100))
    ax.loglog(freq, np.sqrt(SbrZ),
              label=f'Brownian={np.sqrt(SbrZ[idx100])/1e-22:.2f}e-22 @ 100 Hz',
              c='xkcd:Tomato')
    ax.loglog(freq, np.sqrt(subBrown), label='Substrate Brownian',
              c='xkcd:Dusty Blue')
    ax.loglog(freq, np.sqrt(subTE), label='Substrate Thermo-Elastic',
              c='xkcd:Chocolate', alpha=0.3)
    ax.loglog(freq, CTNtot, c='k', lw=3, ls='--',
              label=Rf'CTN total = {CTNtot[idx100]/1e-22:.2f}e-22 @ 100 Hz')

    n_layers = len(stack["Ls"])
    thickness_um = np.sum(stack["Ls"]) * 1e6
    ax.text(80, 11e-21, f'# of layers = {n_layers}', size='x-small')
    ax.text(80, 5e-21, f'Thickness = {thickness_um:.2f} um', size='x-small')

    ax.legend()
    ax.set_xlim([10, 10e3])
    ax.set_ylim([8e-24, 2e-20])
    ax.set_ylabel(R"Displacement Noise $[\mathrm{m} / \sqrt{\mathrm{Hz}}]$")
    ax.set_xlabel(R"Frequency [Hz]")

    if save_path:
        plt.savefig(save_path)
    return fig


# ── Starfish (polar cost chart) ──────────────────────────────────────

def plot_starfish(scalar_costs, scale=None, save_path=None, title=''):
    """Polar projection chart for cost function terms.

    Parameters
    ----------
    scalar_costs : dict
        ``{cost_name: value}`` for each active cost term.
    scale : float, optional
        Radial axis limit. Default: ``max(values)``.
    save_path : str, optional
        If given, save figure.
    title : str, optional
        Figure title.

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    labels = list(scalar_costs.keys())
    r = np.array(list(scalar_costs.values()), dtype=np.float64)

    theta = np.linspace(0, 2 * np.pi, len(r) + 1)
    normalized_list = np.linspace(0, 1, len(r))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='polar')
    for i, cost in enumerate(r):
        color = plt.cm.Spectral_r(normalized_list[i])
        ax.plot(theta[i], cost, marker='o', c=color, markersize=8,
                markeredgewidth=0.8, markeredgecolor='k')
    ax.grid(True, ls='--')
    ax.set_thetagrids(np.degrees(theta[:-1]), labels, c='k')
    ax.fill(theta[:-1], r, color='goldenrod', alpha=0.2)

    if scale is None:
        scale = r.max()
    ax.set_ylim(0, scale)

    gridlines = ax.yaxis.get_gridlines()
    for gl in gridlines:
        gl.get_path()._interpolation_steps = len(r)

    ax.set_title(title)

    if save_path:
        plt.savefig(save_path, transparent=True, dpi=200)
    return fig


# ── Corner plot ──────────────────────────────────────────────────────

# Observable labels matching MCout row order
CORNER_LABELS = [
    r'$T_{1064}$ [ppm]',
    r'$T_{532}$ [%]',
    r'$S_\mathrm{TO}$ [$\times 10^{-21}$ m/$\sqrt{\mathrm{Hz}}$]',
    r'$S_\mathrm{Br}$ [$\times 10^{-21}$ m/$\sqrt{\mathrm{Hz}}$]',
    r'$E_\mathrm{surface}$ [V/m]',
]

CORNER_VAR_NAMES = ['T1064_ppm', 'T532_pct', 'S_TO', 'S_Br', 'E_surf']


def sigma_clip(samples, nsigma=3):
    """Remove samples beyond nsigma from the median in any dimension."""
    mask = np.ones(samples.shape[1], dtype=bool)
    for i in range(samples.shape[0]):
        med = np.median(samples[i])
        sig = np.std(samples[i])
        if sig > 0:
            mask &= np.abs(samples[i] - med) < nsigma * sig
    return samples[:, mask], mask


def plot_corner(mc_samples, mirror_type='ETM', save_path=None):
    """Generate ArviZ corner plot from MC samples.

    Parameters
    ----------
    mc_samples : ndarray
        Shape ``(n_observables, n_samples)`` — the MCout array.
    mirror_type : str
        'ETM' or 'ITM' for title.
    save_path : str, optional
        If given, save figure.

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    import arviz as az

    matplotlib.use('Agg')

    n_obs, n_samples = mc_samples.shape

    # 3-sigma clipping
    samples_clean, mask = sigma_clip(mc_samples, nsigma=3)
    n_clipped = (~mask).sum()
    if n_clipped > 0:
        print(f"Clipped {n_clipped} outliers "
              f"({100*n_clipped/n_samples:.1f}%)")

    # Filter out degenerate dimensions (zero variance)
    keep = []
    for i in range(n_obs):
        if np.std(samples_clean[i]) > 1e-15 * np.abs(np.mean(samples_clean[i])):
            keep.append(i)
        else:
            print(f"Skipping degenerate dimension {i} "
                  f"({CORNER_VAR_NAMES[i]}): zero variance")
    samples_clean = samples_clean[keep]
    used_names = [CORNER_VAR_NAMES[i] for i in keep]
    used_labels = [CORNER_LABELS[i] for i in keep]
    n_plot = len(keep)

    # Build ArviZ InferenceData
    posterior = {name: samples_clean[i] for i, name in enumerate(used_names)}
    idata = az.from_dict(posterior=posterior)

    apply_style()

    ax = az.plot_pair(
        idata, kind='kde', marginals=True, point_estimate='median',
        figsize=(18, 18), textsize=14,
        kde_kwargs={'contourf_kwargs': {'cmap': 'Reds'},
                    'contour_kwargs': {'colors': 'firebrick'}},
        marginal_kwargs={'color': 'firebrick'},
        point_estimate_kwargs={'color': 'k'},
        point_estimate_marker_kwargs={'marker': '+', 's': 144,
                                      'linewidths': 2},
    )

    fig = plt.gcf()
    labels = used_labels
    axes = fig.get_axes()

    # Bottom row x-labels, left column y-labels
    for j in range(n_plot):
        idx = (n_plot - 1) * n_plot + j
        if idx < len(axes):
            axes[idx].set_xlabel(labels[j], fontsize=14, fontweight='bold')
    for i in range(n_plot):
        idx = i * n_plot
        if idx < len(axes):
            axes[idx].set_ylabel(labels[i], fontsize=14, fontweight='bold')

    fig.suptitle(f'{mirror_type} Monte Carlo Sensitivity '
                 f'({samples_clean.shape[1]} samples)',
                 fontsize=18, y=0.99)

    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=150)
        plt.close(fig)
    return fig
