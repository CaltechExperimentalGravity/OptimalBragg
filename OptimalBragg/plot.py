import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import corner
from physunits.frequency import Hz
from physunits.angle import deg, rad
from physunits.length import nm, um
from physunits.other import ppm

from matplotlib.ticker import FormatStrFormatter
from .layers import *


def plot_layers(stack):
    """Stack layers and normalized E

    Args:
        stack (dict): Stack attributes
    """
    lam_ref = stack["lam_ref"]
    ns, Ls, alphas = stack["ns"], stack["Ls"], stack["alphas"]
    films = list(stack["thin_films"].keys())
    film_colors = dict(
        {
            film: colors.to_hex(plt.cm.tab20b(10 * i))
            for i, film in enumerate(films)
        }
    )
    pattern = stack["pattern"]

    # Normalized z-dependent field magnitude and absorption
    N_points = 2**6
    z_arr, Enorm = field_zmag(ns, Ls, n_pts=N_points, lam=lam_ref)
    intAbs = calc_abs(Enorm, Ls, alphas)
    z_layers = np.append(0, np.cumsum(Ls / um))

    # Make the plot of the Layer structure
    fig, ax = plt.subplots(2, 1, sharex=True)
    ax[0].plot(
        z_arr / um,
        Enorm,
        color="xkcd:electric purple",
        alpha=0.97,
        rasterized=False,
    )
    absStr = Rf"$|E_s/E_{{\rm inc}}|^2 = {Enorm[0]/ppm:.0f}$ ppm"
    absStr += "\n"
    absStr += f"Absorption = {intAbs/ppm:.1f} ppm"
    ax[0].text(0.5, 0.7, absStr, transform=ax[0].transAxes, fontsize=14)

    # Add some vlines
    ax[0].vlines(
        np.cumsum(Ls / um),
        10 * um,
        0.55,
        color="xkcd:poop",
        linewidth=0.6,
        linestyle="--",
        alpha=0.75,
        rasterized=False,
    )

    # Also visualize the layer thicknesses
    for i, layer in enumerate(pattern):
        ax[1].bar(
            z_layers[i],
            Ls[i] / nm,
            width=Ls[i] / um,
            align="edge",
            color=film_colors[layer],
            alpha=0.4,
            label=stack["thin_films"][layer].Name,
        )

    ax[1].fill_betweenx(
        [0, 1.3 * Ls.max() / nm],
        z_layers[-1],
        1.5 * z_layers[-1],
        color="purple",
        alpha=0.7,
        label=stack["sub"].Name,
    )
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax[1].legend(by_label.values(), by_label.keys())
    ax[1].set_xlim(-0.5, 1.1 * z_layers[-1])
    ax[0].set_ylabel(R"Normalized $|E(z)|^2$")
    ax[1].set_ylabel(R"Physical layer thickness [nm]")
    ax[1].set_xlabel(R"Distance from air interface, $[\mu \mathrm{m}]$")
    fig.subplots_adjust(hspace=0.01, left=0.09, right=0.95, top=0.92)
    materials = list([mat.Name for mat in stack["thin_films"].values()])
    fig.suptitle(
        str([mat + ":" for mat in materials]) + " coating electric field"
    )


def plot_spectral(wavelengths, stack, markers=[], **multilayer_diel_pars):
    """Spectral transmission and reflection plots

    Args:
        wavelengths (arr): Wavelengths at which to evaluate
                           spectral reflection/transmission [m]
        stack (dict): Stack dictionary containing at least n,L keys
        markers (list, optional): Wavelengths of interest [m]
        **multilayer_diel_pars: Kwargs for multilayer_diel_pars

    """
    ns, Ls = stack["ns"], stack["Ls"]
    design_Rs = np.zeros(len(markers) + 1)
    markers.append(stack["lam_ref"])
    for ii, lam in enumerate(markers):
        design_Rs[ii] = refl(lam, stack, **multilayer_diel_pars)
    design_Ts = 1 - design_Rs

    # Spectral reflection and transmission
    RR = refl(wavelengths, stack, **multilayer_diel_pars)
    TT = trans(wavelengths, stack, **multilayer_diel_pars)
    # lam1, lam2 = qwbandedges(stack)

    fig, ax = plt.subplots(1, 1)
    # HR bandwidth wavelengths, assuming qw stack:
    # ax.axvline((lam1 / um), label="lam1", ls="--")
    # ax.axvline((lam2 / um), label="lam2", ls="--")

    ax.semilogy(
        wavelengths / um,
        TT,
        lw=1.5,
        label="Transmission",
        c="xkcd:Red",
    )
    ax.semilogy(
        wavelengths / um,
        RR,
        lw=1.5,
        label="Reflection",
        c="xkcd:electric blue",
        alpha=0.7,
    )

    colors = plt.cm.Spectral_r(np.linspace(0.1, 2, len(markers)))
    for lam, Tlam, c in zip(markers, design_Ts, colors):
        ax.vlines(
            lam / um,
            Tlam,
            1.0,
            linestyle="--",
            color=c,
            label=f"T={Tlam/ppm:.2f} ppm @ {lam/um:.3f} um",
        )
    ax.set_xlabel(R"Wavelength [$\mu \mathrm{m}$]")
    ax.set_ylabel(R"T or R")
    ax.set_ylim((1 * ppm, 1.0))
    ax.legend(loc="lower left")


def plot_noises(ff, noiselabels, plot_total=False):
    """Coating and substrate thermal noises

    Args:
        ff (array): Fourier frequency
        noiselabels (dict): Precomputed PSDs
        plot_total (bool, optional): Description
        noises (list) = Psd arrays
    """
    total = np.zeros_like(ff)
    fig, ax = plt.subplots(1, 1)
    for label, noise in noiselabels.items():
        str_eval = Rf"={np.sqrt(noise[np.argmin(np.abs(ff-100*Hz))])/1e-22:.2f}e-22 m/rtHz @ 100 Hz"
        ax.loglog(ff, np.sqrt(noise), label=label + str_eval)
        total += noise
    str_eval = Rf"={np.sqrt(total[np.argmin(np.abs(ff-100*Hz))])/1e-22:.2f}e-22 m/rtHz @ 100 Hz"
    if plot_total:
        ax.loglog(
            ff, np.sqrt(total), label="Total" + str_eval, ls="--", lw=3, c="k"
        )
    ax.legend()
    ax.set_xlim([ff[0], ff[-1]])
    ax.set_ylim([8e-24, 2e-20])
    ax.set_ylabel(R"Displacement Noise $[\mathrm{m} / \sqrt{\mathrm{Hz}}]$")
    ax.set_xlabel(R"Frequency [Hz]")


def plot_corner(mc_samples):
    """Corner plot following MCMC analysis

    Args:
        mc_samples (dict): Precomputed MC samples (see analysis)
    """
    fig, ax = plt.subplots(
        np.shape(mc_samples)[0], np.shape(mc_samples)[0], figsize=(18, 18)
    )
    corner.corner(
        mc_samples.T,
        labels=[
            "$\\mathrm{T}_{2050\\mathrm{nm}}$ [ppm]",
            "$\\mathrm{S}_{\\mathrm{TO}} [\\times 10^{-21} \\mathrm{m}/\\sqrt{\\mathrm{Hz}}]$",
            "$\\mathrm{S}_{\\mathrm{Br}} [\\times 10^{-21} \\mathrm{m}/\\sqrt{\\mathrm{Hz}}]$",
            # "$\\vec{E}_{\\mathrm{Surface}}$ [V/m]",
            "${\\alpha}_T$ [ppm]",
        ],
        # quantiles=[0.9, 0.95, 0.98],
        show_titles=True,
        use_math_text=True,
        bins=50,
        range=[
            (
                mc_samples[0].mean() - 3 * mc_samples[0].std(),
                mc_samples[0].mean() + 3 * mc_samples[0].std(),
            ),
            (
                mc_samples[1].mean() - 3 * mc_samples[1].std(),
                mc_samples[1].mean() + 3 * mc_samples[1].std(),
            ),
            (
                mc_samples[2].mean() - 3 * mc_samples[2].std(),
                mc_samples[2].mean() + 3 * mc_samples[2].std(),
            ),
            (
                mc_samples[3].mean() - 3 * mc_samples[3].std(),
                mc_samples[3].mean() + 3 * mc_samples[3].std(),
            ),
            # (
            #     mc_samples[4].mean() - 3 * mc_samples[4].std(),
            #     mc_samples[4].mean() + 3 * mc_samples[4].std(),
            # ),
        ],
        levels=(0.95,),
        color="firebrick",
        hist_kwargs={"linewidth": 2, "alpha": 0.7, "rasterized": True},
        # hist2d_kwargs={"rasterized": True},
        label_kwargs={"fontsize": 16, "fontweight": "bold"},
        title_kwargs={"fontsize": 16, "fontweight": "bold"},
        fig=fig,
    )
    ax[1, 1].text(
        0.1,
        0.5,
        "$N_s=5000$ \n $\\sigma_n^L=0.01$ \n $\\sigma_n^H=0.01$ \n $\\sigma_L=0.01$",
        fontsize=16,
        # transform=ax[1, 3].transAxes,
    )
    # ax[4, 2].xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    # ax[2, 0].yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    # for aa in ax.flatten():
    #     aa.grid(True, which="both", alpha=0.6)
    #     aa.tick_params(labelsize=56)


def plot_score(scores):
    """Polar projection artist for vector cost;
    e.g. rates output of an optimization.

    Args:
        scores (dict): Score for all costs
    """
    labels = list(scores.keys())
    r = np.array(list(scores.values()))
    theta = np.linspace(0, 2 * np.pi, np.size(r) + 1)
    normalized_list = np.linspace(0, 1, np.size(r))

    # Create figure object
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="polar")
    for i, cost in enumerate(r):
        color = plt.cm.Spectral_r(normalized_list[i])
        plt.polar(
            theta[i],
            cost,
            marker="o",
            c=color,
            markersize=8,
            markeredgewidth=0.8,
            markeredgecolor="k",
        )
        ax.grid(True, ls="--")
    ax.set_thetagrids(theta[:-1] / deg, labels, c="k")
    ax.fill(theta[:-1], r, color="goldenrod", alpha=0.2)

    scale = r.max()
    ax.set_ylim(0, scale)

    gridlines = ax.yaxis.get_gridlines()
    for gl in gridlines:
        gl.get_path()._interpolation_steps = np.size(r)
    glob_score = r.mean()

    ax.set_title(rf"Score: {glob_score*100:.2f} %")
