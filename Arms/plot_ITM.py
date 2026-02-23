""" Python script that plots the layer thicknesses and E field (normalized)
Example usage:
    plot_ITM.py
will take the most recent coating design HDF5 and make plots of the
E-field within the dielectric layer structure for aLIGO SiO2/Ta2O5."""

import sys, glob, os
import subprocess
import h5py

from gwinc import noise, Struct

from generic.coatingUtils import *
from generic.reportUtils import generate_run_rst

import numpy as np
import matplotlib.pyplot as plt
from starfish import polar_cost
from matplotlib.ticker import FormatStrFormatter


# setup paths for data and figs; make dirs if do not exist yet
spath = "Data/ITM/"
os.makedirs(spath, exist_ok=True)

fpath = "Figures/ITM/"
os.makedirs(fpath, exist_ok=True)

savePlots = True

_positional = [a for a in sys.argv[1:] if not a.startswith('--')]
if not _positional:
    fname = max(glob.iglob(spath + "*Layers*.hdf5"), key=os.path.getctime)
else:
    fname = str(_positional[0])


def h5read(targets):
    data = {}
    with h5py.File(fname, "r+") as f:
        for target in targets:
            try:
                data[target] = np.array(f["diffevo_output/" + target])
            except:
                data[target] = 0.0
    return data


def plot_layers(save=savePlots):
    # Load ifo
    opt_params = importParams("ITM_params.yml")
    ifo = Struct.from_file(opt_params["misc"]["gwincStructFile"])
    lambdaPSL = ifo.Laser.Wavelength

    data = h5read(targets=["n", "L"])
    L = lambdaPSL * op2phys(data["L"], data["n"][1:-1])

    # Absorption coefficients [1/m] for IBS SiO2/Ta2O5 at 1064 nm
    # Order-of-magnitude estimates consistent with total stack absorption
    # ~0.2-0.5 ppm (Granata et al. CQG 37 095004, 2020; Pinard et al.
    # Appl. Opt. 56 C11, 2017)
    alpha_low = 0.5
    alpha_high = 5.0

    Name_high = ifo.Materials.Coating.Name_high
    Name_low = ifo.Materials.Coating.Name_low
    # Calculate longitudinal field amplitude
    Z, field = fieldDepth(
        L, data["n"], pol="p", nPts=300, lam=ifo.Laser.Wavelength
    )
    intAbs = calcAbsorption(field, L, 300, alpha_low, alpha_high)
    layers = np.cumsum(1e6 * L)
    layers = np.append(0, layers)

    # Make the plot of the Layer structure
    fig2, ax2 = plt.subplots(2, 1, sharex=True)
    ax2[0].plot(
        Z * 1e6,
        field,
        color="xkcd:electric purple",
        alpha=0.97,
        rasterized=False,
    )
    absStr = Rf"$|\vec E_{{\mathrm{{surface}}}}| = {1e6*field[0]:.0f}$ ppm of $\vec |E_{{\mathrm{{inc}}}}|$"
    absStr += "\n"
    absStr += f"Integrated absorption in stack is {intAbs:.1f} ppm"
    ax2[0].text(0.5, 0.7, absStr, transform=ax2[0].transAxes, fontsize=14)

    # Add some vlines
    ax2[0].vlines(
        np.cumsum(L)[1:-1:2] * 1e6,
        1e-5,
        0.55,
        color="xkcd:bright teal",
        linewidth=0.6,
        linestyle="--",
        alpha=0.75,
        rasterized=False,
    )
    ax2[0].vlines(
        np.cumsum(L)[::2] * 1e6,
        1e-5,
        0.55,
        color="xkcd:deep purple",
        linewidth=0.6,
        linestyle="--",
        alpha=0.75,
        rasterized=False,
    )

    # Also visualize the layer thicknesses
    ax2[1].bar(
        layers[:-1:2],
        1e9 * L[::2],
        width=1e6 * L[::2],
        align="edge",
        color="xkcd:bright teal",
        alpha=0.4,
        label=Name_low,
    )
    ax2[1].bar(
        layers[1:-1:2],
        1e9 * L[1::2],
        width=1e6 * L[1::2],
        align="edge",
        color="xkcd:deep purple",
        alpha=0.4,
        label=Name_high,
    )
    ax2[1].legend()
    ax2[1].yaxis.set_major_formatter(FormatStrFormatter("%3d"))
    ax2[0].set_ylabel(R"Normalized $|E(z)|^2$")
    ax2[1].set_ylabel(R"Physical layer thickness [nm]")
    ax2[1].set_xlabel(R"Distance from air interface, $[\mu \mathrm{m}]$")

    fig2.subplots_adjust(hspace=0.01, left=0.09, right=0.95, top=0.92)
    fig2.suptitle(Name_high + ":" + Name_low + " coating electric field")

    plt.savefig(fpath + "/ITM_Layers_" + fname[-16:-5] + ".svg")


def plot_trans(save=savePlots):
    # Load ifo
    opt_params = importParams("ITM_params.yml")
    ifo = Struct.from_file(opt_params["misc"]["gwincStructFile"])
    lambdaPSL = ifo.Laser.Wavelength

    data = h5read(targets=["n", "L", "T1064", "T532", "TOPL"])
    L = lambdaPSL * op2phys(data["L"], data["n"][1:-1])

    # Spectral reflectivities
    T1064 = float(data["T1064"])
    try:
        T532 = float(data["T532"])
    except (TypeError, KeyError):
        T532 = 0.0
    try:
        TOPV = float(data["TOPL"])
    except (KeyError, TypeError):
        TOPV = 0.0

    wavelengths = np.linspace(0.75, 1.25, 512)
    rr, _ = multidiel1(data["n"], data["L"], wavelengths)
    RR = np.abs(rr) ** 2
    TT = 1 - RR

    # Build stats based on various figures for star fish chart
    stats = {}
    for cost, spec in opt_params["costs"].items():
        if spec["weight"]:
            optcost = list(h5read(targets=["vectorCost/" + cost]).values())[0]
            stat = 1 / (1e-11 + optcost)
            if optcost > 1e3 or optcost < 1e-5:
                stat = 1e-2
            stats[cost] = np.abs(np.log(np.abs(stat)))

    Nfixed = opt_params["misc"]["Nfixed"]
    Nlayers = 2 * opt_params["misc"]["Npairs"] + 1
    N_particles = opt_params["misc"]["Nparticles"]

    sfplot_head = fpath + "ITM_SF"
    sfplot_tail = ".png"
    sfplot_name = sfplot_head + fname[-16:-5] + sfplot_tail
    polar_cost(
        stats,
        scale=10,
        fname=sfplot_name,
        figtitle=Rf"Stack with {Nlayers} layers ({Nfixed} fixed bilayers) and {N_particles} particles",
    )
    # Convert from optical thickness to physical thickness
    L = lambdaPSL * op2phys(data["L"], data["n"][1:-1])

    fig, ax = plt.subplots(1, 1)
    ax.semilogy(
        1e6 * wavelengths * lambdaPSL,
        TT,
        lw=1.5,
        label="Transmissivity",
        c="xkcd:Red",
    )
    ax.semilogy(
        1e6 * wavelengths * lambdaPSL,
        RR,
        lw=1.5,
        label="Reflectivity",
        c="xkcd:electric blue",
        alpha=0.7,
    )
    for wvl, trans, c in zip(
        [1.0],
        [T1064],
        ["blue"],
    ):
        if trans:
            ax.vlines(
                wvl * 1e6 * lambdaPSL,
                trans,
                1.0,
                linestyle="--",
                color=c,
                label=f"T={trans*1e6:.2f} ppm @ {1e6*wvl*lambdaPSL:.3f} um",
            )
    ax.set_xlabel(R"Wavelength [$\mu \mathrm{m}$]")
    ax.set_ylabel(R"T or R")
    ax.set_ylim((5e-8, 1.0))

    # Add starfish inset
    im = plt.imread(fpath + "ITM_SF" + fname[-16:-5] + ".png")
    newax = fig.add_axes([0.55, 0.1, 0.4, 0.4], anchor="NE")
    newax.imshow(im)
    newax.axis("off")

    ax.legend(loc="lower left")
    plt.savefig(fpath + "ITM_R" + fname[-16:-5] + ".svg")


def plot_noise(save=savePlots):
    # Load ifo
    opt_params = importParams("ITM_params.yml")
    ifo = Struct.from_file(opt_params["misc"]["gwincStructFile"])
    lambdaPSL = ifo.Laser.Wavelength

    data = h5read(targets=["L", "n"])
    L = lambdaPSL * op2phys(data["L"], data["n"][1:-1])

    # Frequency band
    ff = np.logspace(0, 4, 500)

    # Build up a "mirror" structure as required by pygwinc
    # Use ETM for coating/substrate properties, ITM for beam radius
    mir = ifo.Optics.ETM
    mir.Coating.dOpt = data["L"][:]
    StoZ, SteZ, StrZ, _ = noise.coatingthermal.coating_thermooptic(
        ff, mir, ifo.Laser.Wavelength, ifo.Optics.ITM.BeamRadius
    )
    SbrZ = noise.coatingthermal.coating_brownian(
        ff, mir, ifo.Laser.Wavelength, ifo.Optics.ITM.BeamRadius
    )
    subBrown = noise.substratethermal.substrate_brownian(
        ff, mir, ifo.Optics.ITM.BeamRadius
    )
    subTE = noise.substratethermal.substrate_thermoelastic(
        ff, mir, ifo.Optics.ITM.BeamRadius
    )

    CTNtot = np.sqrt(StoZ + SbrZ)
    SUBtot = np.sqrt(subBrown + subTE)
    Stot = np.sqrt(CTNtot**2 + SUBtot**2)

    Larm = 1
    # Figure
    fig3, ax3 = plt.subplots(1, 1)
    ax3.loglog(
        ff,
        Larm * np.sqrt(np.abs(StoZ)),
        label="Thermo-Optic",
        c="xkcd:Purplish Blue",
    )
    ax3.loglog(
        ff, Larm * np.sqrt(SteZ), label="Thermo-Elastic", c="xkcd:Golden"
    )
    ax3.loglog(
        ff, Larm * np.sqrt(StrZ), label="Thermo-Refractive", c="xkcd:Puke"
    )
    ax3.loglog(
        ff,
        Larm * np.sqrt(SbrZ),
        label=f"Brownian={np.sqrt(SbrZ[np.argmin(np.abs(ff-100))])/1e-22:.2f}e-22 m/rtHz @ 100 Hz",
        c="xkcd:Tomato",
    )
    ax3.loglog(
        ff,
        Larm * np.sqrt(subBrown),
        label="Substrate Brownian",
        c="xkcd:Dusty Blue",
    )
    ax3.loglog(
        ff,
        Larm * np.sqrt(subTE),
        label="Substrate Thermo-Elastic",
        c="xkcd:Chocolate",
        alpha=0.3,
    )
    ax3.loglog(
        ff,
        Larm * CTNtot,
        c="k",
        lw=3,
        ls="--",
        label=Rf"CTN total = {CTNtot[np.argmin(np.abs(ff-100.0))]/1e-22:.2f}e-22 m/rtHz @ 100 Hz",
    )
    ax3.legend()
    ax3.set_xlim([10, 10e3])
    ax3.set_ylim([8e-24, 2e-20])
    ax3.text(80, 11e-21, "# of layers =  {}".format(len(L)), size="x-small")
    ax3.text(
        80,
        5e-21,
        "Thickness = {} um".format(round(sum(L / 1e-6), 2)),
        size="x-small",
    )
    ax3.set_ylabel(R"Displacement Noise $[\mathrm{m} / \sqrt{\mathrm{Hz}}]$")
    ax3.set_xlabel(R"Frequency [Hz]")

    plt.savefig(fpath + "ITM_TN_" + fname[-16:-5] + ".svg")


def main(layers=True, trans=True, noise=True):
    # Setup matplotlib params
    plt.rcParams.update(
        {
            "text.usetex": False,
            "lines.linewidth": 3,
            "font.size": 22,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.labelsize": "medium",
            "ytick.labelsize": "medium",
            "axes.labelsize": "small",
            "axes.titlesize": "medium",
            "axes.grid.axis": "both",
            "axes.grid.which": "both",
            "axes.grid": True,
            "grid.color": "xkcd:Cerulean",
            "grid.alpha": 0.2,
            "lines.markersize": 12,
            "legend.borderpad": 0.2,
            "legend.fancybox": True,
            "legend.fontsize": "small",
            "legend.framealpha": 0.8,
            "legend.handletextpad": 0.5,
            "legend.labelspacing": 0.33,
            "legend.loc": "best",
            "figure.figsize": ((12, 8)),
            "savefig.dpi": 140,
            "savefig.bbox": "tight",
            "pdf.compression": 9,
        }
    )

    if __debug__:
        print("Loading " + "gwinc.ifo" + " " + fname)
    plot_layers(layers)
    plot_trans(trans)
    plot_noise(noise)

    # Generate Sphinx run report
    ts = fname[-16:-5]  # YYMMDD_HHMMSS
    fig_paths = {
        'layers': f'Figures/ITM/ITM_Layers_{ts}.svg',
        'reflectivity': f'Figures/ITM/ITM_R{ts}.svg',
        'starfish': f'Figures/ITM/ITM_SF{ts}.png',
        'thermal_noise': f'Figures/ITM/ITM_TN_{ts}.svg',
    }
    opt_params = importParams("ITM_params.yml")
    ifo = Struct.from_file(opt_params["misc"]["gwincStructFile"])
    opt_params['_wavelength'] = ifo.Laser.Wavelength

    # Check for existing MC data
    mc_path = spath + 'ITM_MC.hdf5'
    mc_exists = os.path.isfile(mc_path)
    if mc_exists:
        corner_path = fpath + 'ITM_nominal_cornerPlt.svg'
        if os.path.isfile(corner_path):
            fig_paths['corner'] = 'Figures/ITM/ITM_nominal_cornerPlt.svg'

    _, rst_path = generate_run_rst(fname, 'ITM', fig_paths, opt_params,
                                   project_dir='Arms',
                                   mc_hdf5_path=mc_path if mc_exists else None)
    print(f"Run report: {rst_path}")

    # Launch MC in background (unless --no-mc flag)
    if '--no-mc' not in sys.argv:
        print("Launching MC in background (5000 samples)...")
        subprocess.Popen(
            [sys.executable, __file__, fname, '--mc-bg'],
            cwd=os.getcwd(),
            stdout=open(spath + 'mc_bg.log', 'w'),
            stderr=subprocess.STDOUT,
        )


def _run_mc_background(layers_hdf5):
    """Background job: run MC, generate corner plot, update report."""
    mc_path = spath + 'ITM_MC.hdf5'
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # 1. Run doMC_ITM as subprocess
    subprocess.run(
        [sys.executable, os.path.join(script_dir, 'doMC_ITM.py'),
         layers_hdf5, mc_path, '5000'],
        check=True, cwd=os.getcwd(),
    )

    # 2. Generate corner plot
    from cornerPlt import make_corner
    make_corner(mc_path, mirror_type='ITM')

    # 3. Regenerate report with MC data
    ts = layers_hdf5[-16:-5]
    fig_paths = {
        'layers': f'Figures/ITM/ITM_Layers_{ts}.svg',
        'reflectivity': f'Figures/ITM/ITM_R{ts}.svg',
        'starfish': f'Figures/ITM/ITM_SF{ts}.png',
        'thermal_noise': f'Figures/ITM/ITM_TN_{ts}.svg',
        'corner': 'Figures/ITM/ITM_nominal_cornerPlt.svg',
    }
    opt_params = importParams("ITM_params.yml")
    ifo = Struct.from_file(opt_params["misc"]["gwincStructFile"])
    opt_params['_wavelength'] = ifo.Laser.Wavelength
    _, rst_path = generate_run_rst(layers_hdf5, 'ITM', fig_paths, opt_params,
                                   project_dir='Arms', mc_hdf5_path=mc_path)
    print(f"MC complete. Updated report: {rst_path}")


if __name__ == "__main__":
    if '--mc-bg' in sys.argv:
        _run_mc_background(sys.argv[1])
    else:
        if __debug__:
            print("Make plots...")
        main(layers=True, trans=True, noise=False)
