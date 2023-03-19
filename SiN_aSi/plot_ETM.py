''' Python script that plots the layer thicknesses and E field (normalized)
Example usage:
	plotLayers.py aLIGO_ETM_20layers.mat
will take the coating design in aLIGO_ETM_20layers.mat and make a plot of the
E-field within the dielectric layer structure.'''

import sys, glob, os
import h5py

# To use pygwinc, first install; 
# >> conda install -c conda-forge gwinc
# import sys
# sys.path.append('../../pygwinc/')
from gwinc import noise, Struct

sys.path.append('../generic/')
from generic_local.coatingUtils import *

import numpy as np
import matplotlib.pyplot as plt
from starfish import polar_cost
from matplotlib.ticker import FormatStrFormatter


# setup paths for data and figs; make dirs if do not exist yet
spath = 'Data/ETM/'
os.makedirs(spath, exist_ok = True)

fpath = 'Figures/ETM/'
os.makedirs(fpath, exist_ok = True)

savePlots = True

if len(sys.argv) == 1:
    fname = max(glob.iglob(spath + '*Layers*.hdf5'), key=os.path.getctime)
    # fname = fname[5:] # rm 'data' from the name
else:
    # For example fname = 'ETM_Layers_190519_1459.hdf5'
    fname = str(sys.argv[1])

def h5read(targets):
    data = {}
    with h5py.File(fname, 'r+') as f:
        for target in targets:
            try:
                data[target] = np.array(f['diffevo_output/' + target])
            except:
                data[target] = 0.
    return data

def plot_layers(save = savePlots):
    # Load ifo
    opt_params = importParams('ETM_params.yml')
    ifo = Struct.from_file(opt_params['misc']['gwincStructFile'])
    lambdaPSL = ifo.Laser.Wavelength

    data = h5read(targets=['n', 'L'])
    L  = lambdaPSL * op2phys(data['L'], data['n'][1:-1])

    # Constants for calculating absorption, we should get a better 
    # number for 2.128 um at 123 K than -->
    #   https://journals.aps.org/prd/pdf/10.1103/PhysRevD.91.042002
    #   https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.120.263602#page=3
    #       -- figures (2-3)
    alpha_low = 1e-3
    alpha_high = 10e-6 / 0.5e-6  # personal comm from Manel Ruiz to RXA 3/2023

    Name_high = ifo.Materials.Coating.Name_high
    Name_low  = ifo.Materials.Coating.Name_low
    # Calculate longitudinal field amplitude
    Z, field = fieldDepth(L, data['n'], pol='p', nPts=300, lam=ifo.Laser.Wavelength)
    intAbs = calcAbsorption(field, L, 300, alpha_low, alpha_high)
    layers = np.cumsum(1e6 * L)
    layers = np.append(0, layers)

    # Make the plotof the Layer structure
    fig2 , ax2 = plt.subplots(2,1, sharex=True)
    ax2[0].plot(Z*1e6,field, color='xkcd:electric purple',
                   alpha=0.97, rasterized=False)
    absStr = fR'$|\vec E_{{\mathrm{{surface}}}}| = {1e6*field[0]:.0f}$ ppm of $\vec |E_{{\mathrm{{inc}}}}|$'
    absStr += '\n'
    absStr += f'Integrated absorption in stack is {intAbs:.1f} ppm'
    ax2[0].text(0.5, 0.7, absStr, transform=ax2[0].transAxes, fontsize=14)

    #Add some vlines
    ax2[0].vlines(np.cumsum(L)[1:-1:2]*1e6, 1e-5, 0.55,
                color='xkcd:bright teal', linewidth=0.6,
                     linestyle='--', alpha=0.75, rasterized=False)
    ax2[0].vlines(np.cumsum(L)[::2]*1e6, 1e-5, 0.55,
                color='xkcd:deep purple', linewidth=0.6,
                     linestyle='--', alpha=0.75,rasterized=False)

    #Also visualize the layer thicknesses
    ax2[1].bar(layers[:-1:2], 1e9*L[::2], width=1e6*L[::2],

            align='edge', color='xkcd:bright teal',
                  alpha=0.4, label=Name_low)
    ax2[1].bar(layers[1:-1:2],  1e9*L[1::2], width=1e6*L[1::2],
            align='edge', color='xkcd:deep purple',
                  alpha=0.4, label=Name_high)
    ax2[1].legend()
    ax2[1].yaxis.set_major_formatter(FormatStrFormatter("%3d"))
    ax2[0].set_ylabel(R'Normalized $|E(z)|^2$')
    ax2[1].set_ylabel(R'Physical layer thickness [nm]')
    ax2[1].set_xlabel(R'Distance from air interface, $[\mu \mathrm{m}]$')

    fig2.subplots_adjust(hspace=0.01,left=0.09,right=0.95,top=0.92)
    fig2.suptitle(Name_high + ':' + Name_low + ' coating electric field')

    plt.savefig(fpath + '/ETM_Layers_' + fname[-16:-5] + '.pdf')
    plt.savefig(fpath + 'ETM_Layers' + '.pdf')


def plot_trans(save = savePlots):
    # Load ifo
    opt_params = importParams('ETM_params.yml')
    ifo = Struct.from_file(opt_params['misc']['gwincStructFile'])
    lambdaPSL = ifo.Laser.Wavelength

    data = h5read(targets=['n', 'L', 'TPSL', 'TAUX', 'TOPL'])
    L  = lambdaPSL * op2phys(data['L'], data['n'][1:-1])

    # Other wavelength ratios to ifo.Laser.Wavelength
    rellambdaAUX = 2/3
    rellambda1550 = 1550e-9/lambdaPSL
    rellambdaHeNe = 632e-9/lambdaPSL

    # Spectral reflectivities
    rr1550, _ = multidiel1(data['n'], data['L'], rellambda1550)
    T1550 = 1 - np.abs(rr1550[0])**2
    TPSL = data['TPSL'][0]
    TAUX = data['TAUX'][0]
    TOPV = data['TOPL'][0]

    wavelengths = np.linspace(0.2, 1.8, 512)
    rr, _   = multidiel1(data['n'], data['L'], wavelengths)
    RR      = np.abs(rr)**2
    TT      = 1 - RR

    # Build stats based on various figures for star fish chart
    stats = {}
    for cost, spec in opt_params['costs'].items():
        if spec['weight']:
            optcost = list(h5read(targets=['vectorCost/' + cost]).values())[0]
            stat = 1 / (1e-11 + optcost)  # 1e-11 is there to avoid div by zero
            if optcost > 1e3 or optcost < 1e-5:
                stat = 1e-2
            stats[cost] = np.abs(np.log(np.abs(stat)))

    Nfixed      = opt_params['misc']['Nfixed']
    Nlayers     = 2*opt_params['misc']['Npairs'] + 1
    N_particles = opt_params['misc']['Nparticles']

    sfplot_head = fpath + 'ETM_SF'
    sfplot_tail = '.png'
    sfplot_name = sfplot_head + fname[-16:-5] + sfplot_tail
    polar_cost(stats, 
               scale = 10,
               fname = sfplot_name,
               figtitle = fR'Stack with {Nlayers} layers ({Nfixed} fixed bilayers) and {N_particles} particles',
                )
    os.system('cp ' + sfplot_name + ' ' + sfplot_head + sfplot_tail)  

    # Convert from optical thickness to physical thickness
    L  = lambdaPSL * op2phys(data['L'], data['n'][1:-1])

    fig, ax = plt.subplots(1,1)
    ax.semilogy(1e6*wavelengths*lambdaPSL, TT, lw=1.5, 
                label='Transmissivity', c='xkcd:Red')
    ax.semilogy(1e6*wavelengths*lambdaPSL, RR, lw=1.5, 
                label='Reflectivity', c='xkcd:electric blue', alpha=0.7)
    for wvl, trans, c in zip([1., rellambda1550, rellambdaAUX, rellambdaHeNe],
                             [TPSL, T1550, TAUX, TOPV],
                             ['brown', 'crimson', 'blue', 'red']):
        ax.vlines(wvl*1e6*lambdaPSL, trans, 1.0, linestyle='--', color=c,
                  label=f'T={trans*1e6:.2f} ppm @ {1e6*wvl*lambdaPSL:.3f} um')
    ax.set_xlabel(R'Wavelength [$\mu \mathrm{m}$]')
    ax.set_ylabel(R'T or R')
    ax.set_ylim((5e-8, 1.0))

    # Add starfish inset
    im = plt.imread(fpath + 'ETM_SF' + fname[-16:-5] + '.png')
    newax = fig.add_axes([0.55, 0.1, 0.4, 0.4], anchor = 'NE')
    newax.imshow(im)
    newax.axis('off')

    ax.legend(loc = 'best')
    plt.savefig(fpath + 'ETM_R' + fname[-16:-5] + '.pdf')
    plt.savefig(fpath + 'ETM_R' + '.pdf')


def plot_noise(save = savePlots):
    # Load ifo
    opt_params = importParams('ETM_params.yml')
    ifo = Struct.from_file(opt_params['misc']['gwincStructFile'])
    lambdaPSL = ifo.Laser.Wavelength

    data = h5read(targets=['L', 'n'])
    L  = lambdaPSL * op2phys(data['L'], data['n'][1:-1])

    # Frequency band
    ff = np.logspace(0, 4, 500)

    # Build up a "mirror" structure as required by pygwinc
    mir = ifo.Optics.ETM
    mir.Coating.dOpt = data['L'][:]
    StoZ, SteZ, StrZ, _ = noise.coatingthermal.coating_thermooptic(ff, 
                                mir, ifo.Laser.Wavelength, ifo.Optics.ETM.BeamRadius)
    SbrZ = noise.coatingthermal.coating_brownian(ff, 
                                mir, ifo.Laser.Wavelength, ifo.Optics.ETM.BeamRadius)
    subBrown = noise.substratethermal.substrate_brownian(ff, 
                                mir, ifo.Optics.ETM.BeamRadius)
    subTE = noise.substratethermal.substrate_thermoelastic(ff, 
                                mir, ifo.Optics.ETM.BeamRadius)

    Larm = 1 #4000

    # Figure
    fig3, ax3 = plt.subplots(1,1)
    ax3.loglog(ff, Larm * np.sqrt(StoZ), label='Thermo-Optic', c='xkcd:Purplish Blue')
    ax3.loglog(ff, Larm * np.sqrt(SteZ), label='Thermo-Elastic', c='xkcd:Golden')
    ax3.loglog(ff, Larm * np.sqrt(StrZ), label='Thermo-Refractive', c='xkcd:Puke')
    ax3.loglog(ff, Larm * np.sqrt(SbrZ), label='Brownian', c='xkcd:Tomato')
    ax3.loglog(ff, Larm * np.sqrt(subBrown), label='Substrate Brownian', c='xkcd:Dusty Blue')
    ax3.loglog(ff, Larm * np.sqrt(subTE), label='Substrate Thermo-Elastic', c='xkcd:Chocolate', alpha=0.3)
    ax3.legend()
    ax3.set_xlim([10, 10e3])
    ax3.set_ylim([8e-24, 2e-20])
    ax3.text(80, 11e-21, '# of layers =  {}'.format(len(L)), size='x-small')
    ax3.text(80, 5e-21, 'Thickness = {} um'.format(round(1e6*sum(L),2)), size='x-small')
    #ax3.grid(which='major', alpha=0.6)
    #ax3.grid(which='minor', alpha=0.4)
    ax3.set_ylabel(R'Displacement Noise $[\mathrm{m} / \sqrt{\mathrm{Hz}}]$')
    ax3.set_xlabel(R'Frequency [Hz]')

    plt.savefig(fpath + 'ETM_TN.pdf')

    

def main(layers = True, trans = True, noise = True):
    # Setup matplotlib params
    plt.rcParams.update({'text.usetex': False,
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
                     'figure.figsize': ((12, 8)),
                     'savefig.dpi': 140,
                     'savefig.bbox': 'tight',
                     'pdf.compression': 9})

    if __debug__:
        print('Loading ' + 'gwinc.ifo' + ' ' + fname)
    plot_layers(layers)
    plot_trans(trans)
    plot_noise(noise)

if __name__ == '__main__':    
    if __debug__:
        print('Make plots...')
    main(layers = True, trans = True, noise = False)
