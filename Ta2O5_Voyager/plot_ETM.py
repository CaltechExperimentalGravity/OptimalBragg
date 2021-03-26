'''
'Python script that plots the layer thicknesses and E field (normalized)
Example usage:
	plotLayers.py aLIGO_ETM_20layers.mat
will take the coating design in aLIGO_ETM_20layers.mat and make a plot of the
E-field within the dielectric layer structure.
'''
import sys,glob,os

# installed through conda rather than direct import
#sys.path.append('../../pygwinc/')

sys.path.append('../generic/')
from generic_local.coatingUtils import *

import numpy as np
from timeit import default_timer
#import matplotlib as mpl
import matplotlib.pyplot as plt

import scipy.io as scio
from matplotlib.ticker import FormatStrFormatter
from gwinc import noise, Struct

#plt.style.use('bmh')
paramfilename = 'ETM_params.yml'
opt_params = importParams(paramfilename)
Npairs = opt_params['Npairs']
Nfixed = opt_params['Nfixed']
Nlayers = 2*Npairs + 1
N_particles = opt_params['Nparticles']

plt.rcParams.update({'text.usetex': False,
                     'lines.linewidth': 3,
                     'font.family': 'serif',
                     'font.serif': 'Georgia',
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


# Constants for calculating absorption
alpha_SiO2 = 1e-3 # (get a better # for 2 um and 123 K) https://journals.aps.org/prd/pdf/10.1103/PhysRevD.91.042002
alpha_aSi = 100e-6 / 1e-6 # Figs 2/3, https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.120.263602#page=3

# load Data file from run of mkITM.py
if len(sys.argv) == 1:
    fname = max(glob.iglob('Data/ETM/*Layers*.mat'), key=os.path.getctime)
    fname = fname[5:] # rm 'data' from the name
    z = scio.loadmat('Data/' + fname, squeeze_me=True)
else:
    fname = str(sys.argv[1])
    # fname = fname[fname.find('ETM/')::]
    z = scio.loadmat(fname, squeeze_me=True)

#fname = 'ETM_Layers_190519_1459.mat'
if __debug__:
    print('Loading ' + fname + '...')

costOut = z['vectorCost']
n       = z['n']
T       = z['T']
Taux    = z['Taux']


if __debug__:
    print('Loading ' + 'gwinc.ifo' + ' ' + fname)
    tic = default_timer()

ifo     = Struct.from_file(z["ifo_name"])
if __debug__:
    dt = default_timer() - tic
    print('Took ' + str(round(dt,3)) + ' sec to load IFO w/ matlab.')

lambda0 = ifo.Laser.Wavelength # meters
lambda1 = ifo.Laser.Wavelength * 2/3
lamb1550 = 1550/2128.2
lamb632 = 632/2128.2

#  plot of the spectral reflectivity
rr1550, _ = multidiel1(n, z['L'], lamb1550)
rr632, _ = multidiel1(n, z['L'], lamb632) 
T1550 = 1 - np.abs(rr1550[0])**2
T632 = 1 - np.abs(rr632[0])**2
lams    = np.linspace(0.2, 1.8, 512)
rr, _   = multidiel1(n, z['L'], lams)
RR      = np.abs(rr)**2
TT      = 1 - RR

# Build stats based on various figures for star fish chart
stats = {}
for c, s, w in zip(opt_params['costs'], 
                 costOut, 
                 opt_params['weights']):
    stat = 1 / s
    if s > 1e2 or s < 1e-5:
        stat = 1e-2
    stats[c] = np.abs(np.log(np.abs(stat)))

print(stats)
from starfish import polar_cost
polar_cost(stats, scale=10,
  fname=r'Figures/ETM/ETM_SF'+ fname[-16:-4]+'.png',
  figtitle=fR'Stack with {Nlayers} layers ({Nfixed} fixed bilayers) and {N_particles} particles')


#  convert from optical thickness to physical thickness
L       = lambda0 * op2phys(z['L'], n[1:-1])

# calculate field v. distance into coating
Z,field = fieldDepth(L, n, pol='p', nPts=300, lam=ifo.Laser.Wavelength)

layers  = np.cumsum(1e6 * L)
layers  = np.append(0, layers)

if __debug__:
    print('Make plots...')
    tic = default_timer()

fig, ax = plt.subplots(1,1)
ax.semilogy(1e6*lams*lambda0, TT,
                lw=1.5, label='Transmissivity', c='xkcd:Red')
ax.semilogy(1e6*lams*lambda0, RR,
                lw=1.5, label='Reflectivity', c='xkcd:electric blue', alpha=0.5)
ax.vlines(lambda0*1e6, T, 1, linestyle='--')
ax.vlines(lambda1*1e6, Taux, 1, linestyle='--')
ax.vlines(lamb632*lambda0*1e6, T632, 1, linestyle='-.')
ax.vlines(lamb1550*lambda0*1e6, T1550, 1, linestyle='-.')
ax.set_xlabel('Wavelength [$\mu \\mathrm{m}$]')
ax.set_ylabel('T or R')
ax.set_ylim((1e-6, 1))
ax.text(lambda0*1.051e6, 1e-2, f'T @ {1e6*lambda0:.4f} um', size='x-small')
ax.text(lambda0*1.051e6, 0.5e-2, f'= {round(1e6*T,1):.2f} ppm', size='x-small')
ax.text(lambda1*0.851e6, 1e-2, f'T @ {lambda1*1e6:.4f} um', size='x-small')
ax.text(lambda1*0.851e6, 0.5e-2, f'= {round(1e6*Taux,1):.2f} ppm', size='x-small')
# Oplevs stuff
ax.text(lamb1550*lambda0*1.051e6, 1e-3, f'T @ {1e6*lambda0*lamb1550:.3f} um', size='x-small')
ax.text(lamb1550*lambda0*1.051e6, 0.5e-3, f'= {1e6*T1550:.2f} ppm', size='x-small')
ax.text(lamb632*lambda0*1e6, 1e-3, f'T @ {lamb632*lambda0*1e6:.3f} um', size='x-small')
ax.text(lamb632*lambda0*1e6, 0.5e-3, f'= {1e6*T632:.2f} ppm', size='x-small')

im = plt.imread('./Figures/ETM/ETM_SF'+ fname[-16:-4]+'.png')

newax = fig.add_axes([0.55, 0.1, 0.4, 0.4], anchor='NE')
newax.imshow(im)
newax.axis('off')

ax.legend()
if __debug__:
    print('Transmission of this coating at {} nm is {} ppm'.format(
        1e9*lambda0, round(1e6*T,2)))


plt.savefig('./Figures/ETM/' + 'ETM_R' + fname[-16:-4] + '.pdf')
plt.savefig('./Figures/ETM/' + 'ETM_R' + '.pdf')

# Make the plotof the Layer structure
fig2 , ax2 = plt.subplots(2,1, sharex=True)
ax2[0].plot(Z*1e6,field, color='xkcd:electric purple',
               alpha=0.97, rasterized=False)
absStr = f'$|\\vec E_{{\mathrm{{surface}}}}| = {1e6*field[0]:.0f}$ ppm of $\\vec |E_{{\mathrm{{inc}}}}|$'
absStr += '\n'
intAbs = calcAbsorption(field, L, 300, alpha_SiO2, alpha_aSi)
absStr += f'Integrated absorption in stack is {intAbs:.1f} ppm'
# print(f'Total integrated absorption {intAbs:.3f} ppm, assuming abs in SiO2 is {alpha_SiO2:.1E}/m and a-Si is {alpha_aSi:.1E}/m.')
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
              alpha=0.4, label='$\mathrm{SiO}_2$')
ax2[1].bar(layers[1:-1:2],  1e9*L[1::2], width=1e6*L[1::2],
        align='edge', color='xkcd:deep purple',
              alpha=0.4, label='$Ta2O5$')
ax2[1].legend()
ax2[1].yaxis.set_major_formatter(FormatStrFormatter("%3d"))
ax2[0].set_ylabel('Normalized $|E(z)|^2$')
ax2[1].set_ylabel('Physical layer thickness [nm]')
ax2[1].set_xlabel('Distance from air interface, $[\mu \mathrm{m}]$')


fig2.subplots_adjust(hspace=0.01,left=0.09,right=0.95,top=0.92)
fig2.suptitle('Ta2O5:SiO$_2$ coating electric field')

plt.savefig('./Figures/ETM/ETM_Layers' + fname[-16:-4] + '.pdf')
plt.savefig('./Figures/ETM/' + 'ETM_Layers' + '.pdf')


# ----  plot the Thermal Noise
ff = np.logspace(0, 4, 500)
fig3, ax3 = plt.subplots(1,1)
# Build up a "mirror" structure as required by pygwinc
mir = ifo.Optics.ETM
mir.Coating.dOpt = z['L'][:]
StoZ, SteZ, StrZ, _ = noise.coatingthermal.coating_thermooptic(ff, 
                                                mir, ifo.Laser.Wavelength, ifo.Optics.ETM.BeamRadius)
SbrZ = noise.coatingthermal.coating_brownian(ff, mir, ifo.Laser.Wavelength, ifo.Optics.ETM.BeamRadius)

subBrown = noise.substratethermal.substrate_brownian(ff, mir, ifo.Optics.ETM.BeamRadius)
subTE    = noise.substratethermal.substrate_thermoelastic(ff, mir, ifo.Optics.ETM.BeamRadius)

Larm = 1 #4000

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
ax3.set_ylabel('Displacement Noise $[\\mathrm{m} / \\sqrt{\\mathrm{Hz}}]$')
ax3.set_xlabel('Frequency [Hz]')

plt.savefig('Figures/ETM/' + 'ETM_TN.pdf')

if __debug__:
    dt = default_timer() - tic
    print('Took ' + str(round(dt, 1)) + ' sec to make the plots and save them.')

