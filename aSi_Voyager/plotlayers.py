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
from coatingUtils import *

import numpy as np
from timeit import default_timer
#import matplotlib as mpl
import matplotlib.pyplot as plt

import scipy.io as scio
from matplotlib.ticker import FormatStrFormatter
from gwinc import noise

#plt.style.use('bmh')


plt.rcParams.update({'text.usetex': False,
                     'lines.linewidth': 4,
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



# load Data file from run of mkMirror.py
fname = max(glob.iglob('Data/*Layers*.mat'), key=os.path.getctime)
fname = fname[5:] # rm 'data' from the name
#fname = 'ETM_Layers_190519_1459.mat'
if __debug__:
    print('Loading ' + fname + '...')

z       = scio.loadmat('Data/' + fname, squeeze_me=True)
n       = z['n']
T       = z['T']


if __debug__:
    print('Loading ' + 'gwinc.ifo' + ' ' + fname)
    tic = default_timer()

ifo     = gwinc.Struct.from_file(z["ifo_name"])
if __debug__:
    dt = default_timer() - tic
    print('Took ' + str(round(dt,3)) + ' sec to load IFO w/ matlab.')

lambda0 = ifo.Laser.Wavelength # meters

#  plot of the spectral reflectivity
lams    = np.linspace(0.4, 1.6, 300)
rr, _   = multidiel1(n, z['L'], lams)
RR      = np.abs(rr)**2
TT      = 1 - RR

#  convert from optical thickness to physical thickness
L       = lambda0 * op2phys(z['L'], n[1:-1])
#  calculate field v. distance into coating
Z,field = fieldDepth(L, n, pol='p', nPts=300)
layers  = np.cumsum(1e6 * L)
layers  = np.append(0, layers)

if __debug__:
    print('Make plots...')
    tic = default_timer()

fig, ax = plt.subplots(1,1)
ax.semilogy(1e6*lams*lambda0, TT,
                lw=3, label='Transmissivity', c='xkcd:Red')
ax.semilogy(1e6*lams*lambda0, RR,
                lw=3, label='Reflectivity', c='xkcd:electric blue', alpha=0.5)
ax.vlines(lambda0*1e6, T, 1, linestyle='--')
ax.set_xlabel('Wavelength [$\mu \\mathrm{m}$]')
ax.set_ylabel('T or R')
ax.set_ylim((1e-6, 1))
ax.text(lambda0*1.051e6, 1e-1, 'T @ {} um'.format(
    1e6*lambda0), size='x-small')
ax.text(lambda0*1.051e6, 0.5e-1, '= {} ppm'.format(
    round(1e6*T,1)), size='x-small')
ax.legend()
if __debug__:
    print('Transmission of this coating at {} nm is {} ppm'.format(
        1e9*lambda0, round(1e6*T,2)))

plt.savefig('Figures/' + 'ETM_R' + '.pdf')


# Make the plotof the Layer structure
fig , ax = plt.subplots(2,1, sharex=True)
ax[0].plot(Z*1e6,field, color='xkcd:electric purple',
               alpha=0.97, rasterized=False)

#Add some vlines
ax[0].vlines(np.cumsum(L)[1:-1:2]*1e6, 1e-5, 0.55,
            color='xkcd:bright teal', linewidth=0.6,
                 linestyle='--', alpha=0.75, rasterized=False)
ax[0].vlines(np.cumsum(L)[::2]*1e6, 1e-5, 0.55,
            color='xkcd:deep purple', linewidth=0.6,
                 linestyle='--', alpha=0.75,rasterized=False)

#Also visualize the layer thicknesses
ax[1].bar(layers[:-1:2], 1e9*L[::2], width=1e6*L[::2],
        align='edge', color='xkcd:bright teal',
              alpha=0.4, label='$\mathrm{SiO}_2$')
ax[1].bar(layers[1:-1:2],  1e9*L[1::2], width=1e6*L[1::2],
        align='edge', color='xkcd:deep purple',
              alpha=0.4, label='$a-Si$')
ax[1].legend()
ax[1].yaxis.set_major_formatter(FormatStrFormatter("%3d"))
ax[0].set_ylabel('Normalized $|E(z)|^2$')
ax[1].set_ylabel('Physical layer thickness [nm]')
ax[1].set_xlabel('Distance from air interface, $[\mu \mathrm{m}]$')

fig.subplots_adjust(hspace=0.01,left=0.09,right=0.95,top=0.92)
plt.suptitle('a-Si:SiO$_2$ coating electric field')


plt.savefig('Figures/' + fname[:-4] + '.pdf')
plt.savefig('Figures/' + 'ETM_Layers' + '.pdf')


# ----  plot the Thermal Noise
ff = np.logspace(0, 4, 500)
fig3, ax3 = plt.subplots(1,1)
# Build up a "mirror" structure as required by pygwinc
mir = ifo.Optics.ETM
mir.Coating.dOpt = z['L'][:]
StoZ, SteZ, StrZ, _ = gwinc.noise.coatingthermal.coating_thermooptic(ff, 
                                                mir, ifo.Laser.Wavelength, ifo.Optics.ETM.BeamRadius)
SbrZ = gwinc.noise.coatingthermal.coating_brownian(ff, mir, ifo.Laser.Wavelength, ifo.Optics.ETM.BeamRadius)

subBrown = gwinc.noise.substratethermal.substrate_brownian(ff, mir, ifo.Optics.ETM.BeamRadius)
subTE    = gwinc.noise.substratethermal.substrate_thermoelastic(ff, mir, ifo.Optics.ETM.BeamRadius)

Larm = 1 #4000

ax3.loglog(ff, Larm * np.sqrt(StoZ), label='Thermo-Optic', c='xkcd:Purplish Blue')
ax3.loglog(ff, Larm * np.sqrt(SteZ), label='Thermo-Elastic', c='xkcd:Golden')
ax3.loglog(ff, Larm * np.sqrt(StrZ), label='Thermo-Refractive', c='xkcd:Puke')
ax3.loglog(ff, Larm * np.sqrt(SbrZ), label='Brownian', c='xkcd:Tomato')
ax3.loglog(ff, Larm * np.sqrt(subBrown), label='Substrate Brownian', c='xkcd:Dusty Blue')
ax3.loglog(ff, Larm * np.sqrt(subTE), label='Substrate Thermo-Elastic', c='xkcd:Chocolate', alpha=0.3)
ax3.legend()
ax3.set_ylim([8e-24, 2e-20])

#ax3.grid(which='major', alpha=0.6)
#ax3.grid(which='minor', alpha=0.4)
ax3.set_ylabel('Displacement Noise $[\\mathrm{m} / \\sqrt{\\mathrm{Hz}}]$')
ax3.set_xlabel('Frequency [Hz]')
plt.savefig('Figures/' + 'ETM_TN.pdf')

if __debug__:
    dt = default_timer() - tic
    print('Took ' + str(round(dt,3)) + ' sec to make the plots and save them.')
