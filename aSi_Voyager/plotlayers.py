'''
'Python script that plots the layer thicknesses and E field (normalized)
Example usage:
	plotLayers.py aLIGO_ETM_20layers.mat
will take the coating design in aLIGO_ETM_20layers.mat and make a plot of the
E-field within the dielectric layer structure.
'''
import sys,glob,os
sys.path.append('../../pygwinc/')
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
                     'lines.linewidth': 2.5,
                     'font.family': 'serif',
                     'font.serif': 'Georgia',
                     'font.size': 18,
                     'xtick.labelsize': 'medium',
                     'ytick.labelsize': 'medium',
                     'axes.labelsize': 'medium',
                     'axes.titlesize': 'small',
                     'axes.grid': True,
                     'grid.alpha': 0.73,
                     'lines.markersize': 12,
                     'legend.borderpad': 0.2,
                     'legend.fancybox': True,
                     'legend.fontsize': 13,
                     'legend.framealpha': 0.7,
                     'legend.handletextpad': 0.1,
                     'legend.labelspacing': 0.2,
                     'legend.loc': 'best',
                     'figure.figsize': (10,8),
                     'savefig.dpi': 80,
                     'pdf.compression': 9})


# Constants for calculating absorption
alpha_SiO2 = 1e-3 # https://journals.aps.org/prd/pdf/10.1103/PhysRevD.91.042002
alpha_aSi = 1e-2 #Guess?


# load Data file from run of mkMirror.py
fname = max(glob.iglob('Data/*Layers*.mat'), key=os.path.getctime)
fname = fname[5:] # rm 'data' from the name
#fname = 'ETM_Layers_190517_1828.mat'
if __debug__:
    print('Loading ' + fname + '...')
   

z       = scio.loadmat('Data/' + fname, squeeze_me=True)
n       = z['n']
T       = z['T']


if __debug__:
    print('Loading ' + 'gwinc.ifo' + '...')
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

#convert from optical thickness to physical thickness
L       = lambda0 * op2phys(z['L'], n[1:-1])
# calculate field v. distance into coating
Z,field = fieldDepth(L, n, pol='p', nPts=300, lam=ifo.Laser.Wavelength)
layers  = np.cumsum(1e6 * L)
layers  = np.append(0, layers)

if __debug__:
    print('Make plots...')
fig, ax = plt.subplots(1,1)
ax.semilogy(1e6*lams*lambda0, TT,
                lw=3, label='Transmissivity', c='xkcd:Red')
ax.semilogy(1e6*lams*lambda0, RR,
                lw=3, label='Reflectivity', c='xkcd:electric blue', alpha=0.5)
ax.vlines(lambda0*1e6, T, 1, linestyle='--')
ax.set_xlabel('Wavelength [$\mu \\mathrm{m}$]')
ax.set_ylabel('T or R')
ax.set_ylim((1e-6, 1))
ax.text(lambda0*1.1e6, 1e-1, 'T @ {} um'.format(
    1e6*lambda0), size='x-small')
ax.text(lambda0*1.1e6, 0.7e-1, '= {} ppm'.format(
    round(1e6*T,1)), size='x-small')
ax.legend()
if __debug__:
    print('Transmission of this coating at {} nm is {} ppm'.format(
        1e9*lambda0, round(1e6*T,2)))

fig.savefig('Figures/' + fname[:-4] + '_R' + '.pdf', bbox_inches='tight')


# Make the plotof the Layer structure
fig2 , ax2 = plt.subplots(2,1, sharex=True)
ax2[0].plot(Z*1e6,field, color='xkcd:electric purple',
               alpha=0.97, rasterized=False)
absStr = f'$|\\vec E_{{\mathrm{{surface}}}}| = {1e6*field[0]:.0f}$ ppm of $\\vec |E_{{\mathrm{{inc}}}}|$'
absStr += '\n'
intAbs = calcAbsorption(field, L, 300, 1e-3, 1e-2)
absStr += f'Integrated absorption in stack is {intAbs:.3f} ppm'
print(f'Total integrated absorption for this stack is {intAbs:.3f} ppm, assuming absorption in SiO2 is {alpha_SiO2:.1E}/m and that in a-Si is {alpha_aSi:.1E}/m.')
ax2[0].text(0.5,0.7,absStr,transform=ax2[0].transAxes, fontsize=14)

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
              alpha=0.4, label='$a-Si$')
ax2[1].legend()
ax2[1].yaxis.set_major_formatter(FormatStrFormatter("%3d"))
ax2[0].set_ylabel('Normalized $|E(z)|^2$')
ax2[1].set_ylabel('Physical layer thickness [nm]')
ax2[1].set_xlabel('Distance from air interface, $[\mu \mathrm{m}]$')

fig2.subplots_adjust(hspace=0.01,left=0.09,right=0.95,top=0.92)
fig2.suptitle('a-Si:SiO$_2$ coating electric field')


fig2.savefig('Figures/' + fname[:-4] + '.pdf', bbox_inches='tight')


# plot the Thermal Noise
ff = np.logspace(0, 4, 500)
fig3, ax3 = plt.subplots(1,1)
# Build up a "mirror" structure as required by pygwinc
mir = ifo.Optics.ETM
mir.Coating.dOpt = z['L'][:]
StoZ, SteZ, StrZ, _ = gwinc.noise.coatingthermal.coating_thermooptic(ff, 
                                                mir, ifo.Laser.Wavelength, ifo.Optics.ETM.BeamRadius)
SbrZ = gwinc.noise.coatingthermal.coating_brownian(ff, mir, ifo.Laser.Wavelength, ifo.Optics.ETM.BeamRadius)
ax3.loglog(ff, np.sqrt(StoZ), label='Thermo-Optic')
ax3.loglog(ff, np.sqrt(SteZ), label='Thermo-Elastic')
ax3.loglog(ff, np.sqrt(StrZ), label='Thermo-Refractive')
ax3.loglog(ff, np.sqrt(SbrZ), label='Brownian')
ax3.legend()
#ax3.set_ylim([8e-22, 2e-20])

ax3.grid(which='major', alpha=0.6)
ax3.grid(which='minor', alpha=0.4)
ax3.set_ylabel('Displacement Noise $[\\mathrm{m} / \\sqrt{\\mathrm{Hz}}]$')
ax3.set_xlabel('Frequency [Hz]')
plt.savefig('Figures/' + fname[:-4]+ '_TN.pdf')
