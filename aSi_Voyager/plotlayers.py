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
#import matplotlib as mpl
import matplotlib.pyplot as plt

import scipy.io as scio
from matplotlib.ticker import FormatStrFormatter

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





# REMOVE hard-coded wavelength variables !!!
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
ifo     = gwinc.load_ifo(z["ifo_name"])
lambda0 = ifo.Laser.Wavelength # meters

#  plot of the spectral reflectivity
lams    = np.linspace(0.4, 1.6, 300)
rr, _   = multidiel1(n, z['L'], lams)
RR      = np.abs(rr)**2
TT      = 1 - RR

#convert from optical thickness to physical thickness
L       = lambda0 * op2phys(z['L'], n[1:-1])
# calculate field v. distance into coating
Z,field = fieldDepth(L, n, pol='p', nPts=100)
layers  = np.cumsum(1e6 * L)
layers  = np.append(0, layers)

if __debug__:
    print('Make plots...')
fig, ax = plt.subplots(1,1)
ax.semilogy(1e6*lams*lambda0, TT,
                lw=3, label='Transmissivity', c='xkcd:Teal')
ax.semilogy(1e6*lams*lambda0, RR,
                lw=3, label='Reflectivity', c='xkcd:Red', alpha=0.5)
ax.vlines(lambda0*1e6, T, 1, linestyle='--')
ax.set_xlabel('Wavelength [$\mu \\mathrm{m}$]')
ax.set_ylabel('T or R')
ax.set_ylim((T,1))
ax.text(lambda0*1.1e6, 1e-1, 'T @ {} um'.format(
    1e6*lambda0), size='x-small')
ax.text(lambda0*1.1e6, 0.7e-1, '= {} ppm'.format(
    round(1e6*T,3)), size='x-small')
ax.legend()
if __debug__:
    print('Transmission of this coating at {} nm is {} ppm'.format(
        1e9*lambda0, round(1e6*T,3)))

plt.savefig('Figures/' + 'ETM_R' + '.pdf', bbox_inches='tight')


# Make the plotof the Layer structure
fig , ax = plt.subplots(2,1, sharex=True)
ax[0].plot(Z*1e6,field, color='xkcd:electric purple',
               alpha=0.97, rasterized=False)

#Add some vlines
ax[0].vlines(np.cumsum(L)[1:-1:2]*1e6, 1e-5, 0.55,
            color='xkcd:bright teal', linewidth=0.6,
                 linestyle='--', alpha=0.75, rasterized=False)
ax[0].vlines(np.cumsum(L)[::2]*1e6, 1e-5, 0.55,
            color='xkcd:deep purple',
            linewidth=0.6, linestyle='--', alpha=0.75,rasterized=False)

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
ax[1].set_xlabel('Distance from air interface, $z [\mu \mathrm{m}]$')

fig.subplots_adjust(hspace=0.01,left=0.09,right=0.95,top=0.92)
plt.suptitle('a-Si:SiO$_2$ coating electric field')


plt.savefig('Figures/' + fname[:-4] + '.pdf', bbox_inches='tight')
plt.savefig('Figures/' + 'ETM_Layers' + '.pdf', bbox_inches='tight')
