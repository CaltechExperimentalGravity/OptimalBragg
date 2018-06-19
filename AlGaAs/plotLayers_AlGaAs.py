'''
'Python script that plots the layer thicknesses and E field (normalized)
Example usage:
	plotLayers.py aLIGO_ETM_20layers.mat
will take the coating design in aLIGO_ETM_20layers.mat and make a plot of the
E-field within the dielectric layer structure.
'''
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
sys.path.append('../generic/pythonAddOns/')
from coatingUtils import *
import scipy.io as scio
from matplotlib.ticker import FormatStrFormatter

if 'gvELOG' in plt.style.available:
	plt.style.use('gvELOG')
else:
	plt.style.use('bmh')

mpl.rcParams.update({'text.usetex': False,
                     'lines.linewidth': 2.5,
                     'font.size': 14,
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
                     'savefig.dpi': 80,
                     'pdf.compression': 9})

# holy shit, please learn to use ARGPARSE
matFileName = sys.argv[1]


# REMOVE hard-coded wavelength variables !!!
data    = scio.loadmat(matFileName, struct_as_record=False, squeeze_me=True)
L       = 1064e-9*op2phys(data['TNout'].L, data['TNout'].n[1:-1])
n       = data['TNout'].n
Z,field = fieldDepth(L,n,pol='p',nPts=100)
layers  = np.cumsum(1e6*L)
layers  = np.append(0.,layers)

#Make the plot
fig , ax = plt.subplots(2,1,figsize=(12,12),sharex=True)
ax[0].plot(Z*1e6,field, color='xkcd:electric purple', alpha=0.97, rasterized=True)

#Add some vlines
ax[0].vlines(np.cumsum(L)[1:-1:2]*1e6, 1e-5, 0.55, color='xkcd:bright teal',
                 linewidth=0.6, linestyle='--', alpha=0.5, rasterized=True)
ax[0].vlines(np.cumsum(L)[::2]*1e6, 1e-5, 0.55, color='xkcd:deep purple',
                 linewidth=0.6, linestyle='--', alpha=0.5,rasterized=True)

#Also visualize the layer thicknesses
ax[1].bar(layers[:-1:2], 1e9*L[::2], width=1e6*L[::2], align='edge',
              color='xkcd:bright teal', alpha=0.4, label='$\mathrm{GaAs}$')
ax[1].bar(layers[1:-1:2],  1e9*L[1::2], width=1e6*L[1::2], align='edge',
              color='xkcd:deep purple', alpha=0.4, label='$\mathrm{AlGaAs}$')
ax[1].legend()
ax[1].yaxis.set_major_formatter(FormatStrFormatter("%3d"))
ax[0].set_ylabel('Normalized $|E(z)|^2$')
ax[1].set_ylabel('Physical layer thickness [nm]')
ax[1].set_xlabel('Distance from air interface, $z [\mu \mathrm{m}]$')

fig.subplots_adjust(hspace=0.01,left=0.09,right=0.95,top=0.92)
plt.suptitle('AlGaAs coating electric field')


plt.savefig('Figures/' + matFileName[5:] + '.pdf', bbox_inches='tight')
