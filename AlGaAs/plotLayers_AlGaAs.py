'''
'Python script that plots the layer thicknesses and E field (normalized)
Example usage:
	plotLayers.py aLIGO_ETM_20layers.mat
will take the coating design in aLIGO_ETM_20layers.mat and make a plot of the
E-field within the dielectric layer structure.
'''
import numpy as np
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

matFileName = sys.argv[1]

data = scio.loadmat(matFileName, struct_as_record=False, squeeze_me=True)
L = 1064e-9*op2phys(data['TNout'].L, data['TNout'].n[1:-1])
n = data['TNout'].n
Z,field = fieldDepth(L,n,pol='p',nPts=100)
layers = np.cumsum(1e6*L)
layers = np.append(0.,layers)

#Make the plot
fig , ax = plt.subplots(2,1,figsize=(12,12),sharex=True)
ax[0].plot(Z*1e6,field, color='xkcd:bright red', alpha=0.7, rasterized=True)
#Add some vlines
ax[0].vlines(np.cumsum(L)[1:-1:2]*1e6, 1e-5, 0.55, color='xkcd:bright teal', linewidth=0.6, linestyle='--', rasterized=True)
ax[0].vlines(np.cumsum(L)[::2]*1e6, 1e-5, 0.55, color='xkcd:deep purple', linewidth=0.6, linestyle='--', rasterized=True)
#Also visualize the layer thicknesses
ax[1].bar(layers[:-1:2], 1e9*L[::2],width=1e6*L[::2],align='edge',color='xkcd:bright teal', alpha=0.4, label='$\mathrm{GaAs}$')
ax[1].bar(layers[1:-1:2], 1e9*L[1::2],width=1e6*L[1::2],align='edge',color='xkcd:deep purple', alpha=0.4, label='$\mathrm{AlGaAs}$')
ax[1].legend(loc='best')
ax[1].yaxis.set_major_formatter(FormatStrFormatter("%3d"))
ax[0].set_ylabel('Normalized E field',fontsize=20,fontweight='extra bold')
ax[1].set_ylabel('Physical layer thickness [nm]',fontsize=20,fontweight='extra bold')
ax[1].set_xlabel('Distance from air interface [um]',fontsize=20,fontweight='extra bold')
fig.subplots_adjust(hspace=0.01,left=0.09,right=0.95,top=0.92)
plt.suptitle('AlGaAs coating electric field')
plt.savefig('Figures/AlGaAs_fieldStrength.pdf')
