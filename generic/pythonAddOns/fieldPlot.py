'''
'Python script that makes a spectral reflectivity plot of a given coating design
Example usage:
	fieldPlot.py aLIGO_ETM_20layers.mat
will take the coating design in aLIGO_ETM_20layers.mat and make a plot of the
E-field within the dielectric layer structure.
'''
import numpy as np
import matplotlib.pyplot as plt
from coatingUtils import *
import sys
import scipy.io as scio

if 'gvELOG' in plt.style.available:
	plt.style.use('gvELOG')
else:
	plt.style.use('bmh')

matFileName = sys.argv[1]

data = scio.loadmat(matFileName, struct_as_record=False, squeeze_me=True)
L = 1064e-9*op2phys(data['costOut'].L, data['costOut'].n_IR[1:-1])
n = data['costOut'].n_IR
Z,field = fieldDepth(L,n,pol='p')

#Make the plot
fig , ax = plt.subplots(nrows=1,ncols=1, figsize=(18,9))
ax.plot(Z*1e6,field, color='navy', alpha=0.7)
#Add some vlines to denote layer boundaries
L0 = 0
'''
for ii,kk in enumerate(L):
	ax.vlines(1e6*(L[ii]+L0),1e-5,2, linestyles='--')
	L0 += L[ii]
'''
ax.vlines(np.cumsum(L)[1:-1:2]*1e6, 1e-5, 2, color='xkcd:olive green', linewidth=0.6, linestyle='--')
ax.vlines(np.cumsum(L)[::2]*1e6, 1e-5, 2, color='xkcd:ruby', linewidth=0.6, linestyle='--')
ax.grid('on', which='both',linestyle='--')
ax.grid(which='major',alpha=0.6)
ax.grid(which='minor',alpha=0.4)
ax.set_ylabel('Normalized E field',fontsize=20,fontweight='extra bold')
ax.set_xlabel('Distance from air interface [um]',fontsize=20,fontweight='extra bold')
#ax.set_yscale('log')
#ax.set_ylim([1e-5,2])

plt.suptitle('aLIGO ETM coating electric field distribution')
plt.savefig('../Figures/aLIGO_ETM_fieldStrength.pdf')
