'''
Python script that makes plots of coating TN (TE + TR + Brownian) of a given coating design
Example usage:
	python coatingNoisePlot.py 'pathToFile1,pathToFile2,...,pathToFileN'
In the above, the list argument is a list of .mat files for which the TN should
be plotted
'''
import numpy as np
import matplotlib.pyplot as plt
from mcutils import *
import sys
import h5py 
import scipy.io as scio
import gwinc

if 'gvELOG' in plt.style.available:
	plt.style.use('gvELOG')
else:
	plt.style.use('bmh')

ff=np.logspace(1,3,100)

fig , ax = plt.subplots(nrows=1,ncols=1, figsize=(18,9))
fileList = sys.argv[1]
for fil in fileList.split(','):
	data = scio.loadmat(fil, struct_as_record=False, squeeze_me=True)
	phi_Ta = data['ifo'].Materials.Coating.Phihighn
	aa,bb,cc,dd = gwinc.noise.coatingthermal.getCoatThermoOptic(ff,data['ifo'],6e-2,np.array([data['costOut'].L]).T)
	SbrZ = gwinc.noise.coatingthermal.getCoatBrownian(ff,data['ifo'],6e-2,np.array([data['costOut'].L]).T)
	ax.loglog(ff,np.sqrt(aa + SbrZ),label='$\\phi_{\\mathrm{Ta}} = %.3E$'%phi_Ta)
ax.legend()
ax.set_ylim([8e-22,2e-20])
ax.grid('on', which='both',linestyle='--')
ax.grid(which='major',alpha=0.6)
ax.grid(which='minor',alpha=0.4)
ax.set_ylabel('Displacement Noise $[\\mathrm{m} / \\sqrt{\\mathrm{Hz}}]$',fontsize=20,fontweight='extra bold')
ax.set_xlabel('Frequency [Hz]',fontsize=24,fontweight='extra bold')
plt.savefig('../Figures/aLIGO_ETM_TN_comparison.pdf')
