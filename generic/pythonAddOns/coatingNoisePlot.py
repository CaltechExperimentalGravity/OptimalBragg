#Python script that makes plots of coating TN (TE + TR + Brownian) of a given coating design

import numpy as np
import matplotlib.pyplot as plt
from mcutils import *
import sys
import h5py 

if 'gvELOG' in plt.style.available:
	plt.style.use('gvELOG')
else:
	plt.style.use('bmh')

MCFileName = sys.argv[1]
QWFileName = sys.argv[2]

f_optimal = h5py.File(MCFileName,'r')
f_QW = h5py.File(QWFileName,'r')

TOnoise_optimal=np.array(f_optimal['TOnoise'][1:,:])
Brnoise_optimal=np.array(f_optimal['Brnoise'][1:,:])
TOnoise_QW=np.array(f_QW['TOnoise'][1:,:])
Brnoise_QW=np.array(f_QW['Brnoise'][1:,:])


ff = np.logspace(1,3,50)

fig , ax = plt.subplots(nrows=1,ncols=1, figsize=(18,9))

for ii in range(np.shape(TOnoise_optimal)[0]):
    ax.loglog(ff,(TOnoise_optimal[ii][:]),alpha=0.005,color='#ff7f0e',linewidth=0.1,rasterized=True)
    ax.loglog(ff,(Brnoise_optimal[ii][:]),alpha=0.005,color='#d62728',linewidth=0.1,rasterized=True)
    ax.loglog(ff,(TOnoise_QW[ii][:]),alpha=0.005,color='#e377c2',linewidth=0.1,linestyle='--',rasterized=True)
    ax.loglog(ff,(Brnoise_QW[ii][:]),alpha=0.005,color='#8c564b',linewidth=0.1,linestyle='--',rasterized=True)
ax.grid('on', which='both',linestyle='--')
ax.grid(which='major',alpha=0.6)
ax.grid(which='minor',alpha=0.4)
ax.set_ylabel('Displacement Noise $[\\mathrm{m} / \\sqrt{\\mathrm{Hz}}]$',fontsize=20,fontweight='extra bold')

plt.savefig('aLIGO_ETM_TN_comparison.pdf')
