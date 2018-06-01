'''
Script to take output of MC analysis and make a corner plot
Example usage:
    python cornerPlt.py aLIGO_ETM_MC.hdf5

The above will make a corner plot from the data in aLIGO_ETM_MC.hdf5

'''

import numpy as np
import h5py
import matplotlib.pyplot as plt
import sys
import matplotlib
import corner

hdfFileName = sys.argv[1]

#if 'gvELOG' in plt.style.available:
#    plt.style.use('gvELOG')
#else:
#    plt.style.use('bmh')

#Open the file, load the data
f = h5py.File(hdfFileName,'r')
samples=np.array(f['MCout'][:])
samples[2,-1] = samples[2,-2]
samples[3,-1] = samples[3,-2]

#Make the plot
fig,ax = plt.subplots(np.shape(samples)[1], np.shape(samples)[1], figsize=(18,18))
#fig.subplots_adjust(wspace=0.5,hspace=0.35)
corner.corner(samples,
        labels=['$\\mathrm{T}_{1064 \\mathrm{nm}}$ [ppm]', 
            '$\\mathrm{S}_{\\mathrm{TO}} [\\times 10^{-21} \\mathrm{m}/\\sqrt{\\mathrm{Hz}}]$',
            '$\\mathrm{S}_{\\mathrm{Br}} [\\times 10^{-21} \\mathrm{m}/\\sqrt{\\mathrm{Hz}}]$',
            '$\\vec{E}_{\\mathrm{Surface}}$ [V/m]'],
            #quantiles=[0.9, 0.95, 0.98],
            show_titles=True, use_math_text=True,
            bins=50,
            range=[(0.,15.),(0.4,1.4),(2.27,2.34),(5,12)],
            #   levels=(0.95,),
	    color='xkcd:bright blue',
            hist_kwargs={'linewidth':2},
            label_kwargs={'fontsize':16, 'fontweight':'bold'},
            title_kwargs={'fontsize':16, 'fontweight':'bold'}, fig=fig)

plt.savefig('Figures/AlGaAs_cornerPlt.pdf')
