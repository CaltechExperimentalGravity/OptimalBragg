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

if 'gvELOG' in plt.style.available:
    plt.style.use('gvELOG')
else:
    plt.style.use('bmh')

#Open the file, load the data
with h5py.File(hdfFileName, 'r') as f:
    samples = np.array(f['MCout'][:])
print(samples.shape)

#Make the plot
fig,ax = plt.subplots(samples.shape[0], samples.shape[0], figsize=(18,18))
corner.corner(samples.T,
        labels=['$\\mathrm{T}_{1064 \\mathrm{nm}}$ [ppm]',
            '$\\mathrm{T}_{532 \\mathrm{nm}}$ [\\%]',
            '$\\mathrm{S}_{\\mathrm{TO}} [\\times 10^{-21} \\mathrm{m}/\\sqrt{\\mathrm{Hz}}]$',
            '$\\mathrm{S}_{\\mathrm{Br}} [\\times 10^{-21} \\mathrm{m}/\\sqrt{\\mathrm{Hz}}]$',
            '$\\vec{E}_{\\mathrm{Surface}}$ [V/m]'],
            quantiles=[0.05, 0.5, 0.95],
            show_titles=True, use_math_text=True,
            bins=50,
            color='firebrick',
            hist_kwargs={'linewidth':2},
            label_kwargs={'fontsize':16, 'fontweight':'bold'},
            title_kwargs={'fontsize':16, 'fontweight':'bold'}, fig=fig)

plt.savefig('Figures/ETM/ETM_nominal_cornerPlt.pdf')
