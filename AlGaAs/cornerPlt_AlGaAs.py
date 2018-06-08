'''
Script to take output of MC analysis and make a corner plot
Example usage:
    python cornerPlt.py aLIGO_ETM_MC.hdf5

The above will make a corner plot from the data in aLIGO_ETM_MC.hdf5

'''

import numpy as np
import h5py
import matplotlib.pyplot as plt
#import sys
import matplotlib
import corner
from scipy.io import loadmat
import os
import glob

import argparse

newest = max(glob.iglob('MCout/*.[Hh]5'), key=os.path.getctime)

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filename", type=str, default=newest,
                    help="file with MCMC data")
args = parser.parse_args()

# use argparse instead
hdfFileName = args.filename

#if 'gvELOG' in plt.style.available:
#    plt.style.use('gvELOG')
#else:
#    plt.style.use('bmh')



#Open the file, load the data
f = h5py.File(hdfFileName,'r')
#f = loadmat(hdfFileName)
samples = np.array(f['MCout'][:])
samples[2,-1] = samples[2,-2]
samples[3,-1] = samples[3,-2]

# TODO: add text about deposition errors on the plot in the unused space
# Make the plot
fig,ax = plt.subplots(np.shape(samples)[1], np.shape(samples)[1], figsize=(12,12))
#fig.subplots_adjust(wspace=0.5,hspace=0.35)
corner.corner(samples,
        labels=['$\\mathrm{T}_{1064} \\mathrm{[ppm]}$', 
            '$\\mathrm{S}_{\\mathrm{TO}} [\\times 10^{-21} \\mathrm{m}/\\sqrt{\\mathrm{Hz}}]$',
            '$\\mathrm{S}_{\\mathrm{Br}} [\\times 10^{-21} \\mathrm{m}/\\sqrt{\\mathrm{Hz}}]$',
            '$\\vec{E}_{\\mathrm{Surf}} \\mathrm{[V/m]}$'],
            #quantiles=[0.9, 0.95, 0.98],
            truths = [5, 1, 2.3, 2],
            show_titles=True, use_math_text=True,
            bins = 50,
            range=[(0,20), (0,1.4), (2.2,2.4), (0,15)],
            #   levels=(0.95,),
            color = 'xkcd:browny orange',
            smooth = 1,
            hist_kwargs  = {'linewidth':2.5},
            label_kwargs = {'fontsize':'large', 'fontweight':'bold'},
            title_kwargs = {'fontsize':'large', 'fontweight':'bold'},
                  fig = fig)

fubu = hdfFileName + '.pdf'
print("File saved as " + fubu)
plt.savefig(fubu)
