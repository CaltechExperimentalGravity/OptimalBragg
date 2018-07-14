import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib as mpl
import corner
from scipy.io import loadmat
import os
import glob
import argparse

mpl.rcParams.update({'text.usetex': False,
                     'lines.linewidth': 2.5,
                     'font.size': 10,
                     'xtick.labelsize': 'small',
                     'ytick.labelsize': 'small',
                     'axes.labelsize': 'x-small',
                     'axes.titlesize': 'x-small',
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


newest = max(glob.iglob('MCout/*.[Hh]5'), key=os.path.getctime)

parser = argparse.ArgumentParser(description='Script to take output of MC analysis and make a corner plot')
parser.add_argument("-f", "--filename", type=str, default=newest,
                    help="HDF5 File with MCMC data (defaults to newest by datestring in MCout directory)")
args = parser.parse_args()

hdfFileName = args.filename

# Open the file, load the data
f = h5py.File(hdfFileName,'r')
samples = np.array(f['MCout'][:])
samples[2,-1] = samples[2,-2]
samples[3,-1] = samples[3,-2]

# Eventually, these parameters will be loaded from the hdf5 file itself
N = 1e5
sigma_nLow = 0.005
sigma_nHigh = 0.005
parText = "$N_s = ${} \n"\
	"$\sigma_nL$ = {} \n"\
	"$\sigma_nH$ = {}".format(N, sigma_nLow, sigma_nHigh)

# TODO: add text about deposition errors on the plot in the unused space
# Make the plot
fig,ax = plt.subplots(np.shape(samples)[1],
                      np.shape(samples)[1], figsize=(12,12))
corner.corner(samples,
        labels=['$\\mathrm{T}_{1064} \\mathrm{[ppm]}$', 
            '$\\mathrm{S}_{\\mathrm{TO}} [\\times 10^{-21} \\mathrm{m}/\\sqrt{\\mathrm{Hz}}]$',
            '$\\mathrm{S}_{\\mathrm{Br}} [\\times 10^{-21} \\mathrm{m}/\\sqrt{\\mathrm{Hz}}]$',
            '$\\vec{E}_{\\mathrm{Surf}} \\mathrm{[V/m]}$',
            '$\\mathrm{Absorption [ppm]}$'],
            #quantiles=[0.9, 0.95, 0.98],
            truths = [5, 1, 2, 1, 0.5],
            truth_color = 'xkcd:duck egg blue',
            show_titles = True,
            use_math_text = True,
            bins = 50,
            #range=[(0,20), (0, 8), (2.0, 2.5), (0, 20), (0, 3)],
            #   levels=(0.95,),
            color = 'xkcd:cerulean blue',
            smooth = 0.5, # smoothing scale in std's
            hist_kwargs  = {'linewidth':2.5},
            label_kwargs = {'fontsize':'large', 'fontweight':'bold'},
            title_kwargs = {'fontsize':'medium', 'fontweight':'bold'},
                  fig = fig)

# Print the MC parameters onto the plot
ax[0,3].text(0.005,0.27,parText,wrap=True,transform=ax[0,3].transAxes)
pdfFile = 'Figures/'+ hdfFileName.strip('.h5').strip('MCout/') + '.pdf'
print("File saved as " + pdfFile)
fig.savefig(pdfFile)
