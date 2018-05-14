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
if len(sys.argv)>2:
	print(len(sys.argv))
	hdf5QW = sys.argv[2]
	f_QW = h5py.File(hdfFileName,'r')
	samples_QW = np.array(f_QW['MCout'][:])

if 'gvELOG' in plt.style.available:
    plt.style.use('gvELOG')
else:
    plt.style.use('bmh')

#Open the file, load the data
f = h5py.File(hdfFileName,'r')
samples=np.array(f['MCout'][:])

#Make the plot
fig,ax = plt.subplots(np.shape(samples)[0], np.shape(samples)[0], figsize=(18,18))
#fig.subplots_adjust(wspace=0.5,hspace=0.35)
corner.corner(samples.T,
        labels=['$\\mathrm{T}_{1064 }$ [ppm]', 
            '$\\mathrm{R}^{p-pol}_{532 }$ [\\%]', 
            '$\\mathrm{R}^{s-pol}_{532}$ [\\%]', 
            '$\\frac{\\vec{E}_{\\mathrm{Surface}}}{1 \\mathrm{V/m}}$'],
            #quantiles=[0.9, 0.95, 0.98],
            show_titles=True, use_math_text=True,
            bins=50,
            range=[(15.,60.),(0.,15),(5.,20.),(0.,0.6)],
            #   levels=(0.95,),
	    color='xkcd:deep blue',
	    hist_kwargs={'linewidth':3},
            label_kwargs={'fontsize':26, 'fontweight':'bold'},
            title_kwargs={'fontsize':26, 'fontweight':'bold'}, fig=fig)
if len(sys.argv)>2:
	corner.corner(samples_QW.T,
	range=[(15.,60.),(0.,15),(5.,20.),(0.,0.6)],
	color='xkcd:olive green',
	hist_kwargs={'linewidth':3,'alpha':0.7},
	label_kwargs={'fontsize':26, 'fontweight':'bold'},
	title_kwargs={'fontsize':26, 'fontweight':'bold'}, fig=fig)
#Make the font size readable
for aa in ax:
	for bb in aa:
		bb.tick_params(labelsize='xx-large')
plt.savefig('../Figures/PR3_cornerPlt.pdf')
