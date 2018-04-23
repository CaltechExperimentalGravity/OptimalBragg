#!/usr/bin/env python

from __future__ import division

import matplotlib as mpl
#mpl.use('Agg')

import matplotlib.pyplot as plt
#plt.style.use('seaborn-paper')
#plt.style.use('fivethirtyeight')
import numpy as np
from scipy.io import loadmat

mpl.rcParams.update({'text.usetex': False,
                     'lines.linewidth': 2.5,
                     'font.size': 14,
                     'xtick.labelsize': 'x-small',
                     'ytick.labelsize': 'x-small',
                     'axes.labelsize': 'x-small',
                     'axes.titlesize': 'large',
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
                     'savefig.dpi': 240,
                     'pdf.compression': 9})



# load the data from tune_doAlGaAs.m
dat = loadmat('Data/ETM_tune_180422_1055.mat')

outs = dat['scores']

SelfAdj = dat['selfies'][0]
SocAdj  = dat['socialies'][0]


# sum over the 2nd axis (repeated trials)
a,b,N,d = outs.shape
zz = np.sum(outs, axis=2)
tocs = zz[:,:,4]/N/60

fig = plt.figure(772)
cmap = plt.cm.get_cmap('plasma', 15)    # 11 discrete colors

z = np.zeros((len(SelfAdj) * len(SocAdj), 3))

zi = 0
for jj in range(N):
    for kk in range(N):
        val = tocs[jj,kk]
        z[zi] = [SelfAdj[jj], SocAdj[kk], val]
        zi += 1

        
plt.scatter(z[:,0], z[:,1], c = z[:,2], alpha=0.97, cmap=cmap)

plt.clim(tocs.min(), tocs.max())
clb = plt.colorbar()
clb.set_label('Time to Complete [min]')

plt.xlabel('Self Adj')
plt.ylabel(r'Social Adj')
plt.title('PSO Tuning for AlGaAs')
#plt.grid(True)
#plt.xlim([0,90])
#plt.ylim([0,1])
#ax[1].yaxis.set_ticks(np.arange(-180, 181, 45))
plt.legend()
#plt.show()
plt.savefig("Figures/hyperTune.pdf", bbox_inches='tight')


