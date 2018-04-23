#!/usr/bin/env python

from __future__ import division

import matplotlib as mpl
#mpl.use('Agg')

import matplotlib.pyplot as plt
#plt.style.use('seaborn-paper')
#plt.style.use('fivethirtyeight')
import numpy as np
from scipy.io import loadmat
from mpl_toolkits.axes_grid1 import make_axes_locatable


mpl.rcParams.update({'text.usetex': False,
                     'lines.linewidth': 2.5,
                     'font.size': 14,
                     'xtick.labelsize': 'x-small',
                     'ytick.labelsize': 'x-small',
                     'axes.labelsize': 'x-small',
                     'axes.titlesize': 'small',
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

def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)

# load the data from tune_doAlGaAs.m
dat = loadmat('Data/ETM_tune_180422_1055.mat')

outs = dat['scores']

SelfAdj = dat['selfies'][0]
SocAdj  = dat['socialies'][0]


# sum over the 2nd axis (repeated trials)
a,b,N,d = outs.shape
zz = np.sum(outs, axis=2)
t_norm = 60  # sec to min

brown = zz[:,:,0]/N/t_norm
TO = zz[:,:,1]/N/t_norm
tot = np.sum(zz[:,:,0:3], axis=2)/N/t_norm
tocs  = zz[:,:,4]/N/t_norm

fig,ax = plt.subplots(2, 2, sharex=True, sharey=True)

cmap = plt.cm.get_cmap('plasma', 15)    # 11 discrete colors

z = np.zeros((len(SelfAdj) * len(SocAdj), 6))

zi = 0
for jj in range(N):
    for kk in range(N):
        
        z[zi] = [SelfAdj[jj], SocAdj[kk],
                     brown[jj,kk], TO[jj,kk], tot[jj,kk], tocs[jj,kk]]
        zi += 1

        
ax00 = ax[0,0].scatter(z[:,0], z[:,1], c = z[:,2], alpha=0.97, cmap=cmap)

ax01 = ax[0,1].scatter(z[:,0], z[:,1], c = z[:,3], alpha=0.97, cmap=cmap)
ax10 = ax[1,0].scatter(z[:,0], z[:,1], c = z[:,4], alpha=0.97, cmap=cmap)
ax11 = ax[1,1].scatter(z[:,0], z[:,1], c = z[:,5], alpha=0.97, cmap=cmap)

#plt.clim(tocs.min(), tocs.max())
#clb = plt.colorbar()
#clb.set_label('Time to Complete [min]')

ax[1,0].set_xlabel('Self Adj')
ax[1,1].set_xlabel('Self Adj')
ax[0,0].set_ylabel('Soc Adj')
ax[1,0].set_ylabel('Soc Adj')

ax[0,0].set_title('Brownian Noise')
ax[0,1].set_title('Thermo-Optic Noise')
ax[1,0].set_title('Total Score')
ax[1,1].set_title('Time to Complete')

colorbar(ax00)
colorbar(ax01)
colorbar(ax10)
colorbar(ax11)

plt.tight_layout(h_pad=1)
#plt.grid(True)
#plt.xlim([0,90])
#plt.ylim([0,1])
#ax[1].yaxis.set_ticks(np.arange(-180, 181, 45))
#plt.legend()
#plt.show()
plt.savefig("Figures/hyperTune.pdf", bbox_inches='tight')


