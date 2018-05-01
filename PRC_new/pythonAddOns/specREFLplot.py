#Python script that makes a spectral reflectivity plot of a given coating design

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import sys
sys.path.append('../../generic/pythonAddOns')
from coatingUtils import *

if 'gvELOG' in plt.style.available:
	plt.style.use('gvELOG')
else:
	plt.style.use('bmh')

dispFileName = '../../generic/dispersion_revised.mat'

R_optimal, T_optimal = specREFL('../Data/PR3_HR_20_layers_180430_2138.mat',dispFileName, aoi=41.1, pol='tm')

lambda_real = np.linspace(0.4,1.6,2200) * 1064 #in nm
fig , ax = plt.subplots(nrows=1,ncols=1, figsize=(18,9))
ax.semilogy(lambda_real, R_optimal, linewidth=3, color='navy', alpha=0.7)
ax.semilogy(lambda_real, T_optimal, linewidth=3, color='firebrick',alpha=0.7)
ax.grid('on', which='both',linestyle='--')
ax.grid(which='major',alpha=0.6)
ax.grid(which='minor',alpha=0.4)
ax.set_ylabel('R or T',fontsize=20,fontweight='extra bold')
ax.axvline(x=532,linewidth=2,linestyle='dashed',color='g',label='532 nm')
ax.axvline(x=1064,linewidth=2,linestyle='dashed',color='orangered',label='1064 nm')
ax.set_xlabel('Wavelength [nm]')
ax.set_xticks(np.linspace(400,1600,13))
ax.xaxis.set_major_formatter(FormatStrFormatter("%3d"))
fig.suptitle('40m PR3 spectral reflectivity plot')
fig.subplots_adjust(left=0.07,right=0.95)

fig.savefig('../Figures/40m_PR3_specREFL.pdf')
