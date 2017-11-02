#Python script that makes a spectral reflectivity plot of a given coating design

import numpy as np
import matplotlib.pyplot as plt
from mcutils import *
import sys

if 'gvELOG' in plt.style.available:
	plt.style.use('gvELOG')
else:
	plt.style.use('bmh')

dispFileName = '../dispersion_revised.mat'

R_optimal, T_optimal = specREFL('../Data/aLIGO_ETM_HR_20_layers_171030_2100.mat',dispFileName)
R_QW, T_QW = specREFL('../Data/aLIGO_ETM_HR_20_layers_QW.mat',dispFileName)

lambda_real = np.linspace(0.4,1.6,2200) * 1064 #in nm
fig , ax = plt.subplots(nrows=1,ncols=1, figsize=(18,9))
ax.semilogy(lambda_real, R_optimal, linewidth=3, color='navy', alpha=0.7)
ax.semilogy(lambda_real, T_optimal, linewidth=3, color='firebrick',alpha=0.7)
ax.semilogy(lambda_real, R_QW, linewidth=1, color='navy', alpha=0.7, linestyle='--')
ax.semilogy(lambda_real, T_QW, linewidth=1, color='firebrick',alpha=0.7, linestyle='--')
ax.grid('on', which='both',linestyle='--')
ax.grid(which='major',alpha=0.6)
ax.grid(which='minor',alpha=0.4)
ax.set_ylabel('R or T',fontsize=20,fontweight='extra bold')
ax.axvline(x=532,linewidth=2,linestyle='dashed',color='g',label='532 nm')
ax.axvline(x=1064,linewidth=2,linestyle='dashed',color='orangered',label='1064 nm')

plt.savefig('../Figures/aLIGO_ETM_specREFL_comparison.pdf')
