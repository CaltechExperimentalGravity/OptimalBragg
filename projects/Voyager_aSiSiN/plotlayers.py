"""Plot layer structure, E-field, spectral reflectivity, and thermal noise.

Loads the most recent optimization result from Data/ and generates plots.

Usage:
    cd projects/Voyager_aSiSiN
    python plotlayers.py
"""
import sys
import glob
import os

import numpy as np
from timeit import default_timer
import matplotlib.pyplot as plt
import scipy.io as scio
from matplotlib.ticker import FormatStrFormatter

from OptimalBragg.layers import multidiel1, op2phys, fieldDepth, calc_abs
from OptimalBragg.noise import (
    coating_thermooptic, coating_brownian,
    substrate_brownian, substrate_thermoelastic,
)
from OptimalBragg import Material
from OptimalBragg.materials import SiN_123, aSi_123, cSi_123

plt.rcParams.update({'text.usetex': False,
                     'lines.linewidth': 4,
                     'font.family': 'serif',
                     'font.serif': 'Georgia',
                     'font.size': 22,
                     'xtick.direction': 'in',
                     'ytick.direction': 'in',
                     'xtick.labelsize': 'medium',
                     'ytick.labelsize': 'medium',
                     'axes.labelsize': 'small',
                     'axes.titlesize': 'medium',
                     'axes.grid.axis': 'both',
                     'axes.grid.which': 'both',
                     'axes.grid': True,
                     'grid.color': 'xkcd:Cerulean',
                     'grid.alpha': 0.2,
                     'lines.markersize': 12,
                     'legend.borderpad': 0.2,
                     'legend.fancybox': True,
                     'legend.fontsize': 'small',
                     'legend.framealpha': 0.8,
                     'legend.handletextpad': 0.5,
                     'legend.labelspacing': 0.33,
                     'legend.loc': 'best',
                     'figure.figsize': ((12, 8)),
                     'savefig.dpi': 140,
                     'savefig.bbox': 'tight',
                     'pdf.compression': 9})

# Absorption coefficients
alpha_SiN = 1e-3  # 1/m
alpha_aSi = 100e-6 / 1e-6  # 1/m

# Material properties
mat_sub = Material(cSi_123)
mat_low = Material(SiN_123)
mat_high = Material(aSi_123)
lambda0 = 2050e-9  # m

# Load most recent optimization result
fname = max(glob.iglob('Data/*Layers*.mat'), key=os.path.getctime)
fname = fname[5:]  # rm 'Data/' prefix
timestamp_str = 'date = ' + fname[-15:-9] + ' time = ' + fname[-8:-3]

print('Loading ' + fname + '...')
z = scio.loadmat('Data/' + fname, squeeze_me=True)
n = z['n']
T = z['T']

# Spectral reflectivity
lams = np.linspace(0.4, 1.6, 300)
rr, _ = multidiel1(n, z['L'], lams)
RR = np.abs(rr)**2
TT = 1 - RR

# Physical thicknesses
L = lambda0 * op2phys(z['L'], n[1:-1])

# E-field depth profile
Z, field = fieldDepth(L, n, pol='p', nPts=300, lam=lambda0)

layers = np.cumsum(1e6 * L)
layers = np.append(0, layers)

tic = default_timer()

# --- Spectral Transmission plot ---
fig, ax = plt.subplots(1, 1)
xx = 1e6 * lams * lambda0
ax.semilogy(xx, TT, lw=3, label='Transmissivity', c='xkcd:Red')
ax.semilogy(xx, RR, lw=3, label='Reflectivity', c='xkcd:electric blue', alpha=0.5)
ax.vlines(lambda0 * 1e6, T, 1, linestyle='--')
ax.set_xlabel(r'Wavelength [$\mu \mathrm{m}$]')
ax.set_ylabel('T or R')
ax.set_ylim((1e-6, 1))
ax.text(lambda0 * 1.25e6, 1e-3, f'T @ {1e9*lambda0:.1f} nm', size='x-small')
ax.text(lambda0 * 1.25e6, 0.5e-3, f'= {1e6*T:.1f} ppm', size='x-small')
ax.legend()
ax.text(1, 0.7, timestamp_str, fontsize=8, color='xkcd:Burple', alpha=0.95,
        rotation='-90', ha='left', va='center', transform=ax.transAxes)
fig.suptitle('a-Si:SiN coating')
plt.savefig('Figures/ETM_R.pdf')

# --- Layer structure + E-field plot ---
fig2, ax2 = plt.subplots(2, 1, sharex=True)
ax2[0].plot(Z * 1e6, field, color='xkcd:electric purple', alpha=0.97, rasterized=False)

absStr = f'$|\\vec E_{{\mathrm{{surface}}}}| = {1e6*field[0]:.0f}$ ppm of $\\vec |E_{{\mathrm{{inc}}}}|$'
absStr += '\n'
Esq = field  # field is already |E|^2 normalized
alphas = np.where(np.arange(len(L)) % 2 == 0, alpha_SiN, alpha_aSi)
intAbs = calc_abs(Esq, L, alphas)
absStr += f'Integrated absorption in stack is {1e6*intAbs:.1f} ppm'
print(f'Total integrated absorption {1e6*intAbs:.1f} ppm.')
ax2[0].text(0.5, 0.7, absStr, transform=ax2[0].transAxes, fontsize=14)

ax2[0].vlines(np.cumsum(L)[1:-1:2] * 1e6, 1e-5, 0.55,
              color='xkcd:bright teal', linewidth=0.6, linestyle='--', alpha=0.75)
ax2[0].vlines(np.cumsum(L)[::2] * 1e6, 1e-5, 0.55,
              color='xkcd:deep purple', linewidth=0.6, linestyle='--', alpha=0.75)

ax2[1].bar(layers[:-1:2], 1e9 * L[::2], width=1e6 * L[::2],
           align='edge', color='xkcd:bright teal', alpha=0.4, label=r'$\mathrm{SiN}$')
ax2[1].bar(layers[1:-1:2], 1e9 * L[1::2], width=1e6 * L[1::2],
           align='edge', color='xkcd:deep purple', alpha=0.4, label=r'$a-Si$')
ax2[1].legend()
ax2[1].yaxis.set_major_formatter(FormatStrFormatter("%3d"))
ax2[0].set_ylabel(r'Normalized $|E(z)|^2$')
ax2[1].set_ylabel('Physical layer thickness [nm]')
ax2[1].set_xlabel(r'Distance from air interface $[\mu \mathrm{m}]$')
ax2[0].text(1, 0.7, timestamp_str, fontsize=8, color='xkcd:Burple', alpha=0.95,
            rotation='-90', ha='left', va='center', transform=ax2[0].transAxes)
fig2.subplots_adjust(hspace=0.01, left=0.09, right=0.95, top=0.92)
fig2.suptitle('a-Si:SiN coating electric field')
plt.savefig('Figures/' + fname[:-4] + '.pdf')
plt.savefig('Figures/ETM_Layers.pdf')

# --- Thermal Noise plot ---
from OptimalBragg import qw_stack

Npairs = len(L) // 2
stack = qw_stack(
    lam_ref=lambda0,
    substrate=mat_sub,
    superstrate=Material({'Properties': {'n': 1.0}}),
    thin_films={'L': mat_low, 'H': mat_high},
    pattern='LH' * Npairs,
)
# Override with actual optimized thicknesses
stack['Ls_opt'] = z['L']
stack['Ls'] = L

w_beam = 0.062  # ETM beam radius [m]

ff = np.logspace(0, 4, 500)
fig3, ax3 = plt.subplots(1, 1)

StoZ = np.array([coating_thermooptic(np.array([f]), stack, w_beam, 0.17, 0.20)[0]
                 for f in ff])
SbrZ = np.array([coating_brownian(f, stack, w_beam) for f in ff])
subBrown = np.array([substrate_brownian(f, stack, w_beam) for f in ff])
subTE = np.array([substrate_thermoelastic(f, stack, w_beam) for f in ff])

ii = np.abs(ff - 100).argmin()
SbrZ100 = SbrZ[ii]

ax3.loglog(ff, np.sqrt(StoZ), label='Thermo-Optic', c='xkcd:Purplish Blue')
ax3.loglog(ff, np.sqrt(SbrZ), label='Brownian', c='xkcd:Tomato')
ax3.loglog(ff, np.sqrt(subBrown), label='Substrate Brownian', c='xkcd:Dusty Blue')
ax3.loglog(ff, np.sqrt(subTE), label='Substrate Thermo-Elastic',
           c='xkcd:Chocolate', alpha=0.3)
ax3.legend()
ax3.set_xlim([10, 10e3])
ax3.set_ylim([8e-24, 2e-20])

ax3.text(80, 11e-21, f'# of layers = {len(L)}', size='x-small')
ax3.text(80, 7e-21, f'Thickness = {1e6*sum(L):.2f} um', size='x-small')
ax3.text(50, 4.5e-21,
         f'$x_{{Brown}}$ @ 100 Hz = {np.sqrt(SbrZ100)*1e21:.2f} zm/$\\sqrt{{\\mathrm{{Hz}}}}$',
         size='x-small')
ax3.set_ylabel(r'Displacement Noise $[\mathrm{m} / \sqrt{\mathrm{Hz}}]$')
ax3.set_xlabel('Frequency [Hz]')
ax3.text(1, 0.7, timestamp_str, fontsize=8, color='xkcd:Burple', alpha=0.9,
         rotation='-90', ha='left', va='center', transform=ax3.transAxes)
fig3.suptitle('a-Si:SiN coating')
plt.savefig('Figures/ETM_TN.pdf')

dt = default_timer() - tic
print(f'Took {dt:.1f} sec to make the plots and save them.')
