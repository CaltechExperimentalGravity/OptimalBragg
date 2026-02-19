'''
Script to do some MC analysis on a given coating design
Example usage:

    python doMC.py Data/ETM/ETM_Layers_240101_120000.hdf5 output_MC.hdf5 5000

The above function call will get the coating layer structure from
the HDF5 file, do an MC with 5000 data points, and save the
output of the MC to output_MC.hdf5
'''
import numpy as np
import emcee
import h5py
import sys
from generic.coatingUtils import *
import gwinc
import tqdm

# User config
hdf5FileName = sys.argv[1]
hdfFileName = sys.argv[2]
nPts = int(sys.argv[3])

# Initialize the emcee sampler
# Emcee sampler object
nWalkers = 20
nDim = 4

means = np.zeros(nDim)   # The means will be 0.
cov = np.diag(0.005*np.ones(nDim))  # Gaussians with width of 0.5%
cov = np.dot(cov, cov)
icov = np.linalg.inv(cov)
p0 = np.random.rand(nDim * nWalkers).reshape((nWalkers, nDim))
sampler = emcee.EnsembleSampler(nWalkers, nDim, lnprob, args=[means,
    icov])
pos, prob, state = sampler.run_mcmc(p0, 1000)
sampler.reset()

pos_final, prob_final, state_final = sampler.run_mcmc(pos, 5000)

# Define a frequency vector for thermal noise eval
ff = np.logspace(1, 3, 50)

# Use the generated perturbations to run the MCMC
nSamples = nPts
Tp_IR = np.ones(nSamples)
Tp_AUX = np.ones(nSamples)
surfField = np.ones(nSamples)
perturbs = sampler.flatchain[0:nSamples, :]

# Read optimizer output from HDF5
with h5py.File(hdf5FileName, 'r') as f:
    n_IR_out = np.array(f['diffevo_output/n'])
    L_opt = np.array(f['diffevo_output/L'])
    gwincFile = f['gwincStructFile'][()]
    if isinstance(gwincFile, bytes):
        gwincFile = gwincFile.decode()
    else:
        gwincFile = str(gwincFile)

n_green_out = np.copy(n_IR_out)
L_out = op2phys(L_opt, n_IR_out[1:-1])

# Now extract the "ifo" variable
ifo = gwinc.Struct.from_file(gwincFile)

# AUX wavelength ratio — read from params if available, else default
lambdaAUX = 2/3  # default; override from params if needed

for jj in tqdm.tqdm(range(nSamples)):
    aoi = 0
    n_IRs = np.copy(n_IR_out)
    n_greens = np.copy(n_green_out)
    # Make it a fractional change
    perturb = 1 + sampler.flatchain[jj, :]
    # Physical lengths
    Ls = np.copy(L_out)
    Ls *= (perturb[3])
    # Low-n refractive indices
    n_IRs[1:-1:2] *= perturb[2]
    n_greens[1:-1:2] *= perturb[2]
    # High-n refractive indices
    n_IRs[2::2] *= perturb[1]
    n_greens[2::2] *= perturb[1]
    # Compute reflectivities
    [Gamma5p, t1] = multidiel1(n_IRs, Ls*n_IRs[1:-1], 1)
    Tp_IR[jj] = (1. - np.abs(Gamma5p)**2)
    surfField[jj] = surfaceField(Gamma5p)
    [Gamma5p, t2] = multidiel1(n_greens, Ls*n_greens[1:-1], lambdaAUX)
    Tp_AUX[jj] = (1 - np.abs(Gamma5p)**2)

idx = np.argmin(np.abs(ff - 100))
samples = np.vstack((1e6*Tp_IR, 1e6*Tp_AUX))

# Save this thing into an hdf5 file.
with h5py.File(hdfFileName, 'w') as f:
    f['MCout'] = samples
