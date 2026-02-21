'''
Script to do some MC analysis on a given ITM coating design
Example usage:

    python doMC_ITM.py Data/ITM/ITM_Layers_240101_120000.hdf5 output_MC.hdf5 5000

The above function call will get the coating layer structure from
the HDF5 file, do an MC with 5000 data points, and save the
output of the MC to output_MC.hdf5
'''
import numpy as np
import emcee
import h5py
import sys
from generic.coatingUtils import *
from generic.coatingUtils import importParams
import gwinc
from gwinc import noise
import tqdm

# User config
hdf5FileName = sys.argv[1]
hdfFileName = sys.argv[2]
nPts = int(sys.argv[3])

# Initialize the emcee sampler
# Emcee sampler object
nWalkers = 64
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

# Frequency vector for thermal noise evaluation
ff = np.logspace(1, 3, 50)

# Use the generated perturbations to run the MCMC
nSamples = nPts
Tp_IR = np.ones(nSamples)
Tp_AUX = np.ones(nSamples)
surfField = np.ones(nSamples)
TOnoise = np.ones(len(ff))
Brnoise = np.ones(len(ff))
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

L_out = op2phys(L_opt, n_IR_out[1:-1])

# Now extract the "ifo" variable
ifo = gwinc.Struct.from_file(gwincFile)

# AUX wavelength ratio from params YAML
opt_params = importParams('ITM_params.yml')
lambdaAUX = opt_params['misc']['lambdaAUX']

# Pre-build all perturbed arrays at once (vectorized, no per-iteration copies)
perturb_factors = 1 + perturbs  # (nSamples, 4)

# Pre-build all perturbed n arrays: (nSamples, len(n_IR_out))
n_IR_all = np.tile(n_IR_out, (nSamples, 1))
n_IR_all[:, 1:-1:2] *= perturb_factors[:, 2:3]   # low-n perturbation
n_IR_all[:, 2::2] *= perturb_factors[:, 1:2]      # high-n perturbation

# Pre-build all perturbed L arrays: (nSamples, len(L_out))
L_all = L_out[None, :] * perturb_factors[:, 3:4]  # thickness perturbation

for jj in tqdm.tqdm(range(nSamples)):
    n_IRs = n_IR_all[jj]
    Ls = L_all[jj]
    # Compute reflectivities
    [Gamma5p, t1] = multidiel1(n_IRs, Ls * n_IRs[1:-1], 1)
    Tp_IR[jj] = (1. - np.abs(Gamma5p)**2)
    surfField[jj] = surfaceField(Gamma5p)
    [Gamma5p, t2] = multidiel1(n_IRs, Ls * n_IRs[1:-1], lambdaAUX)
    Tp_AUX[jj] = (1 - np.abs(Gamma5p)**2)
    # Compute thermal noise for this perturbed stack
    # Use ETM for coating/substrate properties, ITM for beam radius
    mir = ifo.Optics.ETM
    mir.Coating.dOpt = Ls * n_IRs[1:-1]
    aa, bb, cc, dd = noise.coatingthermal.coating_thermooptic(
        ff, mir, ifo.Laser.Wavelength, ifo.Optics.ITM.BeamRadius)
    SbrZ = noise.coatingthermal.coating_brownian(
        ff, mir, ifo.Laser.Wavelength, ifo.Optics.ITM.BeamRadius)
    TOnoise = np.vstack((np.sqrt(aa), TOnoise))
    Brnoise = np.vstack((np.sqrt(SbrZ), Brnoise))

idx = np.argmin(np.abs(ff - 100))
samples = np.vstack((1e6*Tp_IR, 1e2*Tp_AUX,
                      TOnoise[:-1, idx]*1e21, Brnoise[:-1, idx]*1e21,
                      surfField))

# Save this thing into an hdf5 file.
with h5py.File(hdfFileName, 'w') as f:
    f['MCout'] = samples
    f['TOnoise'] = TOnoise
    f['Brnoise'] = Brnoise
