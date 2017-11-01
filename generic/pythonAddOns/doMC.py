'''
Script to do some MC analysis on a given coating design
Example usage:

    python doMC.py aLIGO_ETM_layers_171026_1303.mat aLIGO_ETM_MC.hdf5 5000

The above function call will get the coating layer structure from
aLIGO_ETM_layers_171026_1303.mat, do an MC with 5000 data points, and save the
output of the MC to aLIGO_ETM_MC.hdf5
'''
import numpy as np
import scipy.io as sio
from scipy.interpolate import interp1d, PchipInterpolator
import emcee
import h5py
import sys
from mcutils import *
import gwinc 

#User config
matFileName = sys.argv[1]
hdfFileName = sys.argv[2]
nPts = int(sys.argv[3])

#Initialize the emcee sampler
# Emcee sampler object
nWalkers = 20
nDim = 4

means = np.zeros(nDim)   #The means will be 0.
cov = np.diag(0.005*np.ones(nDim)) #Gaussians with width of 0.5%
cov = np.dot(cov,cov)
icov = np.linalg.inv(cov)
p0 = np.random.rand(nDim * nWalkers).reshape((nWalkers, nDim))
#p0 = np.zeros((nWalkers, nDim))
sampler = emcee.EnsembleSampler(nWalkers, nDim, lnprob, args=[means,
    icov])
pos, prob, state = sampler.run_mcmc(p0, 1000)
sampler.reset()

pos_final, prob_final, state_final = sampler.run_mcmc(pos, 5000)

##Define a frequency vectory for thermal noise eval
ff = np.logspace(1,3,50)

#Use the generated perturbations to run the MCMC
nSamples = nPts
Tp_IR = np.ones(nSamples)
Rp_green = np.ones(nSamples)
surfField = np.ones(nSamples)
TOnoise = np.ones(len(ff))
Brnoise = np.ones(len(ff))
perturbs = sampler.flatchain[0:nSamples,:]

data = sio.loadmat(matFileName, struct_as_record=False, squeeze_me=True) #last two options for convenient access, see scipy.io help page for why,
n_IR_out = np.copy(data['costOut'].n_IR)
n_green_out = np.copy(data['costOut'].n_green)
L_out = np.copy(op2phys(data['costOut'].L, n_IR_out[1:-1]))
aoi_out = float(np.copy(data['costOut'].aoi))

#Now extract the "ifo" variable
ifo = data['ifo']
#set some of the material params in this structure to match that used to optimize coating, i.e. values from Ramin
#ifo.Materials.Coating.Indexhighn = data['costOut'].n_IR[2]
#ifo.Materials.Coating.Indexlown = data['costOut'].n_IR[1]
#ifo.Materials.Substrate.RefractiveIndex = 1.449641 #Corning datasheet


for jj in range(nSamples):
    aoi = np.copy(aoi_out)
    n_IRs = np.copy(n_IR_out)
    n_greens = np.copy(n_green_out)
    #Make it a fractional change
    perturb = 1+sampler.flatchain[jj,:]
    #AoI
    aoi *= (perturb[0])
    #Physical lengths
    Ls = np.copy(L_out)
    Ls *= (perturb[3])
    #SiO2 refractive indices
    n_IRs[1:-1:2] *= perturb[2]
    n_greens[1:-1:2] *= perturb[2]
    #Ta2O5 refractive indices
    n_IRs[2::2] *= perturb[1]
    n_greens[2::2] *= perturb[1]
    #Compute reflectivities
    [Gamma5p, t1] = multidiel1(n_IRs, Ls*n_IRs[1:-1], 1)
    Tp_IR[jj]=(1. - np.abs(Gamma5p)**2)
    surfField[jj]=surfaceField(Gamma5p)
    [Gamma5p, t2] = multidiel1(n_greens, Ls*n_greens[1:-1], 0.5)
    Rp_green[jj] = (np.abs(Gamma5p)**2)
    aa,bb,cc,dd = gwinc.noise.coatingthermal.getCoatThermoOptic(ff, ifo, 6e-2,np.array([Ls*n_IRs[1:-1]]).T)
    TOnoise = np.vstack((np.sqrt(aa),TOnoise))
    SbrZ = gwinc.noise.coatingthermal.getCoatBrownian(ff, ifo, 6e-2,np.array([Ls*n_IRs[1:-1]]).T)
    Brnoise = np.vstack((np.sqrt(SbrZ),Brnoise))

idx = np.argmin(np.abs(ff-100))
samples = np.vstack((1e6*Tp_IR,1e2*Rp_green,TOnoise[1:,idx]*1e21,Brnoise[1:,idx]*1e21,surfField))

#Save this thing into an hdf5 file.
f = h5py.File(hdfFileName,'w')
f['MCout'] = samples
f['TOnoise'] = TOnoise
f['Brnoise'] = Brnoise
f.close()
