r"""
run this code to optimize a layer structure design for a ITM HR coating

this coating is for the LIGO Voyager ITM operating at 123 K

the low index material is SiO2 and the high index material is a-Si

the center wavelength of the laser is 2128 nm

"""

import sys
from datetime import datetime
from timeit import default_timer

# install gwinc with anaconda: conda install -c conda-forge gwinc
#sys.path.append('../../pygwinc/')

from generic_local.optimUtils import *
from generic_local.coatingUtils import importParams

from scipy.optimize import differential_evolution as devo
from scipy.io import loadmat,savemat


paramfilename = 'ITM_params.yml'
opt_params = importParams(paramfilename)

ifo = gwinc.Struct.from_file(opt_params['gwincStructFile'])
voy = gwinc.load_budget('Voyager')

# how many coating Layers?
Npairs = opt_params['Npairs']
Nfixed = opt_params['Nfixed']
Nlayers = 2*Npairs + 1
Ls = 0.75*np.ones(Nlayers - 2*Nfixed) # initial guess
if __debug__:
    print("Shape of Ls array = " + str(np.shape(Ls)))

# this is an approximater for the Brownian noise which is fast to compute
gam = brownianProxy(ifo)

# do Global Optimization
N_particles = opt_params['Nparticles']

#x0 = np.random.uniform(0.05, 0.5, (N_particles, len(Ls)))
bow = ((0.1, 0.4),)
bounds = bow*(len(Ls) - 1)
minThickCap = 20e-9 # min thickness of cap layer
minThick = minThickCap/ifo.Laser.Wavelength * 1.5
bounds = ((minThick, 0.4),) + bounds # make the first layer thin
if __debug__:
    print(np.shape(bounds))
    tic = default_timer()


# minimize by Differential Evolution Optimizer
res = devo(func=getMirrorCost, bounds=bounds, updating='deferred',
                  strategy = 'best1bin', mutation = (0.05, 1.85),
                  popsize=N_particles, workers = -1, maxiter=2000,
                         args=(paramfilename, ifo, gam, False, Nfixed),
                         polish=True, disp=True)

if __debug__:
    print(" ")
    dt = default_timer() - tic
    print('Took ' + str(round(dt,1)) + ' sec to optimize.')
    print(" ")

# run once to get the costs for the final solution
scalarCost, costOut = getMirrorCost(L=res.x, paramFile=paramfilename,
              ifo=ifo, gam=gam, verbose=True, fixed=Nfixed)


# Manually add fixed layers to solution just for saving and plotting
Lres = res.x
fixedLayers = np.tile(Lres[-2:].copy(), Nfixed)
Lres = np.append(Lres, fixedLayers)


# make new dict for saving the data
z = {}
#z["result"]   = res
z["ifo_name"] = opt_params['gwincStructFile']
z["opt_name"] = 'ITM'

z['L'] = Lres
z['T'] = costOut['T']
z['Taux'] = costOut['Taux']
z['Toplv'] = costOut['Toplv']
z['scalarcost'] = scalarCost
z['vectorCost'] = costOut['vectorCost']

z["n"] = costOut['n']

fname   = 'Layers'
tnowstr = datetime.now().strftime('%y%m%d_%H%M')
fname   = z["opt_name"] + '_' + fname + '_' + tnowstr + '.mat'

z['filename'] = fname
# save layer data and also the whole ifo param file
savemat('Data/ITM/' + fname, z, do_compression=True)
