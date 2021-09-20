r"""
run this code to optimize a layer structure design for a ETM HR coating

this coating is for the LIGO Voyager ETM operating at 123 K

the low index material is SiO2 and the high index material is a-Si

the center wavelength of the laser is 2128 nm (set in the aSiModel.m file)

"""

# import some Python libraries
import sys
from datetime import datetime
from timeit import default_timer

from scipy.optimize import differential_evolution as devo
from scipy.io import loadmat,savemat

# install gwinc with anaconda: conda install -c conda-forge gwinc
#sys.path.append('../../pygwinc/')

sys.path.append('../generic/')
from optimUtils import *
from coatingUtils import importParams
from gwinc import Struct



coating_param_file = 'aSiModel.yaml'       # physical parameters of the mirrors
paramfilename = 'params.yml'               # optimization parameters
opt_params = importParams(paramfilename)

#ifo = gwinc.Struct.from_file(opt_params['gwincStructFile'])
voy = gwinc.load_budget('Voyager')
ifo = Struct.from_file(coating_param_file)

# how many coating Layers?
Npairs  = opt_params['Npairs']
Nlayers = 2*Npairs + 1
Ls      = 0.25 * np.ones((Nlayers, 1)) # initial guess
Ls      = Ls[:,0]
#if __debug__:
#    print("Shape of Ls array = " + str(np.shape(Ls)))

# this is an approximater for the Brownian noise which is fast to compute
brown_noise = brownianProxy(ifo)

# do Global Optimization
N_particles = opt_params['Nparticles']

#x0 = np.random.uniform(0.05, 0.5, (N_particles, len(Ls)))
bow = ((0.1, 0.49),)
bounds = bow*(len(Ls)-1)
minThickCap = 20e-9 # min thickness of cap layer
minThick = minThickCap/ifo.Laser.Wavelength * 1.5
bounds = ((minThick, 0.4),) + bounds # make the first layer thin
if __debug__:
#    print(np.shape(bounds))
    tic = default_timer()

the_strats = ['best1bin', 'best1exp', 'rand1exp', 'randtobest1exp', 'currenttobest1exp',
              'best2exp', 'rand2exp', 'randtobest1bin', 'currenttobest1bin', 'best2bin',
              'rand2bin', 'rand1bin']
mystrat = the_strats[np.random.randint(len(the_strats))]
mystrat = the_strats[0]  # use best1bin

# minimize by Differential Evolution Optimizer
res = devo(func=getMirrorCost, bounds=bounds, updating = 'deferred',
                  strategy = mystrat, mutation = (0.1, 1.5),
                  popsize=N_particles, workers = 1,
                         args=(paramfilename, ifo, brown_noise, False),
                         polish=True, disp=True)


if __debug__:
    print(" ")
    dt = default_timer() - tic
    print('Took ' + str(round(dt,1)) + ' sec to optimize with Strategy = ' + mystrat)
    print(" ")


# run once to get the costs for the final solution
vorb = True # False only returns 1 variable, so don't use that

scalarCost, costOut = getMirrorCost(L=res.x, paramFile=paramfilename,
              ifo=ifo, gam=brown_noise, verbose=vorb)



# make new dict for saving the data
z = {}
#z["result"]   = res
z["ifo_name"] = opt_params['gwincStructFile']
z["opt_name"] = 'ETM'

#        costOut = {}
#        costOut['n'] = n
z['L'] = res.x
z['T'] = costOut['T']
#        costOut['R'] = 1 - T
#        costOut['scalarCost'] = scalarCost
#        costOut['brownianProxy'] = cc2
#        costOut['vectorCost'] = cost

z["n"] = costOut['n']

fname   = 'Layers'
tnowstr = datetime.now().strftime('%y%m%d_%H%M')
fname   = z["opt_name"] + '_' + fname + '_' + tnowstr + '.mat'

z['filename'] = fname
# save layer data and also the whole ifo param file
savemat('Data/' + fname, z, do_compression=True)
