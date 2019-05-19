#!/usr/bin/env python

import sys
from datetime import datetime

sys.path.append('../../pygwinc/')

sys.path.append('../generic/')
from optimUtils import *

from scipy.optimize import differential_evolution as diffevo
from scipy.io import loadmat,savemat

ifo = gwinc.load_ifo('aSiModel.m')
voy = gwinc.load_ifo('Voyager')

# how many coating Layers?
Npairs = 8
Nlayers = 2*Npairs
Ls = 0.25 * np.ones((Nlayers, 1))
Ls = Ls[:,0]
gam = brownianProxy(ifo)

getMirrorCost(L=Ls, paramFile='params.yml',
                  ifo=ifo, gam=gam, verbose=False)


# do Global Optimization
N_particles = 15

#x0 = np.random.uniform(0.05, 0.5, (N_particles, len(Ls)))
bow = ((0.1, 0.5),)
bounds = bow*len(Ls-1)
minThickCap = 20e-9 # min thickness of cap layer
minThick = minThickCap/ifo.Laser.Wavelength * 1.5
bounds = ((minThick, 0.4),) + bounds # make the first layer thin

# minimize by Differential Evolution Optimizer
res = diffevo(func=getMirrorCost, bounds=bounds, updating='deferred',
                  popsize=N_particles, workers=-1,
                         args=('params.yml', ifo, gam, False),
                         polish=True, disp=True)

# run once to get the costs for the final solution
scalarCost, costOut = getMirrorCost(L=res.x, paramFile='params.yml',
              ifo=ifo, gam=gam, verbose=True)



# make new dict for saving the data
z = {}
#z["result"]   = res
z["ifo_name"] = 'aSiModel.m'
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

fname = 'Layers'
tnowstr = datetime.now().strftime('%y%m%d_%H%M')
fname = z["opt_name"] + '_' + fname + '_' + tnowstr

z['filename'] = fname
# save layer data and also the whole ifo param file
savemat('Data/' + fname, z, do_compression=True)
