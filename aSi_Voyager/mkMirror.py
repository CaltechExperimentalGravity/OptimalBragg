r"""
run this code to optimize a layer structure design for a ETM HR coating

this coating is for the LIGO Voyager ETM operating at 123 K

the low index material is SiO2 and the high index material is a-Si

the center wavelength of the laser is 2128 nm (set in the aSiModel.m file)

"""

import sys
from datetime import datetime
from timeit import default_timer
#sys.path.append('../../pygwinc/')

sys.path.append('../generic/')
from optimUtils import *

from scipy.optimize import differential_evolution as diffevo
from scipy.io import loadmat,savemat

ifo = gwinc.Struct.from_file('aSiModel.m')
voy = gwinc.load_budget('Voyager')

# how many coating Layers?
Npairs = 8
Nlayers = 2*Npairs + 1
Ls = 0.25 * np.ones((Nlayers, 1)) # initial guess
Ls = Ls[:,0]
if __debug__:
    print("Shape of Ls array = " + str(np.shape(Ls)))
gam = brownianProxy(ifo)

#getMirrorCost(L=Ls, paramFile='params.yml',
#                  ifo=ifo, gam=gam, verbose=False)


# do Global Optimization
N_particles = 15

#x0 = np.random.uniform(0.05, 0.5, (N_particles, len(Ls)))
bow = ((0.1, 0.48),)
bounds = bow*(len(Ls)-1)
minThickCap = 20e-9 # min thickness of cap layer
minThick = minThickCap/ifo.Laser.Wavelength * 1.5
bounds = ((minThick, 0.6),) + bounds # make the first layer thin
if __debug__:
    print(np.shape(bounds))


tic = default_timer()
# minimize by Differential Evolution Optimizer
res = diffevo(func=getMirrorCost, bounds=bounds, updating='deferred',
                  strategy = 'best1bin', mutation = (0.1, 1.5),
                  popsize=N_particles, workers = -1,
                         args=('params.yml', ifo, gam, False),
                         polish=True, disp=True)

dt = default_timer() - tic
if __debug__:
    print(" ")
    print('Took ' + str(round(dt,1)) + ' sec to optimize.')
    print(" ")

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

fname   = 'Layers'
tnowstr = datetime.now().strftime('%y%m%d_%H%M')
fname   = z["opt_name"] + '_' + fname + '_' + tnowstr + '.mat'

z['filename'] = fname
# save layer data and also the whole ifo param file
savemat('Data/' + fname, z, do_compression=True)
