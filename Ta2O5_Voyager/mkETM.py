r"""
run this code to optimize a layer structure design for a ETM HR coating
this coating is for the LIGO Voyager ETM operating at 123 K;
the low index material is SiO2 and the high index material is Ta2O5(tantala)
the center wavelength of the laser is 2128.2 nm, and the AUX wavelength
is 1418.8 nm
"""

from datetime import datetime
from timeit import default_timer

# To use pygwinc, first install; 
# >> conda install -c conda-forge gwinc
# import sys
# sys.path.append('../../pygwinc/')

from generic_local.optimUtils import *
from generic_local.coatingUtils import importParams

from scipy.optimize import differential_evolution as devo
from scipy.io import loadmat, savemat

import gwinc

def main(save=False):
    # Optimization targets, weights, and parameters
    params_file = 'ETM_params.yml'
    opt_params = importParams(params_file)

    # Load other IFO parameters
    ifo = gwinc.Struct.from_file(opt_params['gwincStructFile'])
    # Load noise budget
    voy = gwinc.load_budget('Voyager')
    # This is a fast approximator for Brownian noise
    gam = brownianProxy(ifo)

    # Set up stack
    Npairs = opt_params['Npairs']
    Ncopies = opt_params['Ncopies']
    Nfixed = opt_params['Nfixed']
    Nlayers = 2 * Npairs * Ncopies + 1
    
    # Initial guess
    Ls = 0.75*np.ones(2 * Npairs + 1)
    if __debug__:
        print(f"Shape of Ls array ={Ls.shape}")
 
    bounds = ((0.1, 0.45),)*(len(Ls) - 1)
    # Minimum thickness with 20 nm cap
    minThick = 20e-9 / ifo.Laser.Wavelength * 1.5
    # Make the first layer thin
    bounds = ((minThick, 0.45),) + bounds

    if __debug__:
        print(f'Bounds = {bounds[-1]}')
        tic = default_timer()
    # Do global optimization
    res = devo(func=getMirrorCost, 
            bounds=bounds, 
            updating='deferred',
            strategy='best1bin', 
            mutation=(0.1, 1.5),
            popsize=opt_params['Nparticles'], 
            workers=-1, 
            maxiter=2000,
            args=(params_file, ifo, gam, False, Ncopies, Nfixed),
            polish=True, 
            disp=True)

    if __debug__:
        print(" ")
        dt = default_timer() - tic
        print(f'Took {round(dt, 1)} sec to optimize.')
        print(" ")

    # Construct the optimum stack for saving / plotting
    Lres = res.x
    copiedLayers = np.tile(Lres[1:2*Npairs+1].copy(), Ncopies)
    Lres = np.append(Lres, copiedLayers)

    fixedLayers = np.tile(Lres[-2:].copy(), Nfixed)
    Lres = np.append(Lres, fixedLayers)

    # Run once to get the final costs
    scalarCost, costOut = getMirrorCost(Lres, 
                                        paramFile=params_file,
                                        ifo=ifo, 
                                        gam=gam, 
                                        verbose=True,
                                        copies=0, 
                                        fixed=0)

    if save:
        # make new dict for saving the data
        z = {}
        z["ifo_name"] = opt_params['gwincStructFile']
        z["opt_name"] = 'ETM'
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
        savemat('Data/ETM/' + fname, z, do_compression=True)

if __name__ == '__main__':
    main(save=True)