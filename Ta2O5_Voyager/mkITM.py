r"""
run this code to optimize a layer structure design for a ITM HR coating
this coating is for the LIGO Voyager ITM operating at 123 K;
the low index material is SiO2 and the high index material is Ta2O5(tantala)
the center wavelength of the laser is 2128.2 nm, and the AUX wavelength
is 1418.8 nm
"""

from datetime import datetime
from timeit import default_timer

from scipy.optimize import differential_evolution as devo
import h5py

# To use pygwinc, first install; 
# >> conda install -c conda-forge gwinc
# import sys
# sys.path.append('../../pygwinc/')
import gwinc

from generic_local.optimUtils import *
from generic_local.coatingUtils import importParams

def main(save=False):
    # Optimization targets, weights, and parameters
    opt_params = importParams('ITM_params.yml')

    # Load other IFO parameters
    ifo = gwinc.Struct.from_file(opt_params['misc']['gwincStructFile'])
    # Load noise budget
    voy = gwinc.load_budget('Voyager')
    # This is a fast approximator for Brownian noise
    gam = brownianProxy(ifo)

    # Initial guess
    Ls = 0.75*np.ones(2 * opt_params['misc']['Npairs'] + 1)
    if __debug__:
        print(f"Shape of Ls array ={Ls.shape}")

    bounds = ((0.05, 0.45),)*(len(Ls) - 1)
    # Minimum thickness with 20 nm cap
    minThick = 20e-9 / ifo.Laser.Wavelength * 1.5
    # Make the first layer thin
    bounds = ((minThick, 0.45),) + bounds

    if __debug__:
        print(f'Bounds = {bounds[-1]}')
        tic = default_timer()

    vector_mon, conv_mon = [], []
    def diffevo_monitor(xk, convergence):
        vector_mon.append(xk)
        conv_mon.append(1/convergence)
        return False

    # Do global optimization
    res = devo(func=getMirrorCost, 
            bounds=bounds, 
            updating='deferred',
            strategy='best1bin', 
            mutation=(0.05, 1.5),
            popsize=opt_params['misc']['Nparticles'], 
            workers=-1, 
            maxiter=2000,
            args=(opt_params['costs'], ifo, gam, False, opt_params['misc']),
            polish=True, 
            callback=diffevo_monitor,
            disp=True)

    if __debug__:
        print(" ")
        dt = default_timer() - tic
        print(f'Took {round(dt, 1)} sec to optimize.')
        print(" ")

    # Construct the optimum stack for saving / plotting
    Lres = res.x
    copiedLayers = np.tile(Lres[1:2*opt_params['misc']['Npairs']+1].copy(), 
                            opt_params['misc']['Ncopies'])
    Lres = np.append(Lres, copiedLayers)

    fixedLayers = np.tile(Lres[-2:].copy(), opt_params['misc']['Nfixed'])
    Lres = np.append(Lres, fixedLayers)

    # Run once to get the final costs (but don't use Ncopies, or Nfixed!)
    final_misc = opt_params['misc']
    final_misc.update({'Ncopies':0, 'Nfixed':0})
    scalar_cost, output = getMirrorCost(Lres, 
                                        costs=opt_params['costs'],
                                        ifo=ifo, 
                                        gam=gam, 
                                        verbose=True,
                                        misc=final_misc)
    if save:
        tnowstr = datetime.now().strftime('%y%m%d_%H%M%S')
        fname = 'Data/ITM/ITM_Layers_' + tnowstr + '.hdf5'

        with h5py.File(fname, 'w') as f:
            main_group = 'diffevo_output'
            f.create_group(main_group)
            # How did the optimizer find its optimum?
            f.create_dataset('trajectory', data=np.array(conv_mon))
            f.create_dataset('vec_evol', data=np.array(vector_mon))
            vector_cost = output.pop('vectorCost')
            for cost, optimum in vector_cost.items():
                f.create_dataset(main_group + '/vectorCost/' + cost, data=optimum)
            for key, out in output.items():
                f.create_dataset(main_group + '/' + key, data=out)

if __name__ == '__main__':
    main(save=True)