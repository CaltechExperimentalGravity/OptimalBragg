# Set of functions used to evaluate a cost function for optimization

import numpy as np
from .coatingUtils import *
import os
import gwinc.noise

def transmissionCost(target, n, L, lamb=1, theta=0, pol='te'):
    ''' 
    Function that evaluates the transmission of the coating specified in
    coat, and evaluates a cost based on how close/far it is to the target value

    Parameters:
    -------------
    target: float
        Target transmission
    n: array_like
        Array of refractive indices, including the incident
        and transmitted media. Ordered from incident medium to
        transmitted medium.
    L: array_like
        Array of optical thicknesses comprising the dielectric stack,
        ordered from incident medium to transmitted medium.
        Should have 2 fewer elements than n.
    lamb: float or array_like
        Wavelength(s) at which the reflectivity is to be evaluated,
        in units of some central (design) wavelength.
    theta: float
        Angle of incidence in degrees. Defaults to 0 degrees (normal incidence)
    pol: str, 'te' or 'tm'
        Polarization at which reflectivity is to be evaluated.
        Defaults to 'te' (s-polarization)
    Returns:
    ---------
    cost: array_like
        An array of scalar costs for the transmission
    T: float
        Transmission of the coating.
    '''
    r, _ = multidiel1(n, L, lamb, theta, pol)
    T = 1 - np.abs(r)**2
    return (np.abs((target - T)/target)**2)[0], T


def sensitivityCost(target, n, L, lamb=1, theta=0, pol='te'):
    ''' 
    Function that evaluates the sensitivity cost from a target 
    transmission relative to 1% thickness perturbation

    Parameters:
    -------------
    target: float
        Target transmission
    n: array_like
        Array of refractive indices, including the incident
        and transmitted media. Ordered from incident medium to
        transmitted medium.
    L: array_like
        Array of optical thicknesses comprising the dielectric stack,
        ordered from incident medium to transmitted medium.
        Should have 2 fewer elements than n.
    lamb: float or array_like
        Wavelength(s) at which the reflectivity is to be evaluated,
        in units of some central (design) wavelength.
    theta: float
        Angle of incidence in degrees. Defaults to 0 degrees (normal incidence)
    pol: str, 'te' or 'tm'
        Polarization at which reflectivity is to be evaluated.
        Defaults to 'te' (s-polarization)
    Returns:
    ---------
    cost: array_like
        An array of scalar costs for the transmission
    '''
    return transmissionCost(target, n, 1.01*L, lamb, theta, pol)[0]


def surfEfieldCost(target, n, L, lamb=1, theta=0, pol='te'):
    ''' 
    Function that evaluates the surface E field cost

    Parameters:
    -------------
    target: float
        Target transmission
    n: array_like
        Array of refractive indices, including the incident
        and transmitted media. Ordered from incident medium to
        transmitted medium.
    L: array_like
        Array of optical thicknesses comprising the dielectric stack,
        ordered from incident medium to transmitted medium.
        Should have 2 fewer elements than n.
    lamb: float or array_like
        Wavelength(s) at which the reflectivity is to be evaluated,
        in units of some central (design) wavelength.
    theta: float
        Angle of incidence in degrees. Defaults to 0 degrees (normal incidence)
    pol: str, 'te' or 'tm'
        Polarization at which reflectivity is to be evaluated.
        Defaults to 'te' (s-polarization)
    Returns:
    ---------
    cost: array_like
        An array of scalar costs for the transmission
    '''
    r, _ = multidiel1(n, L, lamb, theta, pol)
    return (50 * np.arcsinh(np.abs(1 + r)**2))[0]


def stdevLCost(target, L,):
    """Get cost of relative variation of thicknesses in the stack
    
    Args:
        target (float): Target standard deviation
        L: array_like
        Array of optical thicknesses comprising the dielectric stack,
        ordered from incident medium to transmitted medium.
        Should have 2 fewer elements than n.
    Returns:
        cost: array_like
            An array of scalar costs for the relative thickness variance
    """
    if np.std(np.array(L)):
        relative_stdev = np.mean(np.array(L)) / np.std(np.array(L))
    else:
        relative_stdev = 0.0
    return np.abs((target - relative_stdev) / target)**2


def brownianProxy(ifo):
    '''
    Evaluate a bunch of material properties for speedy evaluation of Brownian noise.
    Parameters:
    -----------
    ifo: pygwinc struct
        Return object of pygwinc load_ifo function, which has
        all the material properties
    Returns:
    --------
    gam: float
        The pre-factor for use in the proxy function for Brownian noise
        per E0900068 pg4.
    '''
    phi_high = ifo.Materials.Coating.Phihighn
    n_high = ifo.Materials.Coating.Indexhighn
    Y_high = ifo.Materials.Coating.Yhighn
    phi_low = ifo.Materials.Coating.Philown
    n_low = ifo.Materials.Coating.Indexlown
    Y_low = ifo.Materials.Coating.Ylown
    Y_sub = ifo.Materials.Substrate.MirrorY
    a = phi_high / phi_low
    b = n_low / n_high
    c = Y_high / Y_sub + Y_sub / Y_high
    d = Y_low / Y_sub + Y_sub / Y_low
    gam = a*b*c/d
    return(gam)


def brownianCost(target, L, gam,):
    '''
    Calculate a proxy for the brownian noise for a coating specified by L,
    using the parameters specified in ifoModel which is a matlab struct
    or pygwinc yaml file. Formula taken from E0900068 pg4.
    Parameters:
    -----------
    target: float
        Target value for this cost function.
    L: array_like
        Array of optical thicknesses comprising the dielectric stack,
        ordered from incident medium to transmitted medium.
    gam:float
        The pre-factor for use in the proxy function for Brownian noise
        per E0900068 pg4.
    Returns:
    ---------
    cost: float
        Scalar cost which is a proxy for the Brownian noise level
    '''
    zLow = np.sum(L[::2])    # Sum of thicknesses of low index layers
    zHigh = np.sum(L[1::2])   # Sum of thicknesses of high index layers
    SBrZ = zLow + gam*zHigh  # Proxy brownian noise
    return target * SBrZ


def thermoopticCost(target, fTarget, L, ifo):
    '''
    Function to calculate a cost for the Thermo-Optic noise for
    a coating specified by L.
    pygwinc is used to do the TO noise evaluation.

    Parameters:
    -----------
    target: float
        Target value for this cost function
    L: array_like
        Array of optical thicknesses comprising the dielectric stack,
        ordered from incident medium to transmitted medium.
    fTarget: float
        Frequency at which to evaluate TO noise.
    ifo: pygwinc struct
        Return object of pygwinc load_ifo function, which has
        all the material properties
    Returns:
    ---------
    cost: float
        Scalar cost for the TO noise.
    '''
    # Get the TO noise PSD
    # Build up a "mirror" structure as required by pygwinc
    mir = ifo.Optics.ETM
    mir.Coating.dOpt = L
    StoZ, _, _, _ = gwinc.noise.coatingthermal.coating_thermooptic(
        fTarget, mir, ifo.Laser.Wavelength, ifo.Optics.ETM.BeamRadius)
    return target * StoZ


def getMirrorCost(L, costs, ifo, gam, verbose=False, misc={}):
    '''
    Compute the cost function for a coating design specified by L,
    based on the settings in paramFile.

    Parameters:
    ------------
    costs: dict
        Dictionary with scalar cost keys, and {'target':t, 'weight':w} vals
        e.g. {'TransPSL': {'target':500e-6, 'weight': 1.0},
              'TransAUX': {'target':100e-6, 'weight': 0.1},
              'Esurf':    {'target':0,      'weight': 0.0}}
    L: array_like
        Array of optical thicknesses comprising the dielectric stack,
        ordered from incident medium to transmitted medium.
    ifo: pygwinc struct
        Return object of pygwinc load_ifo function, which has
        all the material properties
    gam:
        pre-computed value of little gamma used as a proxy Brownian
        noise cost
    verbose: bool
        Determines level of detail returned by the function.
        Defaults to False, which outputs only the value of the cost function.
    **misc: keyword arguments parsed directly from the 'misc' parameters
    Returns:
    ----------
    scalarCost: float
        Scalar cost function (obtained by multiplying vector cost by weights)
    costOut: dict
        If verbose=True, several sub-properties of the coating are supplied.
    '''

    # Add Ncopies of single variable stack
    copiedLayers = np.tile(L[:-2].copy(), misc['Ncopies'])
    L = np.append(L, copiedLayers)

    # Add fixed layers at end? Default to 0
    fixedLayers = np.tile(L[-2:].copy(), misc['Nfixed'])
    L = np.append(L, fixedLayers)

    # Build up the array of refractive indices
    doublet = np.tile(np.array([ifo.Materials.Coating.Indexlown,
                                ifo.Materials.Coating.Indexhighn]),
                          int(np.floor(len(L)/2)))

    if len(doublet) != len(L):
        # Add another low index layer at the bottom of the stack
        doublet = np.append(doublet, doublet[0])

    # Add air, fixed extension, and substrate
    n = np.append(1, doublet)
    n = np.append(n, ifo.Materials.Substrate.RefractiveIndex)

    # Iterate over scalar costs
    vector_cost, output = {}, {}
    scalar_cost = 0.0
    for cost, specs in costs.items():
        if specs['weight']:
            if cost == 'TransPSL':
                vector_cost[cost], TPSL = transmissionCost(specs['target'],
                                                    n, L, 1.0, 
                                                    misc['aoi'], misc['pol'])
                if verbose:
                    output['TPSL'] = TPSL
                    output['RPSL'] = 1 - TPSL
            if cost == 'TransAUX':
                vector_cost[cost], TAUX = transmissionCost(specs['target'],
                                                    n, L, 2/3, 
                                                    misc['aoi'], misc['pol'])
                if verbose:
                    output['TAUX'] = TAUX
                    output['RAUX'] = 1 - TAUX
            if cost == 'TransOPLEV':
                vector_cost[cost], TOPL = transmissionCost(specs['target'],
                                                    n, L, 0.297,
                                                    misc['aoi'], misc['pol'])
                if verbose:
                    output['TOPL'] = TOPL
                    output['ROPL'] = 1 - TOPL
            if cost == 'Brownian':
                vector_cost[cost] = brownianCost(specs['target'], L, gam)
            if cost == 'Thermooptic':
                vector_cost[cost] = thermoopticCost(specs['target'], misc['fTO'],
                                                   L, ifo)
            if cost == 'Absorption':
                # NOT implemented yet
                pass
            if cost == 'Lstdev':
                vector_cost[cost] = stdevLCost(specs['target'], L)
            if cost == 'Esurf':
                vector_cost[cost] = surfEfieldCost(specs['target'],
                                                    n, L, 1.0, 
                                                    misc['aoi'], misc['pol'])
            if cost == 'Lsens':
                # Quadrature sum of PSL and AUX sensitivities
                sensPSL = sensitivityCost(specs['target'],
                                                    n, L, 1.0, 
                                                    misc['aoi'], misc['pol'])
                sensAUX = sensitivityCost(specs['target'],
                                                    n, L, 2/3, 
                                                    misc['aoi'], misc['pol'])
                vector_cost[cost] = np.sqrt(sensPSL ** 2 + sensAUX ** 2)
            # Make weighted scalar cost
            scalar_cost += specs['weight'] * vector_cost[cost]
    if verbose:
        for cost in vector_cost.keys():
            print(cost + f' cost = {vector_cost[cost]:.4f}')
        output['n'] = n
        output['L'] = L
        output['scalarCost'] = scalar_cost
        output['vectorCost'] = vector_cost
        return scalar_cost, output
    else:
        return scalar_cost
