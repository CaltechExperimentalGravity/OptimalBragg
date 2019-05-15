# Set of functions used to evaluate a cost function for optimization

from coatingUtils import *
import os, logging
import matplotlib.pyplot as plt, matplotlib as mpl

# Logging setup
logging.basicConfig(
    level=os.getenv('LOG_LEVEL', 'INFO'),
    format="%(levelname)s \n%(message)s")


# Some matplotlib setup
mpl.style.use('fivethirtyeight')
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['font.size'] = 22
mpl.rcParams['axes.facecolor']  = 'white'
mpl.rcParams['axes.edgecolor']  = 'black'
mpl.rcParams['axes.titlesize']  = 'x-large'
mpl.rcParams['axes.labelsize']  = 'x-large'
mpl.rcParams['xtick.labelsize'] = 'x-large'
mpl.rcParams['ytick.labelsize'] = 'x-large'
mpl.rcParams['axes.formatter.limits'] = [-2,2]
mpl.rcParams['axes.grid'] = True
mpl.rcParams['grid.linestyle'] = '--'
mpl.rcParams['grid.linewidth'] = 0.7
mpl.rcParams['grid.alpha'] = 0.4
mpl.rcParams['text.usetex'] = False
mpl.rcParams['legend.loc'] = 'best'
mpl.rcParams['legend.fontsize'] = 'small'
mpl.rcParams['figure.figsize'] = [12,9]
mpl.rcParams['figure.facecolor'] = 'white'
mpl.rcParams['figure.subplot.left'] = 0.07
mpl.rcParams['figure.subplot.right'] = 0.95
mpl.rcParams['figure.subplot.top'] = 0.92


def transmissionCost(target, n, L, lamb=1, theta=0, pol='te', sensL=False, surfE=True):
    '''
    Function that evaluates the transmission of the coating specified in
    coat, and evaluates a cost based on how close/far it is to the target value

    Parameters:
    -------------
    target: float
        Target transmission
    n: array_like
        Array of refractive indices, including the incident and transmitted media. Ordered from incident medium to transmitted medium.
    L: array_like
        Array of optical thicknesses comprising the dielectric stack, ordered from incident medium to transmitted medium.. Should have 2 fewer elements than n.
    lamb: float or array_like
        Wavelength(s) at which the reflectivity is to be evaluated, in units of some central (design) wavelength.
    theta: float
        Angle of incidence in degrees. Defaults to 0 degrees (normal incidence)
    pol: str, 'te' or 'tm'
        Polarization at which reflectivity is to be evaluated. Defaults to 'te' (s-polarization)
    sensL: bool
        Option to compute sensitivity cost
    surfE: bool
        Option to compute cost for surface E field.

    Returns:
    ---------
    cost: array_like
        An array of scalar costs for the transmission, and if requested, the sensitivity and the surface E field
    T: float
        Transmission of the coating.
    '''
    r, _ = multidiel1(n, L, lamb, theta, pol)
    T = 1 - np.abs(r)**2
    costT = np.abs((target - T)/target)**2      # The functional form may be easily changed
    cost = costT
    if sensL:
        rPerturb, _ = multidiel1(n, 1.01*L, lamb, theta, pol)
        TPerturb = 1 - np.abs(rPerturb)**2
        costS = np.abs((T-TPerturb)/0.01/T)     # Fractional change in the transmissivity, relative to 1%
        cost = np.append(cost, costS)
    if surfE:
        costE = 50 * np.arcsinh(np.abs(1 + r)**2)
        cost = np.append(cost, costE)
    return(cost, T)

def brownianProxy(ifo):
    '''
    Evaluate a bunch of material properties for speedy evaluation of Brownian noise.
    Parameters:
    -----------
    ifo: pygwinc struct
        Return object of pygwinc load_ifo function, which has all the material properties
    Returns:
    --------
    gam: float
        The pre-factor for use in the proxy function for Brownian noise per E0900068 pg4.
    '''
    phi_high = ifo.Materials.Coating.Phihighn
    n_high   = ifo.Materials.Coating.Indexhighn
    Y_high   = ifo.Materials.Coating.Yhighn
    phi_low  = ifo.Materials.Coating.Philown
    n_low    = ifo.Materials.Coating.Indexlown
    Y_low    = ifo.Materials.Coating.Ylown
    Y_sub    = ifo.Materials.Substrate.MirrorY
    a = phi_high / phi_low
    b = n_low / n_high
    c = Y_high / Y_sub + Y_sub / Y_high
    d = Y_low / Y_sub  + Y_sub / Y_low
    gam = a*b*c/d
    return(gam)

def brownianCost(target, L, gam):
    '''
    Calculate a proxy for the brownian noise for a coating specified by L,
    using the parameters specified in ifoModel which is a matlab struct
    or pygwinc yaml file. Formula taken from E0900068 pg4.
    Parameters:
    -----------
    target: float
        Target value for this cost function.
    L: array_like
        Array of optical thicknesses comprising the dielectric stack, ordered from incident medium to transmitted medium.
    gam:float
        The pre-factor for use in the proxy function for Brownian noise per E0900068 pg4.
    Returns:
    ---------
    cost: float
        Scalar cost which is a proxy for the Brownian noise level
    '''
    zLow = np.sum(L[::2])    # Sum of thicknesses of low index layers
    zHigh = np.sum(L[1::2])  # Sum of thicknesses of high index layers
    SBrZ = zLow + gam*zHigh  # Proxy brownian noise
    cost = target * SBrZ     # Functional form of the cost may be easily changed
    return(cost)

def TOcost(target, L, fTarget, ifo):
    '''
    Function to calculate a cost for the Thermo-Optic noise for a coating specified by L.
    pygwinc is used to do the TO noise evaluation.

    Parameters:
    -----------
    target: float
        Target value for this cost function
    L: array_like
        Array of optical thicknesses comprising the dielectric stack, ordered from incident medium to transmitted medium.
    fTarget: float
        Frequency at which to evaluate TO noise.
    ifo: pygwinc struct
        Return object of pygwinc load_ifo function, which has all the material properties
    Returns:
    ---------
    cost: float
        Scalar cost for the TO noise.
    '''
    # Get the TO noise PSD
    StoZ, _, _, _ = gwinc.noise.coatingthermal.getCoatThermoOptic(fTarget, ifo, ifo.Optics.ETM.BeamRadius, L)
    cost = target * StoZ  # Functional form of the cost may be easily changed
    return(cost)

def getMirrorCost(paramFile, L, ifo, gam, verbose=False):
    '''
    Compute the cost function for a coating design specified by L, based on the settings in paramFile.

    Parameters:
    ------------
    paramFile: str
        Path to the parameter file.
    L: array_like
        Array of optical thicknesses comprising the dielectric stack, ordered from incident medium to transmitted medium.
    ifo: pygwinc struct
        Return object of pygwinc load_ifo function, which has all the material properties
    gam:
        pre-computed value of little gamma used as a proxy Brownian noise cost
    verbose: bool
        Determines level of detail returned by the function. Defaults to False, which outputs only the value of the cost function.

    Returns:
    ----------
    scalarCost: float
        Scalar cost function (obtained by multiplying vector cost by weights)
    costOut: dict
        If verbose=True, several sub-properties of the coating are supplied.
    '''
    # Load the parameters
    par = importParams(paramFile)
    if len(par['costs']) != len(par['weights']):
        logging.critical('Parameter file is not configured correctly. Please check it.')
        return()
    # Build up the array of refractive indices
    doublet = np.tile(np.array([ifo.Materials.Coating.Indexlown, ifo.Materials.Coating.Indexhighn]), round(len(L)/2))
    if len(doublet) != len(L):
        # Add another low index layer at the bottom of the stack
        doublet = np.append(doublet, doublet[0])
    # Add air and substrate
    n = np.append(1, doublet)
    n = np.append(n, ifo.Materials.Substrate.RefractiveIndex)
    # This is set up in this way just for avoiding computing extra costs if so desired. 
    # In general, the weight vector can be used to select only the desired costs by setting the others to 0
    if 'Trans' in par['costs']:
        if 'sensL' in par['costs']:
            if 'surfE' in par['costs']:
                cc1, T = transmissionCost(par['targets'][0], n, L, lamb=1, theta=par['aoi'], pol=par['pol'], sensL=True, surfE=True)
            else:
                cc1, T = transmissionCost(par['targets'][0], n, L, lamb=1, theta=par['aoi'], pol=par['pol'], sensL=True, surfE=False)
        else:
            cc1, T = transmissionCost(par['targets'][0], n, L, lamb=1, theta=par['aoi'], pol=par['pol'], sensL=False, surfE=False)
    if 'coatBr' in par['costs']:
        cc2 = brownianCost(par['targets'][1], L, gam)
    if 'coatTO' in par['costs']:
        cc3 = TOcost(par['targets'][2], L, par['fTO'], ifo)
    # Make the cost
    cost = np.array([cc1[0], cc2, cc3, cc1[1], cc1[2]])
    scalarCost = np.dot(np.array(par['weights']), cost)
    if verbose:
        logging.critical('Cost for Brownian noise      =  {}'.format(round(cost[1],3)))
        logging.critical('Cost for Transmission        =  {}'.format(round(cost[0],3)))
        logging.critical('Cost for Thermo-Optic noise  =  {}'.format(round(cost[2],3)))
        logging.critical('Cost for sensitivity (dT/dL) =  {}'.format(round(cost[3],3)))
        logging.critical('Cost for surface E field     =  {}'.format(round(cost[4],3)))
        costOut = {}
        costOut['n'] = n
        costOut['L'] = L
        costOut['T'] = T
        costOut['R'] = 1. - T
        costOut['scalarCost'] = scalarCost
        costOut['brownianProxy'] = cc2
        costOut['vectorCost'] = cost
        
        # Optionally also make a plot of the spectral reflectivity
        lambda0 = ifo.Laser.Wavelength # meters
        lams = np.linspace(0.4,1.6,501)
        rr, _ = multidiel1(n, L, np.linspace(0.4,1.6,501), par['aoi'], par['pol'])
        RR = np.abs(rr)**2
        TT = 1 - np.abs(rr)**2
        fig, ax = plt.subplots(1,1,figsize=(16,9))
        ax.semilogy(1e6*lams*ifo.Laser.Wavelength, TT, label='Transmissivity')
        ax.semilogy(1e6*lams*ifo.Laser.Wavelength, RR, label='Reflectivity')
        ax.vlines(ifo.Laser.Wavelength*1e6, T/2, 1, linestyle='--')
        ax.set_xlabel('Wavelength [$\mu \\mathrm{m}$]')
        ax.set_ylabel('T or R')
        ax.legend(loc='best')
        logging.critical('Transmission of this coating at {} um is {} ppm'.format(1e6*ifo.Laser.Wavelength, round(1e6*T[0],3)))
        return(scalarCost, costOut)
    else:
        return(scalarCost)
