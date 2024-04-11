import numpy as np
from physunits.frequency import Hz
from .layers import *
from .noise import coating_brownian, coating_thermooptic


def norm(norm_arg):
    def costdecorator(costeval):
        def wrapped(*args, **kwargs):
            ci = costeval(*args, **kwargs)
            match norm_arg:
                case "l1":
                    return ci
                case "l2":
                    return ci**2
                case "arcsinh":
                    return np.arcsinh(ci)

        return wrapped

    return costdecorator


@norm("l2")
def trans_cost(Ls, target, stack, lamb, **multilayer_diel_pars):
    """Evaluates a cost based on stack reflectivity/transmission

    Args:
        Ls (arr): Physical thicknesses in m
        target (float): Target transmission [0 to 1]
        stack (dict): Stack attributes
        lamb (arr): Wavelength(s) at which to evaluate the
                                reflectivity.
        **multilayer_diel_pars: Kwargs for multilayer_diel

    Returns:
        (float): Scalar cost for the evaluated transmission

    Deleted Parameters:
        norm (str, optional): l1, l2, or hyp
        multilayer_diel_pars (dict, optional): Kwargs for multilayer_diel

    """
    stack["Ls"] = Ls
    Tact = trans(lamb, stack, **multilayer_diel_pars)
    return np.abs((target - Tact) / target)


@norm("l1")
def sens_cost(Ls, target, stack, lamb, **multilayer_diel_pars):
    """Evaluates the sensitivity cost of a stack transmission
    at given wavelength relative to a +1% Ls perturbation

    Args:
        Ls (arr): Physical thickness [m]
        target (float): Target spectral sensitivity [ppm/um]
        stack (dict): Container of stack attributes
        lamb (float): Wavelength at which to evaluate sensitivity [m]
        **multilayer_diel_pars: Kwargs for multilayer_diel

    Returns:
        (float): Scalar cost for the evaluated sens
    """
    stack["Ls"] = 1.00 * Ls
    Tref = trans(lamb, stack, **multilayer_diel_pars)
    stack["Ls"] = 1.01 * Ls
    Tper = trans(lamb, stack, **multilayer_diel_pars)
    return np.abs(target - np.abs((Tref - Tper) / Tref)) / target


@norm("arcsinh")
def surfield_cost(Ls, target, stack, lamb, **multilayer_diel_pars):
    """Evaluates the normalized surface E field cost for the wavelength lamb

    Args:
        Ls (arr): Physical thicknesses [m]
        target (float): Target surface E-field (normalized)
        stack (dict): Container of stack attributes
        lamb (float): Wavelength in m
        **multilayer_diel_pars: Kwargs for multilayer_diel

    Returns:
        (float): Scalar cost for the normalized surf E-field
    """
    stack["Ls"] = Ls
    rr = refl(lamb, stack, **multilayer_diel_pars)
    return np.abs((target - np.abs(1 + rr) ** 2) / target)


@norm("l1")
def var_cost(Ls, target):
    """Evaluate physical thicknesses variance

    Args:
        Ls (arr): Physical thickness [m]
        target (float): Target variance [m**2]

    Returns:
        (float): Scalar cost for the thickness variance
    """
    return np.abs((target - Ls.var()) / target)


@norm("l1")
def brownian_cost(Ls, target, stack, **Sbr_pars):
    """Evaluate coating brownian noise according to the
        formula from E0900068 pg4.

    Args:
        Ls (arr): Physical thickness [m]
        target (float): Target equivalent displacement Brownian noise PSD [m**2/Hz]
        stack (dict): Container of stack attributes
        **Sbr_pars: Kwargs for coating_brownian noise calculation

    Returns:
        (float): Scalar cost for Brownian noise
    """
    ################################################################
    ### TODO: How to adapt E0900068 pg 4 for nonbinary coatings?
    ### See below for a naive extension:
    # ns = stack["ns"][1:-1]
    # opt_Ls = ns * Ls / stack["lam_ref"]
    # proxy = np.take(opt_Ls, np.argwhere(ns == ns.min()).flatten()).sum()
    #
    # for X, gam in gams.items():
    #     opt_LXs = np.array(
    #         [opt_Ls[j] for j, Xj in enumerate(stack["pattern"]) if Xj == X]
    #     )
    #     proxy += gam * opt_LXs.sum()
    ################################################################
    stack["Ls"] = Ls
    ff = Sbr_pars.pop("freq")
    Sbr = coating_brownian(
        np.array([ff]),
        stack,
        Sbr_pars.pop("w_beam"),
        Sbr_pars.pop("power"),
        Sbr_pars.pop("m_mirror"),
    )[0]
    return np.abs((target - Sbr) / target)


@norm("l1")
def thermooptic_cost(Ls, target, stack, **Sto_pars):
    """Evaluates thermo-optic noise cost

    Args:
        Ls (arr): Physical thicknesses
        target (float): Inverse target PSD of CTO noise [m^2/Hz]
        stack (dict): Stack attributes
        **Sto_pars: Kwargs for TO calculation

    Returns:
        (float): Scalar cost for the TO noise.
    """
    stack["Ls"] = Ls
    ff = Sto_pars.pop("freq")
    Sto, _, _ = coating_thermooptic(ff, stack, **Sto_pars)
    return np.abs((target - Sto) / target)


@norm("l1")
def absorb_cost(Ls, target, stack):
    """Evaluate cost for integrated absorption

    Args:
        Ls (arr): Physical thickness [m]
        target (float): Target stack absorption
        stack (dict): Stack attributes

    Returns:
        (float): Scalar cost of absorption

    """
    # To speed this up, maybe take the peak |E(z)| at each layer,
    # and reduce the integral to a cheaper dot product ?
    # Also consider using the analytic form in:
    # https://arxiv.org/pdf/1603.02720.pdf

    alphas = stack["alphas"]
    z_arr, Enorm = field_zmag(
        stack["ns"], Ls, n_pts=2**2, lam=stack["lam_ref"]
    )
    absorp = calc_abs(Enorm, Ls, alphas)
    return np.abs((target - absorp) / target)
