import numpy as np
from physunits.frequency import Hz
from .layers import multilayer_diel, field_zmag, calc_abs
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
def trans_cost(Ls, target, stack, **multilayer_diel_pars):
    """Evaluates a cost based on stack reflectivity/transmission

    Args:
        Ls (arr): Physical thicknesses
        target (float): Target transmission [0 to 1]
        stack (dict): Stack attributes
        lamb (arr): Wavelength(s) at which to evaluate the
                                reflectivity.
        norm (str, optional): l1, l2, or hyp
        multilayer_diel_pars (dict, optional): Kwargs for evaluating
                                                the transmission.

    Returns:
        (float): Scalar cost for the evaluated transmission

    Deleted Parameters:

    """
    ns, lamb = stack["ns"], stack["lam_ref"]
    rr, _ = multilayer_diel(ns, Ls, lamb, **multilayer_diel_pars)
    T = 1 - np.abs(rr) ** 2
    return np.abs((target - T) / target)


@norm("l1")
def sens_cost(Ls, target, stack, **multilayer_diel_pars):
    """Evaluates the sensitivity cost of a stack transmission
    at given wavelength relative to a 1% Ls perturbation

    Args:
        target(float): Target sensitivity
        ns (arr): Refractive index
        Ls (arr): Physical thickness
        lamb (float): Wavelength at which to evaluate sens
        multilayer_diel_pars (dict, optional): Other parameters for evaluating
                                               the transmission.
    Returns:
        (float): Scalar cost for the evaluated sens
    """
    ns, lamb = stack["ns"], stack["lam_ref"]
    rr0, _ = multilayer_diel(ns, Ls, lamb, **multilayer_diel_pars)
    T0 = 1 - np.abs(rr0) ** 2
    rr1, _ = multilayer_diel(ns, 1.01 * Ls, lamb, **multilayer_diel_pars)
    T1 = 1 - np.abs(rr1) ** 2
    return np.abs(100 * (T0 - T1) / T0)


@norm("arcsinh")
def surfield_cost(Ls, target, stack, **multilayer_diel_pars):
    """Evaluates the normalized surface E field cost

    Args:
           target(float): Target sensitivity
           ns (arr): Refractive index
           Ls (arr): Physical thickness
           lamb (float): Wavelength at which to evaluate sens
           lamb_0 (float): Reference wavelength
           multilayer_diel_pars (dict, optional): Other parameters for evaluating
                                                  the transmission.
       Returns:
           cost (float): Scalar cost for the normalized surf E-field
    """
    ns, lamb = stack["ns"], stack["lam_ref"]
    rr, _ = multilayer_diel(ns, Ls, lamb, **multilayer_diel_pars)
    return np.abs(1 + rr) ** 2


@norm("l2")
def var_cost(Ls, target):
    """Evaluate relative variation of thicknesses

    Args:
        target (float): Target standard deviation
        ns (arr): Refractive index
        Ls (arr): Physical thickness
    Returns:
        (float): Scalar cost for the thickness variation
    """
    rel_std = Ls.std() / Ls.mean()
    return np.abs((target - rel_std) / target)


@norm("l1")
def brownian_cost(Ls, target, stack, **Sbr_pars):
    """Evaluate coating brownian noise according to the
        formula from E0900068 pg4.

    Args:
        target (float): Target Brownian noise proxy
        ns (arr): Refractive index
        Ls (arr): Physical thickness
        gam (dict): Layer keyed prefactors for Brownian noise proxy
    Returns:
        (float): Scalar cost for Brownian noise
    """
    #####################################################################
    ### TODO: How to adapt E0900068 pg 4 for nonbinary coatings? See below:
    # ns = stack["ns"][1:-1]
    # opt_Ls = ns * Ls / stack["lam_ref"]
    # proxy = np.take(opt_Ls, np.argwhere(ns == ns.min()).flatten()).sum()
    #
    # for X, gam in gams.items():
    #     opt_LXs = np.array(
    #         [opt_Ls[j] for j, Xj in enumerate(stack["pattern"]) if Xj == X]
    #     )
    #     proxy += gam * opt_LXs.sum()
    #####################################################################
    Br_stack = stack
    Br_stack["Ls"] = Ls

    Sbr = coating_brownian(
        np.array([100 * Hz]),
        Br_stack,
        Sbr_pars.pop("w_beam"),
        Sbr_pars.pop("power"),
        Sbr_pars.pop("m_mirror"),
    )[0]
    return np.abs((target - Sbr) / target)


@norm("l1")
def thermooptic_cost(Ls, target, stack, **Sto_pars):
    """Evaluates thermo-optic noise cost

    Args:
        target (float): Inverse target PSD of CTO noise [m^2/Hz]
        freq (float): Frequency [Hz] at which to evaluate noise
        ns (arr): Refractive index
        Ls (arr): Physical thicknesses
        mirror (gwinc.Struct): Gwinc sctructure for the mirror
        ifo (gwinc.Struct): Gwinc structure

    Returns:
        (float): Scalar cost for the TO noise.
    """
    TO_stack = stack
    TO_stack["Ls"] = Ls
    Sto, _, _ = coating_thermooptic(100 * Hz, TO_stack, **Sto_pars)
    return np.abs((target - Sto) / target)


@norm("l1")
def absorb_cost(Ls, target, stack):
    """Evaluate cost for integrated absorption

    Args:
        target (float): Target stack absorption
        ns (arr): Refractive index
        Ls (arr): Physical thickness
        ifo (gwinc.Struct): Gwinc structure

    Returns:
        (float): Scalar cost of absorption
    """
    alphas = stack["alphas"]
    z_arr, Enorm = field_zmag(
        stack["ns"], Ls, n_pts=2**2, lam=stack["lam_ref"]
    )
    absorp = calc_abs(Enorm, Ls, alphas)
    return np.abs((target - absorp) / target)
