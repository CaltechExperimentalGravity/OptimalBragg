import numpy as np
from physunits.length import nm
from physunits.frequency import Hz
from scipy.optimize import differential_evolution as devopt
from scipy.optimize import Bounds
from .layers import multilayer_diel, surfield, calc_abs, field_zmag
from .costs import *
from .noise import brownian_proxy, coating_thermooptic, coating_brownian


def diff_evo(
    stack,
    costs,
    to_pars=None,
    br_pars=None,
    multilayer_diel_pars=None,
    Lmin=20 * nm,
    Lmax=450 * nm,
    verbose=False,
):
    """Wrapped differential evolution optimizer for coatings

    Args:
        stack (dict): Stack attributes (initial coating)
        costs (dict): Target costs and their weights
        to_pars (dict, optional): Kwargs for TO noise calculation
        br_pars (dict, optional): Kwargs for Br noise calculation
        multilayer_diel_pars (dict, optional): Kwargs for multilayer_diel
        Lmin (float, optional): Minimum layer physical thickness [m]
        Lmax (float, optional): Maximum layer physical thickness [m]
        verbose (bool, optional): Optimizer verbosity
    """
    # Unpack stack attributes
    ns, Ls, alphas = stack["ns"], stack["Ls"], stack["alphas"]
    # Global optimization
    vector_mon, conv_mon = [], []

    def diffevo_monitor(xk, convergence):
        """Callback func to track cost: see scipy.differential_evolution
        for more information.

        Args:
            xk (float): Intermediate scalar cost evaluation
            convergence (float): Convergence value

        Returns:
            Bool: Optimizer will halt otherwise
        """
        vector_mon.append(xk)
        conv_mon.append(1 / convergence)
        return False

    opt_res = devopt(
        func=scalar_cost,
        bounds=Bounds(Lmin * np.ones_like(Ls), Lmax * np.ones_like(Ls)),
        mutation=(0.05, 1.5),
        args=(stack, costs, to_pars, br_pars, multilayer_diel_pars, verbose),
        callback=diffevo_monitor,
        updating="deferred",
        strategy="best1bin",
        maxiter=2000,
        popsize=12,
        atol=3e-8,
        tol=1e-3,
        polish=False,
        init="sobol",  # Tryme: "latinhypercube",
        disp=True,
        workers=-1,
    )

    # Optimized stack
    Lres = opt_res.x
    final_vector_cost, vector_weights = vector_cost(
        Ls, stack, costs, to_pars, br_pars, multilayer_diel_pars, True
    )
    final_vector_score = vector_score(
        Lres, stack, costs, to_pars, br_pars, multilayer_diel_pars
    )
    final_scalar_score = np.mean(np.array(list(final_vector_score.values())))
    results = {
        "Ls": Lres,
        "vector_cost": final_vector_cost,
        "vector_score": final_vector_score,
        "scalar_score": final_scalar_score,
        "convergence": np.array(conv_mon),
        "evolution": np.array(vector_mon),
    }
    return results


def vector_cost(
    Ls,
    stack,
    costs,
    Sto_pars=None,
    Sbr_pars=None,
    multilayer_diel_pars=None,
    verbose=False,
):
    """Vector cost function for a coating stack.

    Args:
        Ls (array): Physical thicknesses [m]
        stack (dict): Stack attributes
        costs (dict): Vector cost in the form:
                      {'cost_name':{'target':t, 'weight':w}}
        Sto_pars (None, optional): Kwargs for TO noise calculation
        Sbr_pars (None, optional): Kwargs for Br noise calculation
        multilayer_diel_pars (None, optional): Kwargs for multilayer_diel
        verbose (bool, optional): Verbosity of vector cost

    Returned:
        vector_cost, vector weights (dict, dict): Vector cost and weigths
    """
    if multilayer_diel_pars is None:
        multilayer_diel_pars = {}
    ns = stack["ns"]
    vector_cost = {}
    vector_weights = {}
    for cost, specs in costs.items():
        if specs["weight"]:
            if cost == "T":
                for lam, target in specs["target"].items():
                    if specs["weight"][lam]:
                        vector_cost[cost + f"{int(lam/nm):d}"] = specs[
                            "weight"
                        ][lam] * trans_cost(
                            Ls, target, stack, lam, **multilayer_diel_pars
                        )
                        vector_weights[cost + f"{int(lam/nm):d}"] = specs[
                            "weight"
                        ][lam]
            if cost == "Esurf":
                for lam, target in specs["target"].items():
                    if specs["weight"][lam]:
                        vector_cost[cost + f"{int(lam/nm):d}"] = specs[
                            "weight"
                        ][lam] * surfield_cost(
                            Ls, target, stack, lam, **multilayer_diel_pars
                        )
                        vector_weights[cost + f"{int(lam/nm):d}"] = specs[
                            "weight"
                        ][lam]
            if cost == "Lsens":
                for lam, target in specs["target"].items():
                    if specs["weight"][lam]:
                        vector_cost[cost + f"{int(lam/nm):d}"] = specs[
                            "weight"
                        ][lam] * sens_cost(
                            Ls, target, stack, lam, **multilayer_diel_pars
                        )
                        vector_weights[cost + f"{int(lam/nm):d}"] = specs[
                            "weight"
                        ][lam]
            if cost == "Sbr":
                vector_cost[cost] = specs["weight"] * brownian_cost(
                    Ls, specs["target"], stack, **Sbr_pars
                )
                vector_weights[cost] = specs["weight"]
            if cost == "Sto":
                vector_cost[cost] = specs["weight"] * thermooptic_cost(
                    Ls, specs["target"], stack, **Sto_pars
                )
                vector_weights[cost] = specs["weight"]
            if cost == "abs":
                vector_cost[cost] = specs["weight"] * absorb_cost(
                    Ls, specs["target"], stack
                )
                vector_weights[cost] = specs["weight"]
            if cost == "Lstdev":
                vector_cost[cost] = specs["weight"] * var_cost(
                    Ls, specs["target"]
                )
                vector_weights[cost] = specs["weight"]
    if verbose:
        print(vector_cost)
        for cost in vector_cost.keys():
            print(cost + f" cost = {vector_cost[cost]:.4f}")
    return vector_cost, vector_weights


def scalar_cost(
    Ls,
    stack,
    costs,
    Sto_pars=None,
    Sbr_pars=None,
    multilayer_diel_pars=None,
    verbose=False,
):
    """Scalar cost function for a coating stack.

    Args:
        Ls (arr): Physical thickness array
        stack (dict): Stack attributes
        costs (dict): Vector cost in the form:
                      {'cost_name':{'target':t, 'weight':w}}
        Sto_pars (dict, optional): Kwargs for TO noise calculation
        Sbr_pars (dict, optional): Kwargs for Br noise calculation
        multilayer_diel_pars (dict, optional): Kwargs for multilayer_diel
        verbose (bool, optional): Description

    Returned:
        (float): Weighted scalar cost
    """
    scalar_cost, weight_sum = 0.0, 0.0
    weighted_vector_cost, weight_vector = vector_cost(
        Ls, stack, costs, Sto_pars, Sbr_pars, multilayer_diel_pars, verbose
    )
    for cost, evaluation in weighted_vector_cost.items():
        weight_sum += weight_vector[cost]
        scalar_cost += evaluation
    return scalar_cost / weight_sum


def vector_score(Ls, stack, costs, to_pars, br_pars, multilayer_diel_pars):
    """Score stack relative to an optimization target
    Args:
        Ls (array): Physical thicknesses [m]
        stack (dict): Stack attributes
        costs (dict): Vector cost in the form:
                      {'cost_name':{'target':t, 'weight':w}}
        Sto_pars (None, optional): Kwargs for TO noise calculation
        Sbr_pars (None, optional): Kwargs for Br noise calculation
        multilayer_diel_pars (None, optional): Kwargs for multilayer_diel
        verbose (bool, optional): Verbosity of vector cost

    Returned:
        rel_score (dict): Relative scores for vector cost
    """
    if multilayer_diel_pars is None:
        multilayer_diel_pars = {}
    stack["Ls"] = Ls
    rel_score = {}
    for cost, specs in costs.items():
        if specs["weight"]:
            if cost == "T":
                for lam, target in specs["target"].items():
                    if specs["weight"][lam]:
                        actual = trans(lam, stack, **multilayer_diel_pars)
                        abserr = np.abs(
                            (min(actual, target) - max(actual, target))
                        )
                        rel_score[cost] = np.abs(max(actual, target)) / abserr
            if cost == "Esurf":
                for lam, target in specs["target"].items():
                    if specs["weight"][lam]:
                        rr = refl(lam, stack, **multilayer_diel_pars)
                        actual = surfield(rr, normalized=True)
                        abserr = np.abs(
                            (min(actual, target) - max(actual, target))
                        )
                        rel_score[cost] = np.abs(max(actual, target)) / abserr
            if cost == "Lsens":
                for lam, target in specs["target"].items():
                    if specs["weight"][lam]:
                        Tref = trans(lam, stack, **multilayer_diel_pars)
                        perturb_stack = stack.copy()
                        perturb_stack["Ls"] = 1.01 * Ls
                        Tper = trans(lam, perturb_stack, **multilayer_diel_pars)
                        actual = np.abs((Tper - Tref) / Tref)
                        abserr = np.abs(
                            (min(actual, target) - max(actual, target))
                        )
                        rel_score[cost] = np.abs(max(actual, target)) / abserr
            if cost == "Sbr":
                actual = coating_brownian(
                    np.array([br_pars.pop("freq")]),
                    stack,
                    br_pars["w_beam"],
                    br_pars["power"],
                    br_pars["m_mirror"],
                )[0]
                target = specs["target"]
                abserr = np.abs((min(actual, target) - max(actual, target)))
                rel_score[cost] = np.abs(max(actual, target)) / abserr
            if cost == "Sto":
                actual, _, _ = coating_thermooptic(
                    to_pars.pop("freq"), stack, **to_pars
                )
                target = specs["target"]
                abserr = np.abs((min(actual, target) - max(actual, target)))
                rel_score[cost] = np.abs(max(actual, target)) / abserr
            if cost == "abs":
                _, Enorm = field_zmag(
                    stack["ns"], Ls, n_pts=2**6, lam=stack["lam_ref"]
                )
                actual = calc_abs(Enorm, Ls, stack["alphas"])
                target = specs["target"]
                abserr = np.abs((min(actual, target) - max(actual, target)))
                rel_score[cost] = np.abs(max(actual, target)) / abserr
            if cost == "Lstdev":
                actual = np.var(Ls)
                target = specs["target"]
                abserr = np.abs((min(actual, target) - max(actual, target)))
                rel_score[cost] = np.abs(max(actual, target)) / abserr
    return rel_score
