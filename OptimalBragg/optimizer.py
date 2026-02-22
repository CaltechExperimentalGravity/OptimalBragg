"""Coating optimization via differential evolution.

Extracts the common pattern from project-specific mkETM/mkITM scripts
into a single entry point.  Usage::

    from OptimalBragg.optimizer import run_optimization
    result = run_optimization('projects/aLIGO/ETM_params.yml')

Or from the command line::

    python -m OptimalBragg.optimizer projects/aLIGO/ETM_params.yml
"""

import os
import warnings
from datetime import datetime
from pathlib import Path
from timeit import default_timer

import numpy as np
from scipy.optimize import differential_evolution

from OptimalBragg import Material, qw_stack, load_materials_yaml
from OptimalBragg.io import yamlread, h5write
from OptimalBragg.costs import getMirrorCost, precompute_misc
from OptimalBragg.noise import brownian_proxy


def _build_stack(materials, Npairs, optic=None):
    """Build a QW stack from loaded materials config.

    Parameters
    ----------
    materials : dict
        Output of :func:`load_materials_yaml`.
    Npairs : int
        Number of bilayer pairs.
    optic : str, optional
        Optic name (e.g. ``'ETM'``, ``'ITM'``) to look up beam radius
        in ``materials["optics"]``.

    Returns
    -------
    stack : dict
        Stack dict from :func:`qw_stack`.
    """
    wavelength = materials["laser"]["wavelength"]
    substrate = materials["substrate"]
    thin_films = materials["thin_films"]

    # Superstrate is always vacuum
    from OptimalBragg.materials import air
    superstrate = Material(air)

    stack = qw_stack(
        lam_ref=wavelength,
        substrate=substrate,
        superstrate=superstrate,
        thin_films=thin_films,
        pattern="LH" * Npairs,
    )

    # Attach beam radius and mirror geometry from materials config
    if optic and "optics" in materials:
        optic_cfg = materials["optics"].get(optic, {})
        stack["w_beam"] = optic_cfg.get("beam_radius")

    return stack


def run_optimization(params_path, save=True, optic=None):
    """Run a full coating optimization.

    Parameters
    ----------
    params_path : str or Path
        Path to the cost/optimizer YAML (e.g. ``ETM_params.yml``).
        The ``misc.materials_file`` key in this YAML points to
        the materials config, resolved relative to *params_path*.
    save : bool, optional
        Whether to save HDF5 output.  Default True.
    optic : str, optional
        Optic name for beam radius lookup (``'ETM'`` or ``'ITM'``).
        If None, inferred from the filename (e.g. ``ETM_params.yml``
        → ``'ETM'``).

    Returns
    -------
    dict
        Result dict with keys ``L``, ``scalar_cost``, ``output``,
        ``trajectory``, ``vec_evol``, ``res``, ``stack``, ``gam``.
    """
    params_path = Path(params_path)
    params_dir = params_path.parent

    # Infer optic name from filename if not given
    if optic is None:
        stem = params_path.stem.upper()
        if "ETM" in stem:
            optic = "ETM"
        elif "ITM" in stem:
            optic = "ITM"

    # Load configs
    opt_params = yamlread(str(params_path))
    costs = opt_params["costs"]
    misc = opt_params["misc"]

    # Load materials (relative to params file)
    mat_file = params_dir / misc["materials_file"]
    materials = load_materials_yaml(str(mat_file))

    Npairs = misc["Npairs"]
    stack = _build_stack(materials, Npairs, optic=optic)

    # Attach mirror geometry from materials overrides to misc
    sub_overrides = {}
    raw_mat = yamlread(str(mat_file))
    if "substrate" in raw_mat and "overrides" in raw_mat["substrate"]:
        sub_overrides = raw_mat["substrate"]["overrides"]
    if "MassRadius" in sub_overrides:
        misc.setdefault("r_mirror", sub_overrides["MassRadius"])
    if "MassThickness" in sub_overrides:
        misc.setdefault("d_mirror", sub_overrides["MassThickness"])
    if stack.get("w_beam"):
        misc.setdefault("w_beam", stack["w_beam"])

    # Brownian noise proxy
    gam = brownian_proxy(stack)

    # Set up bounds
    n_layers = 2 * Npairs
    wavelength = materials["laser"]["wavelength"]
    min_thick = 10e-9 * stack["ns"][1] / wavelength  # ~20 nm cap
    bounds = ((min_thick, 0.48),) + ((0.05, 0.48),) * (n_layers - 1)

    print(f"Optimizing {optic or 'coating'}: {Npairs} bilayers, "
          f"{len(bounds)} free layers")

    tic = default_timer()

    # Pre-compute cached values for the hot loop
    precompute_misc(costs, stack, misc)

    # Convergence monitors
    vector_mon, conv_mon = [], []

    def monitor(xk, convergence):
        vector_mon.append(xk.copy())
        conv_mon.append(1 / convergence if convergence else 0)
        return False

    # Global optimization (workers=1: IPC overhead dwarfs per-eval cost)
    res = differential_evolution(
        func=getMirrorCost,
        bounds=bounds,
        updating="deferred",
        strategy="best1bin",
        mutation=(0.05, 1.5),
        popsize=misc.get("Nparticles", 500),
        init=misc.get("init_method", "halton"),
        workers=1,
        maxiter=2000,
        atol=misc.get("atol", 1e-10),
        tol=misc.get("tol", 1e-3),
        args=(costs, stack, gam, False, misc),
        polish=True,
        callback=monitor,
        disp=True,
    )

    if not res.success:
        warnings.warn(f"Optimizer did not converge: {res.message}")

    dt = default_timer() - tic
    print(f"\nOptimization took {dt:.1f} sec")

    # Expand L with Ncopies and Nfixed
    Lres = res.x.copy()
    Ncopies = misc.get("Ncopies", 0)
    Nfixed = misc.get("Nfixed", 0)

    if Ncopies > 0:
        copied = np.tile(Lres[1:2 * Npairs + 1].copy(), Ncopies)
        Lres = np.append(Lres, copied)

    if Nfixed > 0:
        fixed = np.tile(Lres[-2:].copy(), Nfixed)
        Lres = np.append(Lres, fixed)

    # Final cost evaluation (no copies/fixed — they're already in Lres)
    final_misc = dict(misc)
    final_misc.update({"Ncopies": 0, "Nfixed": 0})
    precompute_misc(costs, stack, final_misc)
    scalar_cost, output = getMirrorCost(
        Lres, costs=costs, stack=stack, gam=gam,
        verbose=True, misc=final_misc,
    )

    result = {
        "L": Lres,
        "scalar_cost": scalar_cost,
        "output": output,
        "trajectory": np.array(conv_mon),
        "vec_evol": np.array(vector_mon) if vector_mon else np.array([]),
        "res": res,
        "stack": stack,
        "gam": gam,
    }

    if save:
        _save_hdf5(result, params_path, params_dir)

    return result


def _save_hdf5(result, params_path, params_dir):
    """Save optimization result to HDF5."""
    tnowstr = datetime.now().strftime("%y%m%d_%H%M%S")

    # Infer optic name for directory structure
    stem = params_path.stem.upper()
    if "ETM" in stem:
        optic_dir = "ETM"
    elif "ITM" in stem:
        optic_dir = "ITM"
    else:
        optic_dir = "output"

    spath = params_dir / "Data" / optic_dir
    os.makedirs(spath, exist_ok=True)

    fname = spath / f"{optic_dir}_Layers_{tnowstr}.hdf5"

    output = dict(result["output"])
    vector_cost = output.pop("vectorCost")

    h5_dict = {
        "trajectory": result["trajectory"],
        "vec_evol": result["vec_evol"],
        "diffevo_output": {
            "vectorCost": {k: np.float64(v) for k, v in vector_cost.items()},
            **{k: (np.float64(v) if np.isscalar(v) else v)
               for k, v in output.items()},
        },
        "params_file": str(params_path),
    }
    h5write(str(fname), h5_dict)
    print(f"\nSaved: {fname}")
    return str(fname)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Run coating optimization"
    )
    parser.add_argument("params", help="Path to cost/optimizer YAML")
    parser.add_argument("--optic", help="Optic name (ETM or ITM)")
    parser.add_argument("--no-save", action="store_true",
                        help="Skip HDF5 output")
    args = parser.parse_args()

    run_optimization(args.params, save=not args.no_save, optic=args.optic)
