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
from scipy.optimize import differential_evolution, dual_annealing

from OptimalBragg import Material, qw_stack, load_materials_yaml
from OptimalBragg.io import yamlread, h5write
from OptimalBragg.costs import getMirrorCost, precompute_misc
from OptimalBragg.noise import brownian_proxy


def _build_stack(materials, Npairs, optic=None, hwcap=""):
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
    hwcap : str, optional
        Half-wave cap layer(s) at the superstrate (surface) side.
        ``'L'`` adds a half-wave SiO2 cap for environmental stability
        and low surface E-field.  Default ``""`` (no cap).

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
        hwcap=hwcap,
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
    hwcap = misc.get("hwcap", "")
    stack = _build_stack(materials, Npairs, optic=optic, hwcap=hwcap)

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

    # Set up bounds — includes cap layers + bilayer stack
    n_layers = len(stack["Ls_opt"])  # 2*Npairs + len(hwcap)
    wavelength = materials["laser"]["wavelength"]
    min_thick = 10e-9 * stack["ns"][1] / wavelength  # ~20 nm cap
    bounds = ((min_thick, 0.48),) + ((0.05, 0.48),) * (n_layers - 1)

    cap_str = f" + {hwcap} cap" if hwcap else ""
    print(f"Optimizing {optic or 'coating'}: {Npairs} bilayers{cap_str}, "
          f"{n_layers} free layers")

    tic = default_timer()

    # Pre-compute cached values for the hot loop
    precompute_misc(costs, stack, misc)

    # Convergence monitors
    vector_mon, conv_mon = [], []

    algorithm = misc.get("algorithm", "dual_annealing")

    if algorithm == "dual_annealing":
        # dual_annealing: ~3x faster, lower variance, same cost quality
        def da_callback(x, f, context):
            vector_mon.append(x.copy())
            conv_mon.append(f)

        res = dual_annealing(
            func=getMirrorCost,
            bounds=list(bounds),
            args=(costs, stack, gam, False, misc),
            maxiter=misc.get("maxiter", 1000),
            callback=da_callback,
        )

        if not res.success:
            warnings.warn(f"Optimizer did not converge: {res.message}")

    else:
        # Default: differential_evolution
        def monitor(xk, convergence):
            vector_mon.append(xk.copy())
            conv_mon.append(1 / convergence if convergence else 0)
            return False

        # Convert Nparticles (total population) to popsize (multiplier).
        # scipy DE: total population = popsize * n_vars.
        n_particles = misc.get("Nparticles", 15 * len(bounds))
        popsize = max(15, n_particles // len(bounds))

        res = differential_evolution(
            func=getMirrorCost,
            bounds=bounds,
            updating="deferred",
            strategy=misc.get("strategy", "best1bin"),
            mutation=misc.get("mutation", (0.5, 1.5)),
            recombination=misc.get("recombination", 0.7),
            popsize=popsize,
            init=misc.get("init_method", "halton"),
            workers=1,
            maxiter=misc.get("maxiter", 10000),
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
        "algorithm": algorithm,
        "wall_time": dt,
    }
    if misc.get("_original_params"):
        result["_original_params"] = Path(misc["_original_params"])

    if save:
        hdf5_path = _save_hdf5(result, params_path, params_dir)
        result["hdf5_path"] = hdf5_path

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
        optic_dir = params_path.stem.split('_')[0]

    spath = params_dir / "Data" / optic_dir
    os.makedirs(spath, exist_ok=True)

    nlayers = len(result["L"])
    npairs = (nlayers - 1) // 2 if nlayers % 2 == 1 else nlayers // 2
    fname = spath / f"{optic_dir}_N{npairs}_Layers_{tnowstr}.hdf5"

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
        "params_file": str(result.get("_original_params", params_path)),
        "algorithm": result.get("algorithm", "differential_evolution"),
        "wall_time": np.float64(result.get("wall_time", 0.0)),
    }
    h5write(str(fname), h5_dict)
    print(f"\nSaved: {fname}")
    return str(fname)


def _qw_min_npairs(n_H, n_L, n_sub, n_sup, T_target):
    """Minimum bilayer count for a QW stack to reach *T_target*.

    For (LH)^N on substrate, the high-reflectivity approximation is:
    ``T ≈ 4 n_sup / (n_sub * (n_H / n_L)^{2N})``

    Returns the smallest integer N such that T_QW(N) <= T_target.
    """
    # T = 4 * n_sup / (n_sub * (n_H/n_L)^(2N))
    # (n_H/n_L)^(2N) = 4 * n_sup / (T * n_sub)
    # N = ln(4 * n_sup / (T * n_sub)) / (2 * ln(n_H / n_L))
    rhs = 4 * n_sup / (T_target * n_sub)
    if rhs <= 1.0:
        return 1  # even 1 bilayer suffices
    N = np.log(rhs) / (2 * np.log(n_H / n_L))
    return int(np.ceil(N))


def _run_one_npairs(params_path, N, save, optic):
    """Run a single Npairs optimization.  Designed for process-pool dispatch."""
    import yaml

    params_path = Path(params_path)
    params_dir = params_path.parent

    opt_params = yamlread(str(params_path))
    opt_params["misc"]["Npairs"] = N
    opt_params["misc"]["_original_params"] = str(params_path)

    tmp_yaml = params_dir / f"_sweep_tmp_{optic}_N{N}.yml"
    with open(tmp_yaml, "w") as f:
        yaml.dump(opt_params, f, default_flow_style=False)

    try:
        result = run_optimization(tmp_yaml, save=save, optic=optic)
        entry = {
            "Npairs": N,
            "cost": result["scalar_cost"],
            "T1": result["output"].get("T1"),
            "T2": result["output"].get("T2"),
            "T3": result["output"].get("T3"),
            "result": result,
        }
    except Exception as e:
        print(f"  Npairs={N} FAILED: {e}")
        entry = {"Npairs": N, "cost": np.inf,
                 "T1": None, "T2": None, "T3": None, "result": None}
    finally:
        if tmp_yaml.exists():
            tmp_yaml.unlink()

    return entry


def sweep_nlayers(params_path, n_range=None, save=True, optic=None):
    """Sweep Npairs to find the minimum that hits all targets.

    Runs all Npairs values in parallel using a process pool.

    Parameters
    ----------
    params_path : str or Path
        Path to cost/optimizer YAML.
    n_range : tuple of (int, int), optional
        (min_pairs, max_pairs) to sweep.  If None, uses QW theory
        to estimate a physics-based minimum and sweeps up to +6.
    save : bool, optional
        Whether to save HDF5 for each run.  Default True.
    optic : str, optional
        Optic name (``'ETM'`` or ``'ITM'``).

    Returns
    -------
    list of dict
        One entry per Npairs with keys ``Npairs``, ``cost``,
        ``T1``, ``T2``, ``result``.
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed

    params_path = Path(params_path)
    params_dir = params_path.parent

    if optic is None:
        stem = params_path.stem.upper()
        if "ETM" in stem:
            optic = "ETM"
        elif "ITM" in stem:
            optic = "ITM"

    # Load configs to get material indices and T targets
    opt_params = yamlread(str(params_path))
    costs = opt_params["costs"]
    misc = opt_params["misc"]

    mat_file = params_dir / misc["materials_file"]
    materials = load_materials_yaml(str(mat_file))

    n_H = materials["thin_films"]["H"].Index
    n_L = materials["thin_films"]["L"].Index
    n_sub = materials["substrate"].Index
    n_sup = 1.0  # vacuum

    if n_range is None:
        # Find the tightest T target to set minimum Npairs
        t_targets = []
        if costs.get("Trans1", {}).get("weight", 0):
            t_targets.append(costs["Trans1"]["target"])
        if costs.get("Trans2", {}).get("weight", 0):
            # AUX target is at a different wavelength — QW theory
            # only applies at the design wavelength, so use PSL target
            pass
        if not t_targets:
            raise ValueError("No transmission targets with nonzero weight")

        t_min = min(t_targets)
        n_start = _qw_min_npairs(n_H, n_L, n_sub, n_sup, t_min)
        n_range = (n_start, n_start + 6)

    n_values = list(range(n_range[0], n_range[1] + 1))

    print(f"\n{'='*60}")
    print(f"Nlayers sweep: Npairs = {n_range[0]}..{n_range[1]}  "
          f"({len(n_values)} jobs in parallel)")
    print(f"n_H={n_H:.3f}, n_L={n_L:.3f}, ratio={n_L/n_H:.4f}")
    print(f"{'='*60}\n")

    # Launch all Npairs optimizations in parallel
    results_map = {}
    with ProcessPoolExecutor(max_workers=len(n_values)) as pool:
        futures = {
            pool.submit(_run_one_npairs, str(params_path), N, save, optic): N
            for N in n_values
        }
        for future in as_completed(futures):
            N = futures[future]
            entry = future.result()
            results_map[N] = entry
            if entry["T1"] is not None:
                r2_str = (f"  R2={1e6*(1-entry['T2']):.0f} ppm"
                          if entry["T2"] is not None else "")
                r3_str = (f"  R3={1e6*(1-entry['T3']):.0f} ppm"
                          if entry.get("T3") is not None else "")
                print(f"  Npairs={N:>2} done: cost={entry['cost']:.6f}  "
                      f"T1={entry['T1']:.6e}{r2_str}{r3_str}")
            else:
                print(f"  Npairs={N:>2} done: FAILED")

    # Collect results in Npairs order
    results = [results_map[N] for N in n_values]

    # Print summary table
    has_t3 = any(r.get("T3") is not None for r in results)
    t3_tgt = costs.get("Trans3", {}).get("target")

    header = (f"{'Npairs':>6} | {'Layers':>6} | {'Cost':>10} | "
              f"{'T1':>12} | {'R2 (ppm)':>12}")
    sep = f"{'-'*6}-+-{'-'*6}-+-{'-'*10}-+-{'-'*12}-+-{'-'*12}"
    if has_t3:
        header += f" | {'R3 (ppm)':>12}"
        sep += f"-+-{'-'*12}"

    print(f"\n{'='*60}")
    print(header)
    print(sep)

    t1_tgt = costs.get("Trans1", {}).get("target")
    t2_tgt = costs.get("Trans2", {}).get("target")

    best_n, best_cost = None, np.inf
    for r in results:
        t1_str = f"{r['T1']:.4e}" if r["T1"] is not None else "N/A"
        r2_str = (f"{1e6*(1-r['T2']):.0f}" if r["T2"] is not None
                  else "N/A")
        if r["cost"] < best_cost:
            best_cost = r["cost"]
            best_n = r["Npairs"]
        line = (f"{r['Npairs']:>6} | {2*r['Npairs']:>6} | "
                f"{r['cost']:>10.6f} | {t1_str:>12} | {r2_str:>12}")
        if has_t3:
            r3_str = (f"{1e6*(1-r['T3']):.0f}" if r.get("T3") is not None
                      else "N/A")
            line += f" | {r3_str:>12}"
        print(line)

    r2_tgt_ppm = f"{1e6*(1-t2_tgt):.0f}" if t2_tgt else "N/A"
    targets_str = f"\nTargets:  T1={t1_tgt}  R2<{r2_tgt_ppm} ppm"
    if has_t3:
        r3_tgt_ppm = f"{1e6*(1-t3_tgt):.0f}" if t3_tgt else "N/A"
        targets_str += f"  R3<{r3_tgt_ppm} ppm"
    print(targets_str)
    if best_n:
        print(f"Best:     Npairs={best_n} (cost={best_cost:.6f})")

    return results


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
