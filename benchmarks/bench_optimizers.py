#!/usr/bin/env python
"""Benchmark optimizer algorithms and hyperparameters.

Runs the Voyager aSiSiN ETM optimization with different scipy optimizers
and hyperparameter configurations.  Reports wall time, final cost,
number of evaluations, and achieved transmission.

Usage::

    source ~/.zshrc && conda activate coatingDev
    python benchmarks/bench_optimizers.py
    python benchmarks/bench_optimizers.py --reps 5
    python benchmarks/bench_optimizers.py --quick  # 1 rep, fewer configs
"""

import argparse
import json
import signal
import sys
import warnings
from datetime import datetime
from pathlib import Path
from timeit import default_timer

import numpy as np
from scipy.optimize import (
    differential_evolution,
    dual_annealing,
    shgo,
    basinhopping,
)

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from OptimalBragg import qw_stack, load_materials_yaml, Material
from OptimalBragg.costs import getMirrorCost, precompute_misc
from OptimalBragg.io import yamlread
from OptimalBragg.noise import brownian_proxy
from OptimalBragg.optimizer import _build_stack
from OptimalBragg.materials import air


# ── Test case setup ─────────────────────────────────────────────────

PARAMS_PATH = Path(__file__).resolve().parent.parent / \
    "projects" / "Voyager_aSiSiN" / "ETM_params.yml"
TIMEOUT = 120  # seconds per run


class TimeoutError(Exception):
    pass


def _timeout_handler(signum, frame):
    raise TimeoutError("Run exceeded timeout")


def setup_test_case():
    """Load and return the Voyager aSiSiN ETM test case.

    Returns (bounds, args_tuple, stack, costs, misc) where args_tuple
    is (costs, stack, gam, False, misc) — ready for getMirrorCost.
    """
    params_path = Path(PARAMS_PATH)
    params_dir = params_path.parent

    opt_params = yamlread(str(params_path))
    costs = opt_params["costs"]
    misc = opt_params["misc"]

    mat_file = params_dir / misc["materials_file"]
    materials = load_materials_yaml(str(mat_file))

    Npairs = misc["Npairs"]
    hwcap = misc.get("hwcap", "")
    stack = _build_stack(materials, Npairs, optic="ETM", hwcap=hwcap)

    # Attach mirror geometry
    raw_mat = yamlread(str(mat_file))
    sub_overrides = {}
    if "substrate" in raw_mat and "overrides" in raw_mat["substrate"]:
        sub_overrides = raw_mat["substrate"]["overrides"]
    if "MassRadius" in sub_overrides:
        misc.setdefault("r_mirror", sub_overrides["MassRadius"])
    if "MassThickness" in sub_overrides:
        misc.setdefault("d_mirror", sub_overrides["MassThickness"])
    if stack.get("w_beam"):
        misc.setdefault("w_beam", stack["w_beam"])

    gam = brownian_proxy(stack)

    n_layers = len(stack["Ls_opt"])
    wavelength = materials["laser"]["wavelength"]
    min_thick = 10e-9 * stack["ns"][1] / wavelength
    bounds = ((min_thick, 0.48),) + ((0.05, 0.48),) * (n_layers - 1)

    precompute_misc(costs, stack, misc)

    args_tuple = (costs, stack, gam, False, misc)
    return bounds, args_tuple, stack, costs, misc


def jit_warmup(bounds, args_tuple):
    """Run a few evaluations to warm up Numba JIT."""
    x0 = np.array([(b[0] + b[1]) / 2 for b in bounds])
    for _ in range(5):
        getMirrorCost(x0, *args_tuple)


def extract_transmission(L, costs, stack, gam, misc):
    """Run verbose eval to get T1064."""
    final_misc = dict(misc)
    final_misc.update({"Ncopies": 0, "Nfixed": 0})
    precompute_misc(costs, stack, final_misc)
    _, output = getMirrorCost(
        L, costs=costs, stack=stack, gam=gam,
        verbose=True, misc=final_misc,
    )
    return output.get("T1064")


# ── Optimizer runners ───────────────────────────────────────────────

def run_de(bounds, args_tuple, **kwargs):
    """Run differential_evolution with given kwargs, return result dict."""
    defaults = dict(
        strategy="best1bin",
        popsize=15,
        mutation=(0.5, 1.5),
        recombination=0.7,
        tol=1e-2,
        updating="deferred",
        init="halton",
        workers=1,
        maxiter=10000,
        atol=3e-9,
        polish=True,
        disp=False,
    )
    defaults.update(kwargs)

    n_evals = [0]

    def objective(L, *a):
        n_evals[0] += 1
        return getMirrorCost(L, *a)

    res = differential_evolution(
        func=objective,
        bounds=bounds,
        args=args_tuple,
        **defaults,
    )
    return res.x, res.fun, n_evals[0], res.success


def run_dual_annealing(bounds, args_tuple, **kwargs):
    """Run dual_annealing."""
    defaults = dict(
        maxiter=1000,
    )
    defaults.update(kwargs)

    n_evals = [0]

    def objective(L):
        n_evals[0] += 1
        return getMirrorCost(L, *args_tuple)

    lw = [b[0] for b in bounds]
    up = [b[1] for b in bounds]

    res = dual_annealing(
        func=objective,
        bounds=list(zip(lw, up)),
        **defaults,
    )
    return res.x, res.fun, n_evals[0], res.success


def run_basinhopping(bounds, args_tuple, **kwargs):
    """Run basin hopping."""
    defaults = dict(
        niter=200,
    )
    defaults.update(kwargs)

    n_evals = [0]
    x0 = np.array([(b[0] + b[1]) / 2 for b in bounds])

    def objective(L):
        n_evals[0] += 1
        # Penalize out-of-bounds
        for i, (lo, hi) in enumerate(bounds):
            if L[i] < lo or L[i] > hi:
                return 1e10
        return getMirrorCost(L, *args_tuple)

    res = basinhopping(
        func=objective,
        x0=x0,
        **defaults,
    )
    return res.x, res.fun, n_evals[0], True


# ── Configuration grid ──────────────────────────────────────────────

def build_configs(quick=False):
    """Build list of (name, runner_func, kwargs) tuples.

    Sweeps one parameter at a time from the DE baseline.
    """
    configs = []

    # Baseline
    configs.append(("DE:baseline", run_de, {}))

    # Strategy sweep
    for strat in ["best2bin", "rand1bin", "currenttobest1bin"]:
        configs.append((f"DE:strategy={strat}", run_de, {"strategy": strat}))

    if not quick:
        # Popsize sweep
        for ps in [10, 20, 30]:
            configs.append((f"DE:popsize={ps}", run_de, {"popsize": ps}))

        # Mutation sweep
        for mut in [(0.5, 1.0), (0.7, 1.7)]:
            configs.append((
                f"DE:mutation={mut}",
                run_de,
                {"mutation": mut},
            ))

        # Recombination sweep
        for rc in [0.5, 0.9]:
            configs.append((
                f"DE:recombination={rc}",
                run_de,
                {"recombination": rc},
            ))

        # Tolerance sweep
        for tol in [1e-3, 1e-4]:
            configs.append((f"DE:tol={tol}", run_de, {"tol": tol}))

    # Alternative optimizers
    configs.append(("dual_annealing", run_dual_annealing, {}))
    configs.append(("basin_hopping", run_basinhopping, {}))

    return configs


# ── Main benchmark loop ────────────────────────────────────────────

def run_benchmark(reps=3, quick=False):
    """Run all configurations and return results."""
    print("Setting up Voyager aSiSiN ETM test case...")
    bounds, args_tuple, stack, costs, misc = setup_test_case()

    print("JIT warmup...")
    jit_warmup(bounds, args_tuple)

    configs = build_configs(quick=quick)
    print(f"\n{len(configs)} configurations x {reps} reps = "
          f"{len(configs)*reps} runs\n")

    all_results = []

    for name, runner, kwargs in configs:
        print(f"\n{'='*60}")
        print(f"  {name}")
        print(f"{'='*60}")

        runs = []
        for rep in range(reps):
            print(f"  rep {rep+1}/{reps}...", end=" ", flush=True)

            # Set timeout
            old_handler = signal.signal(signal.SIGALRM, _timeout_handler)
            signal.alarm(TIMEOUT)

            try:
                tic = default_timer()
                x_best, cost, n_evals, success = runner(
                    bounds, args_tuple, **kwargs
                )
                dt = default_timer() - tic

                signal.alarm(0)  # cancel alarm

                t1064 = extract_transmission(
                    x_best, costs, stack,
                    brownian_proxy(stack), misc,
                )

                run_info = {
                    "time": round(dt, 2),
                    "cost": round(float(cost), 8),
                    "n_evals": n_evals,
                    "T1064": float(t1064) if t1064 is not None else None,
                    "success": success,
                    "timed_out": False,
                }
                print(f"cost={cost:.6f}  time={dt:.1f}s  "
                      f"evals={n_evals}  T1064={t1064:.2e}"
                      if t1064 else f"cost={cost:.6f}  time={dt:.1f}s")

            except TimeoutError:
                signal.alarm(0)
                run_info = {
                    "time": TIMEOUT,
                    "cost": float("inf"),
                    "n_evals": 0,
                    "T1064": None,
                    "success": False,
                    "timed_out": True,
                }
                print("TIMEOUT")

            except Exception as e:
                signal.alarm(0)
                run_info = {
                    "time": 0,
                    "cost": float("inf"),
                    "n_evals": 0,
                    "T1064": None,
                    "success": False,
                    "timed_out": False,
                    "error": str(e),
                }
                print(f"ERROR: {e}")

            finally:
                signal.signal(signal.SIGALRM, old_handler)

            runs.append(run_info)

        # Summary for this config
        valid_costs = [r["cost"] for r in runs if r["cost"] < float("inf")]
        valid_times = [r["time"] for r in runs if not r["timed_out"]]

        entry = {
            "name": name,
            "kwargs": {k: str(v) for k, v in kwargs.items()},
            "runs": runs,
            "median_cost": round(float(np.median(valid_costs)), 8)
                if valid_costs else None,
            "median_time": round(float(np.median(valid_times)), 2)
                if valid_times else None,
            "std_cost": round(float(np.std(valid_costs)), 8)
                if len(valid_costs) > 1 else 0.0,
            "n_timeouts": sum(1 for r in runs if r.get("timed_out")),
        }
        all_results.append(entry)

    return all_results


def print_summary(results):
    """Print a sorted summary table."""
    print(f"\n\n{'='*80}")
    print("SUMMARY — sorted by median cost, then median time")
    print(f"{'='*80}")

    # Sort: valid results first (by cost then time), timeouts last
    def sort_key(r):
        mc = r["median_cost"] if r["median_cost"] is not None else 1e20
        mt = r["median_time"] if r["median_time"] is not None else 1e20
        return (mc, mt)

    results_sorted = sorted(results, key=sort_key)

    header = (f"{'Config':<35} {'Med Cost':>10} {'Std Cost':>10} "
              f"{'Med Time':>9} {'Med T1064':>11} {'Timeouts':>8}")
    print(header)
    print("-" * len(header))

    for r in results_sorted:
        mc = f"{r['median_cost']:.6f}" if r["median_cost"] is not None else "N/A"
        sc = f"{r['std_cost']:.6f}" if r["std_cost"] else "0"
        mt = f"{r['median_time']:.1f}s" if r["median_time"] is not None else "N/A"

        # Median T1064
        t_vals = [run["T1064"] for run in r["runs"]
                  if run["T1064"] is not None]
        med_t = f"{np.median(t_vals):.2e}" if t_vals else "N/A"
        to = str(r["n_timeouts"])

        print(f"{r['name']:<35} {mc:>10} {sc:>10} {mt:>9} "
              f"{med_t:>11} {to:>8}")

    # Highlight best
    valid = [r for r in results_sorted if r["median_cost"] is not None]
    if valid:
        best = valid[0]
        print(f"\nBest: {best['name']}  "
              f"(cost={best['median_cost']:.6f}, "
              f"time={best['median_time']:.1f}s)")


def save_results(results, outpath):
    """Save results to JSON (handle inf/nan)."""
    def sanitize(obj):
        if isinstance(obj, float):
            if np.isinf(obj) or np.isnan(obj):
                return None
            return obj
        if isinstance(obj, dict):
            return {k: sanitize(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [sanitize(v) for v in obj]
        return obj

    with open(outpath, "w") as f:
        json.dump(sanitize(results), f, indent=2)
    print(f"\nResults saved to {outpath}")


# ── CLI ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Benchmark optimizer algorithms for coating design"
    )
    parser.add_argument("--reps", type=int, default=3,
                        help="Number of repetitions per config (default: 3)")
    parser.add_argument("--quick", action="store_true",
                        help="Quick mode: 1 rep, fewer configs")
    parser.add_argument("--output", type=str,
                        default="benchmarks/results_optimizers.json",
                        help="Output JSON path")
    args = parser.parse_args()

    if args.quick:
        args.reps = 1

    print(f"Optimizer Benchmark — {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print(f"Test case: Voyager aSiSiN ETM ({PARAMS_PATH.relative_to(PARAMS_PATH.parent.parent.parent)})")
    print(f"Reps: {args.reps}, Quick: {args.quick}")

    results = run_benchmark(reps=args.reps, quick=args.quick)
    print_summary(results)
    save_results(results, args.output)
