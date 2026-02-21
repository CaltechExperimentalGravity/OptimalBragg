"""Apples-to-apples benchmark: master vs refactor branch.

Run from repo root (refactor) or SiN_aSi/ (master).
Auto-detects which branch via import path.
Saves results to a JSON sidecar.

Usage:
    python benchmarks/bench_oldvsnew.py          # refactor branch (from repo root)
    python bench_oldvsnew.py                     # master branch (from SiN_aSi/)
"""
import json
import time
import os
import numpy as np
from scipy.optimize import differential_evolution as devo
import gwinc

# --- Auto-detect branch via imports ---
try:
    from generic.optimUtils import getMirrorCost, brownianProxy
    from generic.coatingUtils import multidiel1
    BRANCH = 'refactor'
    GWINC_PATH = 'SiN_aSi/aSiSiN.yaml'
except ImportError:
    from generic_local.optimUtils import getMirrorCost, brownianProxy
    from generic_local.coatingUtils import multidiel1
    BRANCH = 'master'
    GWINC_PATH = 'aSiSiN.yaml'

# ============================================================
# Fixed config — identical on both branches
# ============================================================
NPAIRS = 14
POPSIZE = 40
NLAYERS = 2 * NPAIRS + 1  # 29 layers
N_MULTIDIEL = 10_000
N_COSTFN = 2_000
SEED = 42

# Build a fixed 14-bilayer stack (quarter-wave)
# n = [air, (low,high)*14, low, substrate] = 31 elements for 29 layers
n_low, n_high = 2.17, 3.65
n_sub = 3.39
n_stack = np.array([1.0] + [n_low, n_high] * NPAIRS + [n_low, n_sub])
L_qw = 0.25 * np.ones(NLAYERS)

# Random but reproducible test stack for getMirrorCost
rng = np.random.default_rng(SEED)
L_rand = rng.uniform(0.05, 0.48, size=NLAYERS)

# Cost dict — use the same active costs as the real ETM optimization
# Only Trans1064, Brownian, Thermooptic have nonzero weight
costs = {
    'Trans1064':    {'target': 5e-6,  'weight': 5},
    'Brownian':    {'target': 0.1,   'weight': 2},
    'Thermooptic': {'target': 1e43,  'weight': 2},
    'Lsens':       {'target': 1e-7,  'weight': 0},
    'Esurf':       {'target': 1e-9,  'weight': 0},
    'Absorption':  {'target': 1e-4,  'weight': 0},
    'Trans532':    {'target': 1000e-6, 'weight': 0},
    'TransOPLEV':  {'target': 0.05,  'weight': 0},
    'Lstdev':      {'target': 0.5,   'weight': 0},
}

misc = {
    'fTO': 100,
    'pol': 'te',
    'aoi': 0,
    'Npairs': NPAIRS,
    'Nfixed': 0,
    'Ncopies': 0,
    'Nparticles': POPSIZE,
    'atol': 3e-9,
    'tol': 1e-2,
    'init_method': 'halton',
    'gwincStructFile': GWINC_PATH,
    'lambdaAUX': 0.756,  # only used by refactor branch
}


def main():
    # Load ifo structure for getMirrorCost and full optimization
    # NOTE: master's thermoopticCost mutates ifo.Optics.ETM.Coating.dOpt
    # in-place (pre-existing bug). The refactor fixed this with copy.copy().
    # This benchmark faithfully reproduces each branch's behavior.
    ifo = gwinc.Struct.from_file(GWINC_PATH)
    gam = brownianProxy(ifo)

    results = {'branch': BRANCH, 'config': {
        'Npairs': NPAIRS, 'Nlayers': NLAYERS, 'popsize': POPSIZE,
    }}

    # ============================================================
    # Benchmark 1: multidiel1
    # ============================================================
    print("\n=== multidiel1 (10k calls) ===")

    # Warm up JIT (no-op on master)
    _ = multidiel1(n_stack, L_qw, np.array([1.0]), 0, 'te')

    # Single wavelength
    t0 = time.perf_counter()
    for _ in range(N_MULTIDIEL):
        multidiel1(n_stack, L_qw, np.array([1.0]), 0, 'te')
    dt1 = time.perf_counter() - t0
    us1 = 1e6 * dt1 / N_MULTIDIEL
    print(f"  1 wavelength:  {us1:.1f} us/call  ({N_MULTIDIEL} calls in {dt1:.2f}s)")

    # 3 wavelengths (consolidated call)
    lamb3 = np.array([1.0, 0.756, 0.297])
    t0 = time.perf_counter()
    for _ in range(N_MULTIDIEL):
        multidiel1(n_stack, L_qw, lamb3, 0, 'te')
    dt3 = time.perf_counter() - t0
    us3 = 1e6 * dt3 / N_MULTIDIEL
    print(f"  3 wavelengths: {us3:.1f} us/call  ({N_MULTIDIEL} calls in {dt3:.2f}s)")

    results['multidiel1'] = {
        '1wl_us': round(us1, 1),
        '3wl_us': round(us3, 1),
        'N_calls': N_MULTIDIEL,
    }

    # ============================================================
    # Benchmark 2: getMirrorCost
    # ============================================================
    print("\n=== getMirrorCost (2k calls) ===")

    # Warm up
    _ = getMirrorCost(L_rand.copy(), costs, ifo, gam, verbose=False, misc=misc)

    t0 = time.perf_counter()
    for _ in range(N_COSTFN):
        getMirrorCost(L_rand.copy(), costs, ifo, gam, verbose=False, misc=misc)
    dt_cost = time.perf_counter() - t0
    ms_cost = 1e3 * dt_cost / N_COSTFN
    print(f"  {ms_cost:.3f} ms/call  ({N_COSTFN} calls in {dt_cost:.2f}s)")

    # Get the scalar value for verification
    cost_val = getMirrorCost(L_rand.copy(), costs, ifo, gam, verbose=False, misc=misc)
    print(f"  scalar cost = {cost_val:.10f}")

    results['getMirrorCost'] = {
        'ms_per_call': round(ms_cost, 3),
        'N_calls': N_COSTFN,
        'scalar_cost': float(cost_val),
    }

    # ============================================================
    # Benchmark 3: Full differential_evolution optimization
    # ============================================================
    print(f"\n=== Full devo optimization ===")
    print(f"  Npairs={NPAIRS}, popsize={POPSIZE}, workers=-1 (all cores)")

    bounds = ((10e-9 / ifo.Laser.Wavelength, 0.48),) + ((0.05, 0.48),) * (NLAYERS - 1)

    t0 = time.perf_counter()
    res = devo(
        func=getMirrorCost,
        bounds=bounds,
        updating='deferred',
        strategy='best1bin',
        mutation=(0.05, 1.5),
        popsize=POPSIZE,
        init='halton',
        workers=-1,
        maxiter=2000,
        atol=3e-9,
        tol=1e-2,
        args=(costs, ifo, gam, False, misc),
        polish=True,
        disp=True,
    )
    dt_opt = time.perf_counter() - t0
    print(f"  Completed in {dt_opt:.1f} sec")
    print(f"  Final cost: {res.fun:.10f}")
    print(f"  Converged: {res.success}")
    print(f"  Iterations: {res.nit}")
    print(f"  Function evals: {res.nfev}")

    results['full_optimization'] = {
        'wall_time_sec': round(dt_opt, 1),
        'final_cost': float(res.fun),
        'converged': bool(res.success),
        'iterations': int(res.nit),
        'nfev': int(res.nfev),
    }

    # ============================================================
    # Save results
    # ============================================================
    outfile = f'results_{BRANCH}.json'
    # If running from SiN_aSi/, save there; otherwise save in benchmarks/
    if os.path.isdir('benchmarks'):
        outfile = os.path.join('benchmarks', outfile)

    with open(outfile, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {outfile}")

    # ============================================================
    # Summary table
    # ============================================================
    print(f"""
{'='*60}
BENCHMARK SUMMARY — {BRANCH}
{'='*60}
multidiel1 (14-bilayer, 10k calls):
  1 wavelength:    {results['multidiel1']['1wl_us']:>8.1f} us/call
  3 wavelengths:   {results['multidiel1']['3wl_us']:>8.1f} us/call

getMirrorCost (14 pairs, 2k calls):
  per call:        {results['getMirrorCost']['ms_per_call']:>8.3f} ms/call
  scalar cost:     {results['getMirrorCost']['scalar_cost']:.10f}

Full optimization (14 pairs, popsize={POPSIZE}, workers=-1):
  wall time:       {results['full_optimization']['wall_time_sec']:>8.1f} sec
  final cost:      {results['full_optimization']['final_cost']:.10f}
  iterations:      {results['full_optimization']['iterations']}
  func evals:      {results['full_optimization']['nfev']}
{'='*60}
""")


if __name__ == '__main__':
    print(f"Detected branch: {BRANCH}")
    main()
