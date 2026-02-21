"""Systematic profiling of getMirrorCost and the full optimization pipeline.

Three levels of analysis:
  Level 1: Manual ns-resolution instrumentation of getMirrorCost sections
  Level 2: cProfile call tree for getMirrorCost
  Level 3: devo pipeline comparison (workers=1 vs workers=-1)

Usage:
    source ~/.zshrc && conda activate coatingDev
    python benchmarks/profile_getMirrorCost.py
"""
import cProfile
import io
import pstats
import time
import numpy as np
from scipy.optimize import differential_evolution as devo
import gwinc

from generic.optimUtils import getMirrorCost, brownianProxy, brownianCost, thermoopticCost
from generic.coatingUtils import multidiel1
from generic.thermoopticUtils import coating_thermooptic_fast, extract_ifo_params

# ============================================================
# Shared config (matches bench_oldvsnew.py)
# ============================================================
NPAIRS = 14
NLAYERS = 2 * NPAIRS + 1  # 29
GWINC_PATH = 'SiN_aSi/aSiSiN.yaml'

n_low, n_high = 2.17, 3.65
n_sub = 3.39

costs = {
    'Trans1064':    {'target': 5e-6,    'weight': 5},
    'Brownian':    {'target': 0.1,     'weight': 2},
    'Thermooptic': {'target': 1e43,    'weight': 2},
    'Lsens':       {'target': 1e-7,    'weight': 0},
    'Esurf':       {'target': 1e-9,    'weight': 0},
    'Absorption':  {'target': 1e-4,    'weight': 0},
    'Trans532':    {'target': 1000e-6, 'weight': 0},
    'TransOPLEV':  {'target': 0.05,    'weight': 0},
    'Lstdev':      {'target': 0.5,     'weight': 0},
}

misc = {
    'fTO': 100, 'pol': 'te', 'aoi': 0,
    'Npairs': NPAIRS, 'Nfixed': 0, 'Ncopies': 0,
    'lambdaAUX': 0.756,
}

SEED = 42


def level1_instrumented(ifo, gam, ifo_params, N=10_000):
    """Manual ns-resolution profiling of getMirrorCost internals."""
    print("\n" + "=" * 70)
    print("LEVEL 1: Instrumented getMirrorCost breakdown")
    print(f"  {N} iterations, perf_counter_ns resolution")
    print("=" * 70)

    rng = np.random.default_rng(SEED)
    L_rand = rng.uniform(0.05, 0.48, size=NLAYERS)

    # Pre-compute the n array once (to time it separately)
    aoi, pol = misc['aoi'], misc['pol']
    fTO = misc['fTO']

    # Storage for per-section timings
    sections = [
        'n_array_build', 'active_set_build', 'multidiel1_call',
        'Trans1064_math', 'brownian_cost', 'thermooptic_jit',
        'weighted_sum', 'TOTAL_getMirrorCost',
    ]
    timings = {s: np.empty(N, dtype=np.int64) for s in sections}

    # Warm up JIT
    n_pre = np.array([1.0] + [n_low, n_high] * NPAIRS + [n_low, n_sub])
    _ = multidiel1(n_pre, L_rand, np.array([1.0]), aoi, pol)
    _ = coating_thermooptic_fast(fTO, L_rand, ifo.Laser.Wavelength,
                                  ifo.Optics.ETM.BeamRadius, ifo_params)

    # Also time the unmodified getMirrorCost for cross-check
    misc_with_params = {**misc, '_ifo_params': ifo_params}
    t_total_start = time.perf_counter_ns()
    for _ in range(N):
        getMirrorCost(L_rand.copy(), costs, ifo, gam, verbose=False,
                      misc=misc_with_params)
    t_total_end = time.perf_counter_ns()
    ref_us = (t_total_end - t_total_start) / N / 1000

    # Instrumented loop — reproduces getMirrorCost logic with timing probes
    for i in range(N):
        L = L_rand.copy()

        t0 = time.perf_counter_ns()

        # --- n array construction (lines 271-281) ---
        t1 = time.perf_counter_ns()
        doublet = np.tile(np.array([n_low, n_high]), NPAIRS)
        if len(doublet) != len(L):
            doublet = np.append(doublet, doublet[0])
        n = np.append(1, doublet)
        n = np.append(n, n_sub)
        t2 = time.perf_counter_ns()
        timings['n_array_build'][i] = t2 - t1

        # --- Active set + wavelength list (lines 287-304) ---
        t1 = time.perf_counter_ns()
        active = {c for c, s in costs.items() if s['weight']}
        wl_list = []
        wl_map = {}
        if 'Trans1064' in active or 'Esurf' in active:
            wl_map['PSL'] = len(wl_list)
            wl_list.append(1.0)
        t2 = time.perf_counter_ns()
        timings['active_set_build'][i] = t2 - t1

        # --- multidiel1 call (lines 306-308) ---
        t1 = time.perf_counter_ns()
        r_main, _ = multidiel1(n, L, np.array(wl_list), aoi, pol)
        T_main = 1 - np.abs(r_main)**2
        t2 = time.perf_counter_ns()
        timings['multidiel1_call'][i] = t2 - t1

        # --- Trans1064 math (lines 311-315) ---
        t1 = time.perf_counter_ns()
        idx = wl_map['PSL']
        T1064 = T_main[idx]
        cost_Trans1064 = np.abs((costs['Trans1064']['target'] - T1064)
                               / costs['Trans1064']['target'])**2
        t2 = time.perf_counter_ns()
        timings['Trans1064_math'][i] = t2 - t1

        # --- Brownian cost (line 344) ---
        t1 = time.perf_counter_ns()
        cost_Brownian = brownianCost(costs['Brownian']['target'], L, gam)
        t2 = time.perf_counter_ns()
        timings['brownian_cost'][i] = t2 - t1

        # --- Thermooptic JIT (lines 347-353) ---
        t1 = time.perf_counter_ns()
        cost_TO = thermoopticCost(costs['Thermooptic']['target'], fTO, L, ifo,
                                   ifo_params=ifo_params)
        t2 = time.perf_counter_ns()
        timings['thermooptic_jit'][i] = t2 - t1

        # --- Weighted sum (lines 373-374) ---
        t1 = time.perf_counter_ns()
        scalar_cost = (costs['Trans1064']['weight'] * cost_Trans1064
                       + costs['Brownian']['weight'] * cost_Brownian
                       + costs['Thermooptic']['weight'] * cost_TO)
        t2 = time.perf_counter_ns()
        timings['weighted_sum'][i] = t2 - t1

        timings['TOTAL_getMirrorCost'][i] = time.perf_counter_ns() - t0

    # Print results
    print(f"\n  Reference (unmodified getMirrorCost): {ref_us:.1f} us/call\n")
    print(f"  {'Section':<22s}  {'Mean (us)':>10s}  {'Median':>8s}  "
          f"{'Std':>8s}  {'% total':>8s}")
    print("  " + "-" * 62)

    total_mean = np.mean(timings['TOTAL_getMirrorCost']) / 1000
    for s in sections:
        arr = timings[s] / 1000  # ns -> us
        mean = np.mean(arr)
        med = np.median(arr)
        std = np.std(arr)
        pct = 100 * mean / total_mean if s != 'TOTAL_getMirrorCost' else 100.0
        print(f"  {s:<22s}  {mean:>10.2f}  {med:>8.2f}  {std:>8.2f}  {pct:>7.1f}%")

    # Compute unaccounted overhead
    accounted = sum(np.mean(timings[s]) for s in sections
                    if s != 'TOTAL_getMirrorCost')
    unaccounted_us = (np.mean(timings['TOTAL_getMirrorCost']) - accounted) / 1000
    pct_un = 100 * unaccounted_us / total_mean
    print(f"  {'(unaccounted)':<22s}  {unaccounted_us:>10.2f}  {'':>8s}  "
          f"{'':>8s}  {pct_un:>7.1f}%")

    return timings


def level2_cprofile(ifo, gam, ifo_params, N=10_000):
    """cProfile call tree for getMirrorCost."""
    print("\n" + "=" * 70)
    print("LEVEL 2: cProfile call tree")
    print(f"  {N} calls to getMirrorCost")
    print("=" * 70)

    rng = np.random.default_rng(SEED)
    L_rand = rng.uniform(0.05, 0.48, size=NLAYERS)
    misc_with_params = {**misc, '_ifo_params': ifo_params}

    # Warm up
    _ = getMirrorCost(L_rand.copy(), costs, ifo, gam, verbose=False,
                      misc=misc_with_params)

    pr = cProfile.Profile()
    pr.enable()
    for _ in range(N):
        getMirrorCost(L_rand.copy(), costs, ifo, gam, verbose=False,
                      misc=misc_with_params)
    pr.disable()

    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
    ps.print_stats(25)
    print(s.getvalue())

    # Also print by tottime
    s2 = io.StringIO()
    ps2 = pstats.Stats(pr, stream=s2).sort_stats('tottime')
    ps2.print_stats(25)
    print("\n  --- Sorted by total time ---")
    print(s2.getvalue())


def level3_devo_overhead(ifo, gam, ifo_params, popsize=200):
    """Compare devo with workers=1 vs workers=-1 to isolate MP overhead."""
    print("\n" + "=" * 70)
    print("LEVEL 3: differential_evolution pipeline overhead")
    print(f"  Npairs={NPAIRS}, popsize={popsize}, maxiter=5")
    print("=" * 70)

    bounds = ((10e-9 / ifo.Laser.Wavelength, 0.48),) + ((0.05, 0.48),) * (NLAYERS - 1)
    misc_with_params = {**misc, '_ifo_params': ifo_params}

    common_kw = dict(
        func=getMirrorCost,
        bounds=bounds,
        updating='deferred',
        strategy='best1bin',
        mutation=(0.05, 1.5),
        popsize=popsize,
        init='halton',
        maxiter=5,
        atol=0,
        tol=0,
        args=(costs, ifo, gam, False, misc_with_params),
        seed=SEED,
        polish=False,
    )

    # --- workers=1 (serial, no multiprocessing) ---
    # Warm up JIT
    rng_tmp = np.random.default_rng(0)
    L_tmp = rng_tmp.uniform(0.05, 0.48, size=NLAYERS)
    _ = getMirrorCost(L_tmp, costs, ifo, gam, verbose=False,
                      misc=misc_with_params)

    t0 = time.perf_counter()
    res1 = devo(**common_kw, workers=1)
    dt1 = time.perf_counter() - t0
    us_per_eval_1 = 1e6 * dt1 / res1.nfev

    print(f"\n  workers=1 (serial):")
    print(f"    Wall time:     {dt1:.3f} sec")
    print(f"    Func evals:    {res1.nfev}")
    print(f"    Per eval:      {us_per_eval_1:.1f} us/eval (effective)")

    # --- workers=-1 (all cores) ---
    t0 = time.perf_counter()
    res_mp = devo(**common_kw, workers=-1)
    dt_mp = time.perf_counter() - t0
    us_per_eval_mp = 1e6 * dt_mp / res_mp.nfev

    import os
    ncpu = os.cpu_count()
    print(f"\n  workers=-1 ({ncpu} cores):")
    print(f"    Wall time:     {dt_mp:.3f} sec")
    print(f"    Func evals:    {res_mp.nfev}")
    print(f"    Per eval:      {us_per_eval_mp:.1f} us/eval (effective)")

    # --- Analysis ---
    ideal_mp = dt1 / ncpu
    efficiency = ideal_mp / dt_mp * 100
    overhead_per_eval = us_per_eval_mp - us_per_eval_1 / ncpu

    print(f"\n  Analysis:")
    print(f"    Ideal MP time:       {ideal_mp:.3f} sec ({ncpu} cores, zero overhead)")
    print(f"    Actual MP time:      {dt_mp:.3f} sec")
    print(f"    Parallelization eff: {efficiency:.0f}%")
    print(f"    Speedup:             {dt1/dt_mp:.2f}x (ideal: {ncpu}x)")
    print(f"    Overhead per eval:   {overhead_per_eval:.1f} us "
          f"(serialization + IPC)")

    # --- cProfile of devo with workers=1 to see scipy overhead ---
    print(f"\n  cProfile of devo (workers=1, maxiter=3):")
    common_kw_short = {**common_kw, 'maxiter': 3}
    pr = cProfile.Profile()
    pr.enable()
    _ = devo(**common_kw_short, workers=1)
    pr.disable()

    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
    ps.print_stats(20)
    print(s.getvalue())


def main():
    print("Loading ifo structure...")
    ifo = gwinc.Struct.from_file(GWINC_PATH)
    gam = brownianProxy(ifo)
    ifo_params = extract_ifo_params(ifo)

    # Pre-populate misc with cached ifo_params
    misc['_ifo_params'] = ifo_params

    level1_instrumented(ifo, gam, ifo_params)
    level2_cprofile(ifo, gam, ifo_params)
    level3_devo_overhead(ifo, gam, ifo_params, popsize=200)

    print("\n" + "=" * 70)
    print("PROFILING COMPLETE")
    print("=" * 70)


if __name__ == '__main__':
    main()
