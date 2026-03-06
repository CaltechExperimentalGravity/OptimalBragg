#!/usr/bin/env python
"""Run 10 optimization trials with varied hyperparameters, then MC each.
Goal: find the design with the best (lowest) 95th percentile R@700nm.
"""
import copy, tempfile, yaml, numpy as np
from pathlib import Path
from OptimalBragg.optimizer import run_optimization
from OptimalBragg.mc import run_mc
from OptimalBragg.io import yamlread

BASE_PARAMS = yamlread("projects/SFG/SFG_params.yml")

# 10 trial configs: vary Npairs, sensitivity weights, Trans2 weight, maxiter
TRIALS = [
    # trial 0: baseline N=14, current weights
    dict(Npairs=14, dTdnH=4, dTdnL=10, dTdd=2, Trans2_w=250, maxiter=8000),
    # trial 1: N=14, boost dTdnH sensitivity penalty
    dict(Npairs=14, dTdnH=10, dTdnL=10, dTdd=2, Trans2_w=250, maxiter=8000),
    # trial 2: N=14, boost all sensitivity weights
    dict(Npairs=14, dTdnH=10, dTdnL=15, dTdd=5, Trans2_w=250, maxiter=8000),
    # trial 3: N=14, relax Trans2 weight, boost sensitivity
    dict(Npairs=14, dTdnH=10, dTdnL=15, dTdd=5, Trans2_w=100, maxiter=8000),
    # trial 4: N=15, boost sensitivity
    dict(Npairs=15, dTdnH=10, dTdnL=15, dTdd=5, Trans2_w=250, maxiter=8000),
    # trial 5: N=15, relax Trans2
    dict(Npairs=15, dTdnH=8, dTdnL=12, dTdd=4, Trans2_w=100, maxiter=8000),
    # trial 6: N=16, high sensitivity penalty
    dict(Npairs=16, dTdnH=12, dTdnL=18, dTdd=6, Trans2_w=250, maxiter=8000),
    # trial 7: N=14, very high sensitivity, longer run
    dict(Npairs=14, dTdnH=15, dTdnL=20, dTdd=8, Trans2_w=250, maxiter=12000),
    # trial 8: N=14, moderate sensitivity, lower Trans2
    dict(Npairs=14, dTdnH=8, dTdnL=12, dTdd=4, Trans2_w=150, maxiter=10000),
    # trial 9: N=13, fewer layers, high sensitivity
    dict(Npairs=13, dTdnH=10, dTdnL=15, dTdd=5, Trans2_w=250, maxiter=8000),
]

results = []

for i, trial in enumerate(TRIALS):
    print(f"\n{'='*60}")
    print(f"Trial {i}: Npairs={trial['Npairs']}, dTdnH={trial['dTdnH']}, "
          f"dTdnL={trial['dTdnL']}, dTdd={trial['dTdd']}, "
          f"Trans2_w={trial['Trans2_w']}, maxiter={trial['maxiter']}")
    print('='*60)

    # Build modified params
    p = copy.deepcopy(BASE_PARAMS)
    p['misc']['Npairs'] = trial['Npairs']
    p['misc']['maxiter'] = trial['maxiter']
    p['costs']['dTdnH']['weight'] = trial['dTdnH']
    p['costs']['dTdnL']['weight'] = trial['dTdnL']
    p['costs']['dTdd']['weight'] = trial['dTdd']
    p['costs']['Trans2']['weight'] = trial['Trans2_w']

    # Write temp YAML
    tmp = Path(f"projects/SFG/_trial{i}_params.yml")
    with open(tmp, 'w') as f:
        yaml.dump(p, f, default_flow_style=False)

    # Run optimization
    try:
        res = run_optimization(str(tmp), save=True, optic='SFG')
        opt_hdf5 = res.get('hdf5_path')

        scalar_cost = res['scalar_cost']
        print(f"  Scalar cost: {scalar_cost:.6f}")

        # Run MC
        mc_result = run_mc(opt_hdf5, n_samples=1000)
        mc = mc_result['MCout']
        # rows: 0=T@1064[ppm], 1=R@700[ppm], 2=R@2050[ppm],
        #        3=S_TO[x1e-21], 4=S_Br[x1e-21], 5=Esurf
        r700 = mc[1]
        r700_median = np.median(r700)
        r700_95 = np.percentile(r700, 95)
        r700_5 = np.percentile(r700, 5)
        t1064 = mc[0]
        t1064_median = np.median(t1064)
        t1064_95 = np.percentile(t1064, 95)
        r2050 = mc[2]
        r2050_median = np.median(r2050)
        r2050_95 = np.percentile(r2050, 95)

        results.append({
            'trial': i,
            'Npairs': trial['Npairs'],
            'cost': scalar_cost,
            'r700_median': r700_median,
            'r700_5': r700_5,
            'r700_95': r700_95,
            't1064_median': t1064_median,
            't1064_95': t1064_95,
            'r2050_median': r2050_median,
            'r2050_95': r2050_95,
            'hdf5': opt_hdf5,
            'params': trial,
        })
        print(f"  R@700:  median={r700_median:.0f} ppm, "
              f"[5th={r700_5:.0f}, 95th={r700_95:.0f}] ppm")
        print(f"  T@1064: median={t1064_median:.0f} ppm, 95th={t1064_95:.0f} ppm")
        print(f"  R@2050: median={r2050_median:.0f} ppm, 95th={r2050_95:.0f} ppm")

    except Exception as e:
        print(f"  FAILED: {e}")
        results.append({'trial': i, 'Npairs': trial['Npairs'],
                        'cost': np.inf, 'r700_95': np.inf,
                        'params': trial, 'error': str(e)})

# Summary table sorted by R@700 95th percentile
print(f"\n\n{'='*80}")
print("SUMMARY — sorted by R@700 95th percentile (lower is better)")
print('='*80)
print(f"{'Trial':>5} {'N':>3} {'Cost':>8} {'T1064_med':>10} {'T1064_95':>9} "
      f"{'R700_med':>9} {'R700_95':>8} {'R2050_95':>9}  Params")
print('-'*80)

for r in sorted(results, key=lambda x: x.get('r700_95', np.inf)):
    if 'error' in r:
        print(f"{r['trial']:>5} {r['Npairs']:>3}  FAILED: {r['error']}")
        continue
    p = r['params']
    print(f"{r['trial']:>5} {r['Npairs']:>3} {r['cost']:>8.4f} "
          f"{r['t1064_median']:>10.0f} {r['t1064_95']:>9.0f} "
          f"{r['r700_median']:>9.0f} {r['r700_95']:>8.0f} {r['r2050_95']:>9.0f}  "
          f"dTdnH={p['dTdnH']} dTdnL={p['dTdnL']} dTdd={p['dTdd']} Tw2={p['Trans2_w']}")
