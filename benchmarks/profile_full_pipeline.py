"""End-to-end profiling: full optimization + MCMC pipeline.

Run from the repo root:
    cd SiN_aSi && python ../benchmarks/profile_full_pipeline.py

Instruments wall-clock time for each phase and prints a summary table.
"""
import os
import sys
import warnings
from timeit import default_timer

import numpy as np

# Ensure repo root is on path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
os.chdir(os.path.join(os.path.dirname(__file__), '..', 'SiN_aSi'))

from scipy.optimize import differential_evolution as devo
import gwinc
import emcee

from generic.optimUtils import getMirrorCost, brownianProxy, precompute_misc
from generic.coatingUtils import importParams, multidiel1, lnprob, op2phys, surfaceField


def run_optimization():
    """Run full differential_evolution optimization, return timing + results."""
    opt_params = importParams('ETM_params.yml')
    ifo = gwinc.Struct.from_file(opt_params['misc']['gwincStructFile'])
    gam = brownianProxy(ifo)

    Ls = np.ones(2 * opt_params['misc']['Npairs'] + 1)
    Ls *= np.random.rand(len(Ls))

    bounds = ((0.05, 0.48),) * (len(Ls) - 1)
    minThick = 10e-9 / ifo.Laser.Wavelength
    bounds = ((minThick, 0.48),) + bounds

    # Pre-compute cached values
    precompute_misc(opt_params['costs'], ifo, opt_params['misc'])

    # Warm up Numba JIT
    _ = getMirrorCost(Ls, opt_params['costs'], ifo, gam,
                      verbose=False, misc=opt_params['misc'])

    tic = default_timer()
    res = devo(
        func=getMirrorCost,
        bounds=bounds,
        updating='deferred',
        strategy='best1bin',
        mutation=(0.05, 1.5),
        popsize=opt_params['misc']['Nparticles'],
        init=opt_params['misc']['init_method'],
        workers=1,
        maxiter=2000,
        atol=opt_params['misc']['atol'],
        tol=opt_params['misc']['tol'],
        args=(opt_params['costs'], ifo, gam, False, opt_params['misc']),
        polish=True,
        disp=False,
    )
    dt_opt = default_timer() - tic

    if not res.success:
        warnings.warn(f"Optimizer did not converge: {res.message}")

    # Get final cost
    Lres = res.x
    copiedLayers = np.tile(
        Lres[1:2*opt_params['misc']['Npairs']+1].copy(),
        opt_params['misc']['Ncopies'])
    Lres = np.append(Lres, copiedLayers)
    fixedLayers = np.tile(Lres[-2:].copy(), opt_params['misc']['Nfixed'])
    Lres = np.append(Lres, fixedLayers)

    final_misc = dict(opt_params['misc'], Ncopies=0, Nfixed=0)
    precompute_misc(opt_params['costs'], ifo, final_misc)
    scalar_cost, output = getMirrorCost(
        Lres, costs=opt_params['costs'], ifo=ifo, gam=gam,
        verbose=True, misc=final_misc)

    return {
        'dt': dt_opt,
        'nfev': res.nfev,
        'nit': res.nit,
        'cost': scalar_cost,
        'output': output,
        'ifo': ifo,
        'L_opt': output['L'],
        'n_IR': output['n'],
        'gwincFile': opt_params['misc']['gwincStructFile'],
        'lambdaAUX': opt_params['misc'].get('lambdaAUX', 1550/2050),
    }


def run_mcmc(opt_result, nSamples=5000):
    """Run MCMC perturbation analysis, return timing breakdown."""
    n_IR_out = opt_result['n_IR']
    L_opt = opt_result['L_opt']
    L_out = op2phys(L_opt, n_IR_out[1:-1])
    lambdaAUX = opt_result['lambdaAUX']

    # --- Sampler setup ---
    tic_sampler = default_timer()
    nWalkers, nDim = 20, 4
    means = np.zeros(nDim)
    cov = np.diag(0.005 * np.ones(nDim))
    cov = np.dot(cov, cov)
    icov = np.linalg.inv(cov)
    p0 = np.random.rand(nDim * nWalkers).reshape((nWalkers, nDim))
    sampler = emcee.EnsembleSampler(nWalkers, nDim, lnprob, args=[means, icov])
    pos, prob, state = sampler.run_mcmc(p0, 1000)
    sampler.reset()
    pos_final, prob_final, state_final = sampler.run_mcmc(pos, 5000)
    dt_sampler = default_timer() - tic_sampler

    # --- Array pre-computation ---
    tic_precomp = default_timer()
    perturbs = sampler.flatchain[:nSamples, :]
    perturb_factors = 1 + perturbs

    n_IR_all = np.tile(n_IR_out, (nSamples, 1))
    n_IR_all[:, 1:-1:2] *= perturb_factors[:, 2:3]
    n_IR_all[:, 2::2] *= perturb_factors[:, 1:2]
    L_all = L_out[None, :] * perturb_factors[:, 3:4]
    dt_precomp = default_timer() - tic_precomp

    # --- Evaluation loop ---
    Tp_IR = np.empty(nSamples)
    Tp_AUX = np.empty(nSamples)
    surfField_arr = np.empty(nSamples)

    tic_eval = default_timer()
    for jj in range(nSamples):
        n_IRs = n_IR_all[jj]
        Ls = L_all[jj]
        Gamma5p, _ = multidiel1(n_IRs, Ls * n_IRs[1:-1], 1)
        Tp_IR[jj] = 1.0 - np.abs(Gamma5p)**2
        surfField_arr[jj] = surfaceField(Gamma5p)
        Gamma5p, _ = multidiel1(n_IRs, Ls * n_IRs[1:-1], lambdaAUX)
        Tp_AUX[jj] = 1 - np.abs(Gamma5p)**2
    dt_eval = default_timer() - tic_eval

    return {
        'dt_sampler': dt_sampler,
        'dt_precomp': dt_precomp,
        'dt_eval': dt_eval,
        'dt_total': dt_sampler + dt_precomp + dt_eval,
        'nSamples': nSamples,
        'Tp_IR_mean': np.mean(Tp_IR),
        'Tp_AUX_mean': np.mean(Tp_AUX),
    }


def main():
    np.random.seed(42)
    print("=" * 60)
    print("Full Pipeline Profiling")
    print("=" * 60)

    # --- Optimization ---
    print("\n[1/2] Running full optimization (popsize=200, workers=1)...")
    opt = run_optimization()
    throughput = opt['nfev'] / opt['dt']
    print(f"  Wall time:    {opt['dt']:.1f} sec")
    print(f"  Func evals:   {opt['nfev']}")
    print(f"  Iterations:   {opt['nit']}")
    print(f"  Throughput:   {throughput:,.0f} evals/sec")
    print(f"  Final cost:   {opt['cost']:.6f}")

    # --- MCMC ---
    print("\n[2/2] Running MCMC (5000 samples)...")
    mc = run_mcmc(opt, nSamples=5000)
    print(f"  Sampler setup: {mc['dt_sampler']:.2f} sec")
    print(f"  Array precomp: {mc['dt_precomp']:.4f} sec")
    print(f"  Eval loop:     {mc['dt_eval']:.2f} sec")
    print(f"  Total MCMC:    {mc['dt_total']:.2f} sec")
    print(f"  Per-sample:    {mc['dt_eval']/mc['nSamples']*1e6:.0f} us")

    # --- Summary ---
    total = opt['dt'] + mc['dt_total']
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"{'Phase':<25} {'Time (s)':>10} {'Fraction':>10}")
    print("-" * 45)
    print(f"{'Optimization':<25} {opt['dt']:>10.2f} {opt['dt']/total:>10.1%}")
    print(f"{'MCMC sampler setup':<25} {mc['dt_sampler']:>10.2f} {mc['dt_sampler']/total:>10.1%}")
    print(f"{'MCMC array precomp':<25} {mc['dt_precomp']:>10.4f} {mc['dt_precomp']/total:>10.1%}")
    print(f"{'MCMC eval loop':<25} {mc['dt_eval']:>10.2f} {mc['dt_eval']/total:>10.1%}")
    print("-" * 45)
    print(f"{'TOTAL':<25} {total:>10.2f}")
    print(f"\nOptimizer: {throughput:,.0f} evals/sec, "
          f"{opt['dt']/opt['nfev']*1e6:.1f} us/eval")
    print(f"MCMC eval: {mc['dt_eval']/mc['nSamples']*1e6:.0f} us/sample "
          f"({2*mc['nSamples']} multidiel1 calls)")


if __name__ == '__main__':
    main()
