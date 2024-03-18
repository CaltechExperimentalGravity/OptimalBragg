""" MC analysis on a given coating design"""
import numpy as np
import emcee
import h5py
import tqdm

from physunits.other import ppm
from .layers import stack_refl, stack_R, surfield, field_zmag, calc_abs
from .noise import stack_noise


def lnprob(x, mu, icov):
    diff = x - mu
    return -np.dot(diff, np.dot(icov, diff)) / 2.0


def mc_sens(stack, n_points=5000, n_walkers=20, n_dim=3, **noise_pars):
    # Initialize the emcee sampler, zero mean Gaussians with 0.5% stdev
    means = np.zeros(n_dim)
    cov_0 = np.diag(0.01 * np.ones(n_dim))
    cov = np.dot(cov_0, cov_0)
    inv_cov = np.linalg.inv(cov)
    p0 = np.random.rand(n_dim * n_walkers).reshape((n_walkers, n_dim))
    sampler = emcee.EnsembleSampler(
        n_walkers, n_dim, lnprob, args=[means, inv_cov]
    )
    pos, prob, state = sampler.run_mcmc(p0, 1000)
    sampler.reset()
    pos_final, prob_final, state_final = sampler.run_mcmc(pos, 5000)

    # Use the generated perturbations to run the MCMC
    T_refs = np.ones(n_points)
    Esurfs = np.ones(n_points)
    to_noises = np.ones(n_points)
    br_noises = np.ones(n_points)
    absorpts = np.ones(n_points)
    perturbs = sampler.flatchain[0:n_points, :]

    ns_out = np.copy(stack["ns"])
    Ls_out = np.copy(stack["Ls"])
    als_out = np.copy(stack["alphas"])

    for jj in tqdm.tqdm(range(n_points)):
        # Fractional perturbation
        perturb = 1 + sampler.flatchain[jj, :]

        # Compute reflectivities, surface field and noises @ 100 Hz
        stack["ns"] = np.copy(ns_out) * perturb[0]
        stack["Ls"] = np.copy(Ls_out) * perturb[1]
        stack["alphas"] = np.copy(als_out) * perturb[2]

        rr = stack_refl(stack["lam_ref"], stack)
        T_refs[jj] = 1 - stack_R(stack["lam_ref"], stack)
        Esurfs[jj] = surfield(rr)
        _, Enorm = field_zmag(
            stack["ns"], stack["Ls"], n_pts=2**6, lam=stack["lam_ref"]
        )
        absorpts[jj] = calc_abs(Enorm, stack["Ls"], stack["alphas"])
        Sbr, _, _, Sto = stack_noise(
            freq=np.array([100, 100]), stack=stack, **noise_pars
        )
        to_noises[jj] = np.sqrt(Sto[0])
        br_noises[jj] = np.sqrt(Sbr[0])

    samples = np.vstack(
        (
            (T_refs) / ppm,
            # stack["Ls"].std(),
            to_noises / 1e-21,
            br_noises / 1e-21,
            # Esurfs,
            absorpts / ppm,
        )
    )
    return samples
