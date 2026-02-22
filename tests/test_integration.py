"""Integration test: short optimization run using OptimalBragg."""
import numpy as np
import pytest
from scipy.optimize import differential_evolution as devo

from OptimalBragg import Material, qw_stack
from OptimalBragg.materials import SiO2, TiTa2O5, FusedSilica, air
from OptimalBragg.costs import getMirrorCost, precompute_misc
from OptimalBragg.noise import brownian_proxy


@pytest.mark.slow
def test_short_optimization():
    """Run a tiny differential_evolution to verify end-to-end pipeline."""
    Npairs = 3
    stack = qw_stack(
        lam_ref=1064e-9,
        substrate=Material(FusedSilica),
        superstrate=Material(air),
        thin_films={"L": Material(SiO2), "H": Material(TiTa2O5)},
        pattern="LH" * Npairs,
    )
    gam = brownian_proxy(stack)

    costs = {
        'Trans1064': {'target': 5e-6, 'weight': 5},
        'Brownian': {'target': 20.0, 'weight': 2},
    }
    misc = {
        'pol': 'te', 'aoi': 0,
        'Npairs': Npairs, 'Nfixed': 0, 'Ncopies': 0,
        'lambdaAUX': 0.5,
    }

    precompute_misc(costs, stack, misc)

    n_vars = 2 * Npairs
    bounds = ((0.05, 0.45),) * n_vars

    res = devo(
        func=getMirrorCost,
        bounds=bounds,
        updating='deferred',
        strategy='best1bin',
        popsize=5,
        maxiter=10,
        workers=1,
        args=(costs, stack, gam, False, misc),
        polish=False,
        disp=False,
    )

    assert np.isfinite(res.fun)
    assert len(res.x) == n_vars
    assert all(np.isfinite(res.x))
