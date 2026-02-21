"""Integration test: short optimization run."""
import numpy as np
import pytest
from scipy.optimize import differential_evolution as devo
import gwinc

from generic.optimUtils import getMirrorCost, brownianProxy


@pytest.mark.slow
def test_short_optimization():
    """Run a tiny differential_evolution to verify end-to-end pipeline."""
    ifo = gwinc.Struct.from_file('SiN_aSi/aSiSiN.yaml')
    gam = brownianProxy(ifo)

    Npairs = 3
    costs = {
        'Trans1064':    {'target': 5e-6,    'weight': 5},
        'Brownian':    {'target': 0.1,     'weight': 2},
        'Thermooptic': {'target': 1e43,    'weight': 0},
        'Lsens':       {'target': 1e-7,    'weight': 0},
        'Esurf':       {'target': 1e-9,    'weight': 0},
        'Absorption':  {'target': 1e-4,    'weight': 0},
        'Trans532':    {'target': 1000e-6, 'weight': 0},
        'TransOPLEV':  {'target': 0.05,    'weight': 0},
        'Lstdev':      {'target': 0.5,     'weight': 0},
    }
    misc = {
        'fTO': 100, 'pol': 'te', 'aoi': 0,
        'Npairs': Npairs, 'Nfixed': 0, 'Ncopies': 0,
        'lambdaAUX': 0.756,
    }

    n_vars = 2 * Npairs + 1
    bounds = ((0.05, 0.45),) * n_vars

    res = devo(
        func=getMirrorCost,
        bounds=bounds,
        updating='deferred',
        strategy='best1bin',
        popsize=5,
        maxiter=10,
        workers=1,
        args=(costs, ifo, gam, False, misc),
        polish=False,
        disp=False,
    )

    assert np.isfinite(res.fun)
    assert len(res.x) == n_vars
    assert all(np.isfinite(res.x))
