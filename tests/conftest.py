import pytest
import numpy as np


@pytest.fixture
def quarter_wave_stack():
    """10-bilayer quarter-wave stack with known analytical reflectivity."""
    n_low, n_high, n_sub = 1.45, 3.5, 1.45
    N = 10
    n = np.array([1.0] + [n_low, n_high] * N + [n_sub])
    L = 0.25 * np.ones(2 * N)
    return n, L, n_low, n_high, n_sub, N


@pytest.fixture
def sample_costs():
    """Minimal cost dict for testing getMirrorCost."""
    return {
        'Trans1':      {'target': 5e-6,  'weight': 5},
        'Brownian':    {'target': 0.1,   'weight': 2},
        'Thermooptic': {'target': 1e43,  'weight': 0},
        'Lsens':       {'target': 1e-7,  'weight': 0},
        'Esurf':       {'target': 1e-9,  'weight': 0},
        'Absorption':  {'target': 1e-4,  'weight': 0},
        'Trans2':      {'target': 1000e-6, 'weight': 0},
        'Trans3':      {'target': 0.05,  'weight': 0},
        'Lstdev':      {'target': 0.5,   'weight': 0},
    }


@pytest.fixture
def sample_misc():
    """Minimal misc dict for testing getMirrorCost."""
    return {
        'fTO': 100, 'pol': 'te', 'aoi': 0,
        'Npairs': 5, 'Nfixed': 0, 'Ncopies': 0,
        'lambda2': 0.756,
    }
