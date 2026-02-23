"""Tests for OptimalBragg.costs — cost functions and master evaluator."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from OptimalBragg import Material, qw_stack
from OptimalBragg.materials import SiO2, TiTa2O5, FusedSilica, air


# ── Fixtures ─────────────────────────────────────────────────────────

@pytest.fixture
def aLIGO_stack():
    """18-bilayer SiO2/TiTa2O5 QW stack (aLIGO-like)."""
    return qw_stack(
        lam_ref=1064e-9,
        substrate=Material(FusedSilica),
        superstrate=Material(air),
        thin_films={"L": Material(SiO2), "H": Material(TiTa2O5)},
        pattern="LH" * 18,
    )


# ── Normalization decorator ──────────────────────────────────────────

class TestNormDecorator:
    def test_l1_identity(self):
        from OptimalBragg.costs import norm

        @norm("l1")
        def f():
            return 3.0
        assert f() == 3.0

    def test_l2_squares(self):
        from OptimalBragg.costs import norm

        @norm("l2")
        def f():
            return 3.0
        assert f() == 9.0

    def test_arcsinh(self):
        from OptimalBragg.costs import norm

        @norm("arcsinh")
        def f():
            return 3.0
        assert_allclose(f(), np.arcsinh(3.0))


# ── Individual cost functions ────────────────────────────────────────

class TestTransmissionCost:
    def test_returns_scalar(self):
        from OptimalBragg.costs import transmissionCost
        n = np.array([1.0, 1.45, 2.06, 1.45])
        L = np.array([0.25, 0.25])
        cost, T = transmissionCost(0.5, n, L, 1.0, 0, 'te')
        assert np.isscalar(cost)

    def test_zero_cost_at_target(self):
        """Cost should be near zero when T matches target."""
        from OptimalBragg.costs import transmissionCost
        n = np.array([1.0, 1.45, 2.06, 1.45])
        L = np.array([0.25, 0.25])
        _, T = transmissionCost(0.5, n, L, 1.0, 0, 'te')
        cost, _ = transmissionCost(float(T), n, L, 1.0, 0, 'te')
        assert cost < 1e-10


class TestSensitivityCost:
    def test_nonnegative(self):
        from OptimalBragg.costs import sensitivityCost
        n = np.array([1.0, 1.45, 2.06, 1.45])
        L = np.array([0.25, 0.25])
        cost = sensitivityCost(0.5, n, L, 1.0)
        assert cost >= 0


class TestSurfEfieldCost:
    def test_positive(self):
        from OptimalBragg.costs import surfEfieldCost
        n = np.array([1.0, 1.45, 2.06, 1.45])
        L = np.array([0.25, 0.25])
        cost = surfEfieldCost(1.0, n, L, 1.0)
        assert cost > 0


class TestStdevLCost:
    def test_uniform_layers(self):
        from OptimalBragg.costs import stdevLCost
        L = np.array([0.25, 0.25, 0.25, 0.25])
        cost = stdevLCost(10.0, L)
        # Uniform → std=0 → relative_stdev=0, so cost = |10-0|^2/10^2 = 1
        assert_allclose(cost, 1.0)

    def test_varied_layers(self):
        from OptimalBragg.costs import stdevLCost
        L = np.array([0.20, 0.30, 0.20, 0.30])
        cost = stdevLCost(10.0, L)
        assert cost > 0


class TestBrownianCost:
    def test_zero_below_budget(self):
        from OptimalBragg.costs import brownianCost
        L = np.array([0.25, 0.25] * 5)
        cost = brownianCost(100.0, L, 1.5)  # well under budget
        assert cost == 0.0

    def test_positive_above_budget(self):
        from OptimalBragg.costs import brownianCost
        L = np.array([0.25, 0.25] * 18)
        cost = brownianCost(1.0, L, 1.5)  # tight budget
        assert cost > 0

    def test_proportional_to_thickness(self):
        from OptimalBragg.costs import brownianCost
        # Both above budget so we can compare
        L1 = np.array([0.25, 0.25] * 10)
        L2 = np.array([0.25, 0.25] * 20)
        c1 = brownianCost(1.0, L1, 1.5)
        c2 = brownianCost(1.0, L2, 1.5)
        assert c2 > c1

    def test_accepts_dict_gam(self):
        from OptimalBragg.costs import brownianCost
        L = np.array([0.25, 0.25] * 10)
        gam_dict = {"H": 1.5}
        cost = brownianCost(1.0, L, gam_dict)  # tight budget
        assert cost > 0


class TestThermoopticCost:
    def test_positive_above_budget(self, aLIGO_stack):
        from OptimalBragg.costs import thermoopticCost
        from OptimalBragg.noise import extract_stack_params
        params = extract_stack_params(aLIGO_stack, 0.17, 0.20)
        L = aLIGO_stack["Ls_opt"]
        # Very tight target to ensure above budget
        cost = thermoopticCost(1e-50, 100.0, L, aLIGO_stack,
                               stack_params=params, w_beam=0.062)
        assert cost > 0


# ── Precompute misc ──────────────────────────────────────────────────

class TestPrecomputeMisc:
    def test_populates_cache(self, aLIGO_stack):
        from OptimalBragg.costs import precompute_misc
        costs = {
            'Trans1064': {'target': 5e-6, 'weight': 15},
            'Trans532': {'target': 0.032, 'weight': 5},
            'Brownian': {'target': 20.0, 'weight': 2},
            'Thermooptic': {'target': 1.6e-42, 'weight': 2},
        }
        misc = {
            'Npairs': 18, 'aoi': 0, 'pol': 'te',
            'lambdaAUX': 0.5, 'fTO': 100,
            'r_mirror': 0.17, 'd_mirror': 0.20,
        }
        precompute_misc(costs, aLIGO_stack, misc)
        assert '_ifo_n' in misc
        assert '_active_costs' in misc
        assert '_wl_array' in misc
        assert '_stack_params' in misc
        assert len(misc['_ifo_n']) == 2 * 18 + 2  # layers + vac + sub


# ── Master evaluator ─────────────────────────────────────────────────

class TestGetMirrorCost:
    @pytest.fixture
    def setup(self, aLIGO_stack):
        from OptimalBragg.costs import precompute_misc
        from OptimalBragg.noise import brownian_proxy
        costs = {
            'Trans1064': {'target': 5e-6, 'weight': 15},
            'Brownian': {'target': 20.0, 'weight': 2},
        }
        gam = brownian_proxy(aLIGO_stack)
        misc = {
            'Npairs': 18, 'aoi': 0, 'pol': 'te',
            'lambdaAUX': 0.5, 'fTO': 100,
            'r_mirror': 0.17, 'd_mirror': 0.20,
            'w_beam': 0.062,
        }
        precompute_misc(costs, aLIGO_stack, misc)
        return costs, aLIGO_stack, gam, misc

    def test_returns_scalar(self, setup):
        from OptimalBragg.costs import getMirrorCost
        costs, stack, gam, misc = setup
        L = stack["Ls_opt"]
        cost = getMirrorCost(L, costs, stack, gam, misc=misc)
        assert np.isscalar(cost)
        assert cost > 1.0  # multiplicative product always >= 1

    def test_verbose_returns_dict(self, setup):
        from OptimalBragg.costs import getMirrorCost
        costs, stack, gam, misc = setup
        L = stack["Ls_opt"]
        cost, out = getMirrorCost(L, costs, stack, gam, verbose=True,
                                  misc=misc)
        assert 'vectorCost' in out
        assert 'Trans1064' in out['vectorCost']
        assert 'Brownian' in out['vectorCost']
        assert 'T1064' in out
        assert 'R1064' in out

    def test_zero_weights_unity_cost(self, aLIGO_stack):
        """All zero weights → scalar cost should be 1.0."""
        from OptimalBragg.costs import getMirrorCost, precompute_misc
        costs = {
            'Trans1064': {'target': 5e-6, 'weight': 0},
            'Brownian': {'target': 20.0, 'weight': 0},
        }
        misc = {
            'Npairs': 18, 'aoi': 0, 'pol': 'te',
            'lambdaAUX': 0.5,
        }
        precompute_misc(costs, aLIGO_stack, misc)
        L = aLIGO_stack["Ls_opt"]
        cost = getMirrorCost(L, costs, aLIGO_stack, 1.5, misc=misc)
        assert_allclose(cost, 1.0)

    def test_with_thermooptic(self, aLIGO_stack):
        from OptimalBragg.costs import getMirrorCost, precompute_misc
        from OptimalBragg.noise import brownian_proxy
        costs = {
            'Trans1064': {'target': 5e-6, 'weight': 15},
            'Thermooptic': {'target': 1.6e-42, 'weight': 2},
        }
        gam = brownian_proxy(aLIGO_stack)
        misc = {
            'Npairs': 18, 'aoi': 0, 'pol': 'te',
            'lambdaAUX': 0.5, 'fTO': 100,
            'r_mirror': 0.17, 'd_mirror': 0.20, 'w_beam': 0.062,
        }
        precompute_misc(costs, aLIGO_stack, misc)
        L = aLIGO_stack["Ls_opt"]
        cost = getMirrorCost(L, costs, aLIGO_stack, gam, misc=misc)
        assert cost > 1.0

    def test_ncopies_extends_stack(self, aLIGO_stack):
        """Ncopies should tile extra layers."""
        from OptimalBragg.costs import getMirrorCost, precompute_misc
        costs = {'Trans1064': {'target': 5e-6, 'weight': 15}}
        misc_base = {
            'Npairs': 18, 'aoi': 0, 'pol': 'te',
            'lambdaAUX': 0.5,
        }
        misc_copy = {
            'Npairs': 18, 'Ncopies': 1, 'aoi': 0, 'pol': 'te',
            'lambdaAUX': 0.5,
        }
        precompute_misc(costs, aLIGO_stack, misc_base)
        precompute_misc(costs, aLIGO_stack, misc_copy)
        L = aLIGO_stack["Ls_opt"]
        c1 = getMirrorCost(L, costs, aLIGO_stack, 1.5, misc=misc_base)
        c2 = getMirrorCost(L, costs, aLIGO_stack, 1.5, misc=misc_copy)
        # Different number of layers → different cost
        assert c1 != c2
