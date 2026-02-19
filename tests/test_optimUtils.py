"""Unit tests for generic/optimUtils.py"""
import numpy as np
import pytest
from generic.optimUtils import (
    transmissionCost, brownianProxy, brownianCost,
    stdevLCost, getMirrorCost, surfEfieldCost,
)


class TestTransmissionCost:
    """Tests for the transmission cost function."""

    def test_returns_scalar(self, ifo_voyager):
        """Cost should be a scalar float, not an array."""
        n_low = ifo_voyager.Materials.Coating.Indexlown
        n_high = ifo_voyager.Materials.Coating.Indexhighn
        n_sub = ifo_voyager.Materials.Substrate.RefractiveIndex
        n = np.array([1.0, n_low, n_high, n_low, n_high, n_sub])
        L = 0.25 * np.ones(4)
        cost, T = transmissionCost(0.5, n, L, 1.0, 0, 'te')
        assert np.isscalar(cost) or cost.ndim == 0

    def test_zero_cost_at_target(self):
        """If the stack T matches target exactly, cost ~ 0."""
        # Single slab: compute actual T, then use it as target
        n = np.array([1.0, 1.45, 1.0])
        L = np.array([0.25])
        _, T = transmissionCost(0.5, n, L, 1.0, 0, 'te')
        cost, _ = transmissionCost(float(T), n, L, 1.0, 0, 'te')
        assert cost < 1e-20


class TestBrownianProxy:
    """Tests for the Brownian noise proxy."""

    def test_returns_positive(self, ifo_voyager):
        """Brownian proxy gamma should be positive."""
        gam = brownianProxy(ifo_voyager)
        assert gam > 0

    def test_proportional_to_thickness(self, ifo_voyager):
        """Doubling layer thicknesses should roughly double Brownian cost."""
        gam = brownianProxy(ifo_voyager)
        L = 0.25 * np.ones(10)
        cost1 = brownianCost(1.0, L, gam)
        cost2 = brownianCost(1.0, 2 * L, gam)
        assert abs(cost2 / cost1 - 2.0) < 0.01


class TestStdevLCost:
    """Tests for the layer thickness standard deviation cost."""

    def test_uniform_layers(self):
        """Uniform thickness -> std=0 -> relative_stdev=0 -> high cost."""
        L = 0.25 * np.ones(10)
        cost = stdevLCost(1.0, L)
        # std is zero, relative_stdev = 0, cost = |1-0|/1 = 1
        assert abs(cost - 1.0) < 1e-10

    def test_varied_layers(self):
        """Varied thickness -> nonzero relative_stdev -> cost depends on target."""
        L = np.array([0.1, 0.3, 0.1, 0.3, 0.1, 0.3])
        # mean/std = 0.2/~0.1095 = ~1.826; target=1.826 -> cost~0
        rel_std = np.mean(L) / np.std(L)
        cost = stdevLCost(rel_std, L)
        assert cost < 1e-10


class TestGetMirrorCost:
    """Tests for the master cost function."""

    def test_smoke_test(self, ifo_voyager, sample_costs, sample_misc):
        """Basic call returns a finite scalar."""
        gam = brownianProxy(ifo_voyager)
        L = 0.25 * np.ones(10)
        cost = getMirrorCost(L, sample_costs, ifo_voyager, gam,
                             verbose=False, misc=sample_misc)
        assert np.isfinite(cost)

    def test_verbose_returns_dict(self, ifo_voyager, sample_costs, sample_misc):
        """Verbose mode returns (scalar, dict) with expected keys."""
        gam = brownianProxy(ifo_voyager)
        L = 0.25 * np.ones(10)
        result = getMirrorCost(L, sample_costs, ifo_voyager, gam,
                               verbose=True, misc=sample_misc)
        assert isinstance(result, tuple) and len(result) == 2
        scalar_cost, output = result
        assert np.isfinite(scalar_cost)
        assert 'n' in output
        assert 'L' in output
        assert 'vectorCost' in output
        assert 'TPSL' in output

    def test_zero_weights_zero_cost(self, ifo_voyager, sample_misc):
        """All weights=0 should give cost=0."""
        gam = brownianProxy(ifo_voyager)
        zero_costs = {
            'TransPSL':    {'target': 5e-6,   'weight': 0},
            'Brownian':    {'target': 0.1,    'weight': 0},
            'Thermooptic': {'target': 1e43,   'weight': 0},
            'Lsens':       {'target': 1e-7,   'weight': 0},
            'Esurf':       {'target': 1e-9,   'weight': 0},
            'Absorption':  {'target': 1e-4,   'weight': 0},
            'TransAUX':    {'target': 1000e-6, 'weight': 0},
            'TransOPLEV':  {'target': 0.05,   'weight': 0},
            'Lstdev':      {'target': 0.5,    'weight': 0},
        }
        L = 0.25 * np.ones(10)
        cost = getMirrorCost(L, zero_costs, ifo_voyager, gam,
                             verbose=False, misc=sample_misc)
        assert cost == 0.0

    def test_ncopies_extends_stack(self, ifo_voyager, sample_costs, sample_misc):
        """Ncopies > 0 should produce a longer L in verbose output."""
        gam = brownianProxy(ifo_voyager)
        L = 0.25 * np.ones(10)
        misc_copy = dict(sample_misc, Ncopies=1)
        _, output = getMirrorCost(L, sample_costs, ifo_voyager, gam,
                                  verbose=True, misc=misc_copy)
        # With Ncopies=1, the L in output should be longer than input
        assert len(output['L']) > len(L)
