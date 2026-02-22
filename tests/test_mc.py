"""Tests for OptimalBragg.mc — Monte Carlo sensitivity analysis."""

import numpy as np
import pytest


class TestGeneratePerturbations:
    def test_shape(self):
        from OptimalBragg.mc import _generate_perturbations
        perturbs = _generate_perturbations(100, n_dim=3)
        assert perturbs.shape == (100, 3)

    def test_centered_near_zero(self):
        from OptimalBragg.mc import _generate_perturbations
        perturbs = _generate_perturbations(500, n_dim=3)
        # Mean should be near zero (within a few percent)
        assert np.abs(perturbs.mean()) < 0.05

    def test_width(self):
        from OptimalBragg.mc import _generate_perturbations
        perturbs = _generate_perturbations(1000, n_dim=3, width=0.005)
        # Std should be approximately 0.005
        for dim in range(3):
            assert 0.001 < np.std(perturbs[:, dim]) < 0.02


class TestSaveMC:
    def test_roundtrip(self, tmp_path):
        from OptimalBragg.mc import save_mc
        import h5py

        result = {
            'MCout': np.random.randn(5, 100),
            'TOnoise': np.random.randn(100, 50),
            'Brnoise': np.random.randn(100, 50),
            'freq': np.logspace(1, 3, 50),
            'perturbs': np.random.randn(100, 3),
        }
        path = str(tmp_path / 'mc.hdf5')
        save_mc(result, path)

        with h5py.File(path, 'r') as f:
            assert 'MCout' in f
            assert f['MCout'].shape == (5, 100)
            assert 'perturbs' in f
            assert f['perturbs'].shape == (100, 3)
