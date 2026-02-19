"""Unit tests for generic/coatingUtils.py"""
import numpy as np
import pytest
from generic.coatingUtils import (
    multidiel1, op2phys, surfaceField, sellmeier, _transfer_matrix_loop,
)


class TestMultidiel1:
    """Tests for the transfer matrix method implementation."""

    def test_quarter_wave_analytical(self, quarter_wave_stack):
        """Quarter-wave stack reflectivity matches analytical formula."""
        n, L, n_low, n_high, n_sub, N = quarter_wave_stack
        r, Z = multidiel1(n, L, np.array([1.0]), 0, 'te')
        R_computed = np.abs(r[0])**2

        # Analytical: R = ((1 - (nH/nL)^(2N) * nS) / (1 + (nH/nL)^(2N) * nS))^2
        ratio = (n_high / n_low)**(2 * N) * n_sub
        R_analytical = ((1 - ratio) / (1 + ratio))**2
        assert abs(R_computed - R_analytical) < 1e-7

    def test_scalar_vs_array_lambda(self, quarter_wave_stack):
        """Scalar and array lambda give the same result."""
        n, L, *_ = quarter_wave_stack
        r_scalar, _ = multidiel1(n, L, 1.0, 0, 'te')
        r_array, _ = multidiel1(n, L, np.array([1.0]), 0, 'te')
        assert np.allclose(r_scalar, r_array)

    def test_no_layers(self):
        """Air-substrate interface (M=0) gives Fresnel reflection."""
        n_sub = 1.45
        n = np.array([1.0, n_sub])
        L = np.array([])
        r, _ = multidiel1(n, L, np.array([1.0]), 0, 'te')
        r_fresnel = (1.0 - n_sub) / (1.0 + n_sub)
        assert abs(r[0] - r_fresnel) < 1e-12

    def test_multi_wavelength(self, quarter_wave_stack):
        """Multi-wavelength call produces correct-length output."""
        n, L, *_ = quarter_wave_stack
        lambs = np.array([0.8, 1.0, 1.2])
        r, Z = multidiel1(n, L, lambs, 0, 'te')
        assert r.shape == (3,)
        assert Z.shape == (3,)

    def test_high_reflectivity_at_design(self, quarter_wave_stack):
        """10-bilayer QW stack should have very high reflectivity at design wavelength."""
        n, L, *_ = quarter_wave_stack
        r, _ = multidiel1(n, L, np.array([1.0]), 0, 'te')
        R = np.abs(r[0])**2
        assert R > 0.9999

    def test_lower_reflectivity_off_design(self, quarter_wave_stack):
        """Reflectivity should drop significantly away from design wavelength."""
        n, L, *_ = quarter_wave_stack
        r_off, _ = multidiel1(n, L, np.array([0.5]), 0, 'te')
        R_off = np.abs(r_off[0])**2
        # Far from design wavelength, R should be lower
        assert R_off < 0.99

    def test_te_tm_differ_at_angle(self, quarter_wave_stack):
        """TE and TM polarizations give different results at oblique incidence."""
        n, L, *_ = quarter_wave_stack
        r_te, _ = multidiel1(n, L, np.array([1.0]), 30, 'te')
        r_tm, _ = multidiel1(n, L, np.array([1.0]), 30, 'tm')
        assert not np.allclose(r_te, r_tm)

    def test_normal_incidence_te_tm_equal(self, quarter_wave_stack):
        """TE and TM should be identical at normal incidence."""
        n, L, *_ = quarter_wave_stack
        r_te, _ = multidiel1(n, L, np.array([1.0]), 0, 'te')
        r_tm, _ = multidiel1(n, L, np.array([1.0]), 0, 'tm')
        assert np.allclose(r_te, r_tm, atol=1e-12)


class TestOp2Phys:
    """Tests for optical-to-physical thickness conversion."""

    def test_basic_conversion(self):
        """Known optical thickness -> physical thickness."""
        L = np.array([0.25, 0.25])
        n = np.array([1.5, 3.5])
        phys = op2phys(L, n)
        expected = np.array([0.25 / 1.5, 0.25 / 3.5])
        assert np.allclose(phys, expected)

    def test_dimension_mismatch(self):
        """Mismatched L and n should raise ValueError."""
        with pytest.raises(ValueError):
            op2phys(np.array([0.25, 0.25]), np.array([1.5]))

    def test_uniform_index(self):
        """With n=1, physical = optical."""
        L = np.array([0.3, 0.5, 0.7])
        n = np.ones(3)
        assert np.allclose(op2phys(L, n), L)


class TestSurfaceField:
    """Tests for surface E-field calculation."""

    def test_perfect_ar(self):
        """For r=0 (perfect AR), surface field = Ei * |1+0| = Ei."""
        Ei = 27.46
        sf = surfaceField(0.0, Ei)
        assert abs(sf - Ei) < 1e-10

    def test_perfect_hr(self):
        """For r=-1 (perfect HR), surface field = Ei * |1-1| = 0."""
        sf = surfaceField(-1.0, 27.46)
        assert abs(sf) < 1e-10

    def test_partial_reflector(self):
        """For r=0.5, surface field = Ei * |1.5| = 1.5*Ei."""
        Ei = 27.46
        sf = surfaceField(0.5, Ei)
        assert abs(sf - 1.5 * Ei) < 1e-10


class TestSellmeier:
    """Tests for Sellmeier dispersion equation."""

    def test_fused_silica_1064nm(self):
        """SiO2 at 1064nm should give n ~ 1.4496."""
        n = sellmeier(lam=1064e-9)
        assert abs(n - 1.4496) < 0.001

    def test_wavelength_dependence(self):
        """Index should decrease with increasing wavelength (normal dispersion)."""
        n_short = sellmeier(lam=500e-9)
        n_long = sellmeier(lam=1500e-9)
        assert n_short > n_long

    def test_coefficient_mismatch(self):
        """Mismatched B and C should raise ValueError."""
        with pytest.raises(ValueError):
            sellmeier(B=[0.5, 0.3], C=[0.01])
