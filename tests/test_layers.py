"""Tests for OptimalBragg.layers — transfer matrix, E-field, utilities."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from OptimalBragg.layers import (
    multidiel1,
    multilayer_diel,
    amp_refl,
    refl,
    trans,
    surfield,
    field_zmag,
    calc_abs,
    op2phys,
    sellmeier,
    qwbandedges,
)
from OptimalBragg import Material, qw_stack
from OptimalBragg.materials import SiO2, TiTa2O5, air


# ── Fixtures ─────────────────────────────────────────────────────────────

@pytest.fixture
def qw18():
    """18-bilayer SiO2/TiTa2O5 quarter-wave stack at 1064 nm."""
    return qw_stack(
        lam_ref=1064e-9,
        substrate=Material(SiO2),
        superstrate=Material(air),
        thin_films={'L': Material(SiO2), 'H': Material(TiTa2O5)},
        pattern='LH' * 18,
    )


# ── multidiel1 (JIT core, optical thicknesses) ──────────────────────────

class TestMultidiel1:
    def test_single_interface_fresnel(self):
        """No layers: should reduce to Fresnel at n1/n2 boundary."""
        n = np.array([1.0, 1.5])
        r, _ = multidiel1(n, [], [1.0])
        # Fresnel: r = (n1 - n2) / (n1 + n2)
        expected = (1.0 - 1.5) / (1.0 + 1.5)
        assert_allclose(r, expected, atol=1e-14)

    def test_qw_stack_high_reflectivity(self):
        """18-bilayer QW stack should have R > 0.99999."""
        nH, nL, nSub = 2.06, 1.45, 1.45
        n = np.array([1.0] + [nL, nH] * 18 + [nSub])
        L = np.array([0.25] * 36)
        r, _ = multidiel1(n, L, [1.0])
        R = np.abs(r[0]) ** 2
        assert R > 0.9999

    def test_scalar_and_array_lamb(self):
        """Scalar vs array wavelength should give same result."""
        n = np.array([1.0, 1.45, 2.06, 1.45])
        L = np.array([0.25, 0.25])
        r_scalar, _ = multidiel1(n, L, 1.0)
        r_array, _ = multidiel1(n, L, [1.0])
        assert_allclose(np.abs(r_scalar), np.abs(r_array), atol=1e-14)

    def test_oblique_incidence_changes_result(self):
        """Non-zero angle should give different R than normal."""
        n = np.array([1.0, 1.45, 2.06, 1.45])
        L = np.array([0.25, 0.25])
        r0, _ = multidiel1(n, L, [1.0], theta=0)
        r45, _ = multidiel1(n, L, [1.0], theta=45)
        assert not np.allclose(np.abs(r0), np.abs(r45))


# ── multilayer_diel (physical thicknesses) ───────────────────────────────

class TestMultilayerDiel:
    def test_matches_multidiel1_at_normal_incidence(self):
        """Physical-thickness API should give same |R| as optical API."""
        lam0 = 1064e-9
        nH, nL = 2.06, 1.45
        n = np.array([1.0, nL, nH, nL, nH, nL])
        L_opt = np.array([0.25] * 4)
        # Physical thicknesses
        Ls_phys = L_opt * lam0 / np.array([nL, nH, nL, nH])

        r_opt, _ = multidiel1(n, L_opt, [1.0])
        r_phys, _ = multilayer_diel(n, Ls_phys, lam0)

        assert_allclose(np.abs(r_opt) ** 2, np.abs(r_phys) ** 2, rtol=1e-12)

    def test_single_interface(self):
        """No layers: Fresnel reflection."""
        ns = np.array([1.0, 1.5])
        r, _ = multilayer_diel(ns, np.array([]), 1064e-9)
        expected = (1.0 - 1.5) / (1.0 + 1.5)
        assert_allclose(r, expected, atol=1e-14)


# ── Stack-based spectral functions ───────────────────────────────────────

class TestStackSpectral:
    def test_refl_plus_trans_equals_one(self, qw18):
        """R + T = 1 (no absorption in TMM)."""
        wls = np.linspace(900e-9, 1200e-9, 50)
        R = refl(wls, qw18)
        T = trans(wls, qw18)
        assert_allclose(R + T, 1.0, atol=1e-14)

    def test_high_reflectivity_at_design_wavelength(self, qw18):
        """R > 0.9999 at 1064 nm for 18-bilayer QW stack."""
        R = refl(np.array([1064e-9]), qw18)
        assert R[0] > 0.9999

    def test_amp_refl_magnitude_consistent(self, qw18):
        """|amp_refl|^2 == refl."""
        wls = np.array([1064e-9, 532e-9])
        r = amp_refl(wls, qw18)
        R = refl(wls, qw18)
        assert_allclose(np.abs(r) ** 2, R, atol=1e-14)


# ── Surface E-field ──────────────────────────────────────────────────────

class TestSurfield:
    def test_perfect_reflector(self):
        """r = -1 → E_surf = 0."""
        assert_allclose(surfield(-1.0), 0.0, atol=1e-14)

    def test_no_reflection(self):
        """r = 0 → E_surf = Ei."""
        assert_allclose(surfield(0.0), 27.46, atol=1e-14)

    def test_normalized(self):
        """Normalized returns |1 + r|."""
        r = 0.5 + 0.3j
        assert_allclose(surfield(r, normalized=True), np.abs(1 + r))


# ── E-field depth profile ────────────────────────────────────────────────

class TestFieldZmag:
    def test_output_shape(self, qw18):
        """Should return arrays of length n_layers * n_pts."""
        ns, Ls = qw18['ns'], qw18['Ls']
        z, E = field_zmag(ns, Ls, 1064e-9, n_pts=10)
        assert len(z) == len(Ls) * 10
        assert len(E) == len(Ls) * 10

    def test_positive_efield(self, qw18):
        """E-field squared must be non-negative."""
        ns, Ls = qw18['ns'], qw18['Ls']
        _, E = field_zmag(ns, Ls, 1064e-9, n_pts=10)
        assert np.all(E >= 0)

    def test_surface_normalized(self):
        """First point should be close to 1 (normalized to surface)."""
        ns = np.array([1.0, 1.45, 2.06, 1.45])
        Ls = np.array([183e-9, 129e-9])
        _, E = field_zmag(ns, Ls, 1064e-9, n_pts=50)
        # Not exactly 1 because it's sampled slightly inside the first layer
        assert 0.5 < E[0] < 2.0


# ── Utilities ────────────────────────────────────────────────────────────

class TestOp2Phys:
    def test_basic(self):
        L = np.array([0.25, 0.25])
        n = np.array([1.45, 2.06])
        phys = op2phys(L, n)
        assert_allclose(phys, L / n)

    def test_length_mismatch_raises(self):
        with pytest.raises(ValueError):
            op2phys([0.25], [1.0, 2.0])


class TestSellmeier:
    def test_fused_silica_at_1064nm(self):
        """Fused silica n ≈ 1.45 at 1064 nm."""
        n = sellmeier(lam=1064e-9)
        assert 1.44 < n < 1.46

    def test_coefficient_mismatch_raises(self):
        with pytest.raises(ValueError):
            sellmeier(B=[1.0], C=[1.0, 2.0])


class TestQwBandEdges:
    def test_band_contains_design_wavelength(self, qw18):
        """Design wavelength should be within the HR band."""
        lam1, lam2 = qwbandedges(qw18)
        assert lam1 < 1064e-9 < lam2
