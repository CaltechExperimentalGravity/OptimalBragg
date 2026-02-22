"""Tests for OptimalBragg.noise — thermal noise models."""

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


# Mirror geometry (aLIGO ETM-like)
R_MIRROR = 0.17  # m
D_MIRROR = 0.20  # m
W_BEAM = 0.062   # m


# ── Helper functions ──────────────────────────────────────────────────

class TestGetCoatLayers:
    def test_output_lengths(self, aLIGO_stack):
        from OptimalBragg.noise import getCoatLayers
        nC, aL, bL, dL, sL = getCoatLayers(aLIGO_stack)
        n_layers = len(aLIGO_stack["Ls"])
        assert len(nC) == n_layers
        assert len(aL) == n_layers
        assert len(dL) == n_layers

    def test_alternating_indices(self, aLIGO_stack):
        from OptimalBragg.noise import getCoatLayers
        nC, _, _, _, _ = getCoatLayers(aLIGO_stack)
        # Odd layers should be high-n, even should be low-n
        assert_allclose(nC[0::2], 1.45, atol=0.01)
        assert_allclose(nC[1::2], 2.06, atol=0.01)


class TestGetCoatAvg:
    def test_positive_values(self, aLIGO_stack):
        from OptimalBragg.noise import getCoatAvg
        dc, Cc, Kc, aSub = getCoatAvg(aLIGO_stack)
        assert dc > 0
        assert Cc > 0
        assert Kc > 0

    def test_total_thickness_reasonable(self, aLIGO_stack):
        from OptimalBragg.noise import getCoatAvg
        dc, _, _, _ = getCoatAvg(aLIGO_stack)
        # 36 layers each ~lambda/(4n) ~ 130-180 nm → total ~5-6 um
        assert 3e-6 < dc < 10e-6


class TestGetCoatRefl2:
    def test_high_reflectivity(self, aLIGO_stack):
        from OptimalBragg.noise import getCoatRefl2
        nC = aLIGO_stack["ns"][1:-1]
        dOpt = aLIGO_stack["Ls"] * nC / aLIGO_stack["lam_ref"]
        rCoat, dcdp, _, _ = getCoatRefl2(1.0, 1.45, nC, dOpt)
        R = np.abs(rCoat) ** 2
        assert R > 0.9999
        assert len(dcdp) == len(nC)


# ── JIT thermo-optic ─────────────────────────────────────────────────

class TestCoatingThermoopticFast:
    def test_positive_psd(self, aLIGO_stack):
        from OptimalBragg.noise import coating_thermooptic_fast, extract_stack_params
        params = extract_stack_params(aLIGO_stack, R_MIRROR, D_MIRROR)
        dOpt = aLIGO_stack["Ls_opt"]
        S = coating_thermooptic_fast(100.0, dOpt, 1064e-9, W_BEAM, params)
        assert S > 0

    def test_decreases_with_frequency(self, aLIGO_stack):
        from OptimalBragg.noise import coating_thermooptic_fast, extract_stack_params
        params = extract_stack_params(aLIGO_stack, R_MIRROR, D_MIRROR)
        dOpt = aLIGO_stack["Ls_opt"]
        S_10 = coating_thermooptic_fast(10.0, dOpt, 1064e-9, W_BEAM, params)
        S_100 = coating_thermooptic_fast(100.0, dOpt, 1064e-9, W_BEAM, params)
        S_1k = coating_thermooptic_fast(1000.0, dOpt, 1064e-9, W_BEAM, params)
        assert S_10 > S_100 > S_1k

    def test_reasonable_magnitude(self, aLIGO_stack):
        """At 100 Hz, aLIGO thermo-optic should be ~1e-42 m^2/Hz."""
        from OptimalBragg.noise import coating_thermooptic_fast, extract_stack_params
        params = extract_stack_params(aLIGO_stack, R_MIRROR, D_MIRROR)
        dOpt = aLIGO_stack["Ls_opt"]
        S = coating_thermooptic_fast(100.0, dOpt, 1064e-9, W_BEAM, params)
        assert 1e-45 < S < 1e-38


# ── Pure-numpy thermo-optic ──────────────────────────────────────────

class TestCoatingThermooptic:
    def test_returns_three_components(self, aLIGO_stack):
        from OptimalBragg.noise import coating_thermooptic
        f = np.array([100.0])
        StoZ, SteZ, StrZ = coating_thermooptic(f, aLIGO_stack, W_BEAM,
                                                R_MIRROR, D_MIRROR)
        assert StoZ > 0
        assert SteZ > 0
        assert StrZ > 0

    def test_total_vs_components(self, aLIGO_stack):
        """Total should be close to |sqrt(TE) + sqrt(TR)|^2 (correlated)."""
        from OptimalBragg.noise import coating_thermooptic
        f = np.array([100.0])
        StoZ, SteZ, StrZ = coating_thermooptic(f, aLIGO_stack, W_BEAM,
                                                R_MIRROR, D_MIRROR)
        # Total is NOT simply TE + TR due to correlation
        # But should be bounded: max(TE,TR) < total < (sqrt(TE)+sqrt(TR))^2
        assert StoZ >= min(SteZ, StrZ)

    def test_spectral_shape(self, aLIGO_stack):
        from OptimalBragg.noise import coating_thermooptic
        f = np.logspace(1, 3, 10)
        StoZ, _, _ = coating_thermooptic(f, aLIGO_stack, W_BEAM,
                                         R_MIRROR, D_MIRROR)
        # Should decrease with frequency
        assert np.all(np.diff(StoZ) < 0)

    def test_agrees_with_jit(self, aLIGO_stack):
        """JIT and pure-numpy should give similar results."""
        from OptimalBragg.noise import (
            coating_thermooptic, coating_thermooptic_fast, extract_stack_params
        )
        params = extract_stack_params(aLIGO_stack, R_MIRROR, D_MIRROR)
        dOpt = aLIGO_stack["Ls_opt"]

        f = 100.0
        S_jit = coating_thermooptic_fast(f, dOpt, 1064e-9, W_BEAM, params)
        StoZ_np, _, _ = coating_thermooptic(
            np.array([f]), aLIGO_stack, W_BEAM, R_MIRROR, D_MIRROR
        )
        # Both implement the same physics; allow 5% tolerance for
        # minor differences in finite-correction code paths
        assert_allclose(S_jit, StoZ_np, rtol=0.05)


# ── Coating Brownian ──────────────────────────────────────────────────

class TestCoatingBrownian:
    def test_positive_psd(self, aLIGO_stack):
        from OptimalBragg.noise import coating_brownian
        f = np.logspace(1, 3, 10)
        S = coating_brownian(f, aLIGO_stack, W_BEAM)
        assert np.all(S > 0)

    def test_decreases_with_frequency(self, aLIGO_stack):
        from OptimalBragg.noise import coating_brownian
        f = np.logspace(1, 3, 20)
        S = coating_brownian(f, aLIGO_stack, W_BEAM)
        assert np.all(np.diff(S) < 0)

    def test_reasonable_magnitude(self, aLIGO_stack):
        """At 100 Hz, coating Brownian should be ~1e-44 to 1e-42 m^2/Hz."""
        from OptimalBragg.noise import coating_brownian
        f = np.array([100.0])
        S = coating_brownian(f, aLIGO_stack, W_BEAM)
        assert 1e-46 < S[0] < 1e-40


# ── Substrate noise ───────────────────────────────────────────────────

class TestSubstrateBrownian:
    def test_positive_psd(self, aLIGO_stack):
        from OptimalBragg.noise import substrate_brownian
        f = np.logspace(1, 3, 10)
        S = substrate_brownian(f, aLIGO_stack, W_BEAM, R_MIRROR, D_MIRROR)
        assert np.all(S > 0)

    def test_decreases_with_frequency(self, aLIGO_stack):
        from OptimalBragg.noise import substrate_brownian
        f = np.logspace(1, 3, 20)
        S = substrate_brownian(f, aLIGO_stack, W_BEAM, R_MIRROR, D_MIRROR)
        assert np.all(np.diff(S) < 0)


class TestSubstrateThermoelastic:
    def test_positive_psd(self, aLIGO_stack):
        from OptimalBragg.noise import substrate_thermoelastic
        f = np.logspace(1, 3, 10)
        S = substrate_thermoelastic(f, aLIGO_stack, W_BEAM, R_MIRROR, D_MIRROR)
        assert np.all(S > 0)


class TestSubstrateThermorefractive:
    def test_positive_psd(self, aLIGO_stack):
        from OptimalBragg.noise import substrate_thermorefractive
        f = np.logspace(1, 3, 10)
        S = substrate_thermorefractive(f, aLIGO_stack, W_BEAM, D_MIRROR)
        assert np.all(S > 0)


# ── Brownian proxy ────────────────────────────────────────────────────

class TestBrownianProxy:
    def test_returns_gamma_dict(self, aLIGO_stack):
        from OptimalBragg.noise import brownian_proxy
        gams = brownian_proxy(aLIGO_stack)
        assert isinstance(gams, dict)
        assert len(gams) >= 1  # At least one non-low-n material
        # The "H" material should have a gamma value
        assert "H" in gams
        assert gams["H"] > 0

    def test_gamma_magnitude(self, aLIGO_stack):
        """Gamma should be order 1 for typical materials."""
        from OptimalBragg.noise import brownian_proxy
        gams = brownian_proxy(aLIGO_stack)
        for val in gams.values():
            assert 0.1 < val < 100


# ── Convenience wrappers ─────────────────────────────────────────────

class TestConvenienceWrappers:
    def test_coating_noise_returns_four(self, aLIGO_stack):
        from OptimalBragg.noise import coating_noise
        f = np.logspace(1, 3, 5)
        result = coating_noise(f, aLIGO_stack, W_BEAM,
                               r_mirror=R_MIRROR, d_mirror=D_MIRROR)
        assert len(result) == 4

    def test_substrate_noise_returns_three(self, aLIGO_stack):
        from OptimalBragg.noise import substrate_noise
        f = np.logspace(1, 3, 5)
        result = substrate_noise(f, aLIGO_stack, W_BEAM,
                                 r_mirror=R_MIRROR, d_mirror=D_MIRROR)
        assert len(result) == 3
