"""Verify JIT thermooptic matches gwinc to machine precision."""
import copy
import numpy as np
import pytest
import gwinc
import gwinc.noise.coatingthermal

from generic.thermoopticUtils import coating_thermooptic_fast, extract_ifo_params


@pytest.fixture
def ifo():
    return gwinc.Struct.from_file('SiN_aSi/aSiSiN.yaml')


@pytest.fixture
def ifo_params(ifo):
    return extract_ifo_params(ifo)


def _gwinc_thermooptic(f, L, ifo):
    """Call gwinc's coating_thermooptic the same way thermoopticCost did."""
    mir = copy.copy(ifo.Optics.ETM)
    mir.Coating = copy.copy(ifo.Optics.ETM.Coating)
    mir.Coating.dOpt = L
    StoZ, _, _, _ = gwinc.noise.coatingthermal.coating_thermooptic(
        f, mir, ifo.Laser.Wavelength, ifo.Optics.ETM.BeamRadius)
    return StoZ


def test_quarter_wave_stack(ifo, ifo_params):
    """Quarter-wave 14-bilayer stack at 100 Hz."""
    L = 0.25 * np.ones(29)
    f = 100.0

    StoZ_gwinc = _gwinc_thermooptic(f, L, ifo)
    StoZ_jit = coating_thermooptic_fast(
        f, L, ifo.Laser.Wavelength, ifo.Optics.ETM.BeamRadius, ifo_params)

    assert np.isclose(StoZ_jit, StoZ_gwinc, rtol=1e-12), \
        f"JIT={StoZ_jit:.15e}, gwinc={StoZ_gwinc:.15e}"


def test_random_stack(ifo, ifo_params):
    """Random layer thicknesses at 100 Hz."""
    rng = np.random.default_rng(42)
    L = rng.uniform(0.05, 0.48, size=29)
    f = 100.0

    StoZ_gwinc = _gwinc_thermooptic(f, L, ifo)
    StoZ_jit = coating_thermooptic_fast(
        f, L, ifo.Laser.Wavelength, ifo.Optics.ETM.BeamRadius, ifo_params)

    assert np.isclose(StoZ_jit, StoZ_gwinc, rtol=1e-12), \
        f"JIT={StoZ_jit:.15e}, gwinc={StoZ_gwinc:.15e}"


def test_multiple_frequencies(ifo, ifo_params):
    """Test at several frequencies."""
    L = 0.25 * np.ones(29)

    for f in [10.0, 50.0, 100.0, 500.0, 1000.0]:
        StoZ_gwinc = _gwinc_thermooptic(f, L, ifo)
        StoZ_jit = coating_thermooptic_fast(
            f, L, ifo.Laser.Wavelength, ifo.Optics.ETM.BeamRadius, ifo_params)
        assert np.isclose(StoZ_jit, StoZ_gwinc, rtol=1e-12), \
            f"f={f}: JIT={StoZ_jit:.15e}, gwinc={StoZ_gwinc:.15e}"


def test_different_stack_sizes(ifo, ifo_params):
    """Test with different number of layer pairs."""
    f = 100.0
    for npairs in [5, 10, 14, 20]:
        nlayers = 2 * npairs + 1
        L = 0.25 * np.ones(nlayers)
        StoZ_gwinc = _gwinc_thermooptic(f, L, ifo)
        StoZ_jit = coating_thermooptic_fast(
            f, L, ifo.Laser.Wavelength, ifo.Optics.ETM.BeamRadius, ifo_params)
        assert np.isclose(StoZ_jit, StoZ_gwinc, rtol=1e-12), \
            f"npairs={npairs}: JIT={StoZ_jit:.15e}, gwinc={StoZ_gwinc:.15e}"
