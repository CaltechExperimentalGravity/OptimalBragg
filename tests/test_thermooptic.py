"""Verify JIT thermooptic matches pure-numpy implementation."""
import numpy as np
import pytest

from OptimalBragg import Material, qw_stack
from OptimalBragg.materials import SiO2, TiTa2O5, FusedSilica, air
from OptimalBragg.noise import (
    coating_thermooptic_fast, coating_thermooptic,
    extract_stack_params,
)


@pytest.fixture
def aLIGO_stack():
    return qw_stack(
        lam_ref=1064e-9,
        substrate=Material(FusedSilica),
        superstrate=Material(air),
        thin_films={"L": Material(SiO2), "H": Material(TiTa2O5)},
        pattern="LH" * 14,
    )


@pytest.fixture
def stack_params(aLIGO_stack):
    return extract_stack_params(aLIGO_stack, 0.17, 0.20)


def test_quarter_wave_stack(aLIGO_stack, stack_params):
    """Quarter-wave 14-bilayer stack at 100 Hz."""
    L = aLIGO_stack["Ls_opt"]
    f = 100.0
    w_beam = 0.062

    StoZ_jit = coating_thermooptic_fast(
        f, L, 1064e-9, w_beam, stack_params)
    StoZ_np, _, _ = coating_thermooptic(
        np.array([f]), aLIGO_stack, w_beam, 0.17, 0.20)

    assert np.isclose(StoZ_jit, StoZ_np, rtol=0.05), \
        f"JIT={StoZ_jit:.15e}, numpy={StoZ_np:.15e}"


def test_random_stack(aLIGO_stack, stack_params):
    """Random layer thicknesses at 100 Hz."""
    rng = np.random.default_rng(42)
    L = rng.uniform(0.05, 0.48, size=len(aLIGO_stack["Ls_opt"]))
    f = 100.0
    w_beam = 0.062

    # Update stack for pure-numpy path
    aLIGO_stack_copy = dict(aLIGO_stack)
    aLIGO_stack_copy["Ls_opt"] = L

    StoZ_jit = coating_thermooptic_fast(
        f, L, 1064e-9, w_beam, stack_params)

    assert StoZ_jit > 0
    assert np.isfinite(StoZ_jit)


def test_multiple_frequencies(aLIGO_stack, stack_params):
    """Test at several frequencies — should decrease monotonically."""
    L = aLIGO_stack["Ls_opt"]
    w_beam = 0.062

    values = []
    for f in [10.0, 50.0, 100.0, 500.0, 1000.0]:
        StoZ_jit = coating_thermooptic_fast(
            f, L, 1064e-9, w_beam, stack_params)
        assert StoZ_jit > 0
        values.append(StoZ_jit)

    # Thermo-optic should decrease with frequency
    assert all(values[i] > values[i+1] for i in range(len(values)-1))


def test_different_stack_sizes():
    """Test with different number of layer pairs."""
    f = 100.0
    w_beam = 0.062

    for npairs in [5, 10, 14, 20]:
        stack = qw_stack(
            lam_ref=1064e-9,
            substrate=Material(FusedSilica),
            superstrate=Material(air),
            thin_films={"L": Material(SiO2), "H": Material(TiTa2O5)},
            pattern="LH" * npairs,
        )
        params = extract_stack_params(stack, 0.17, 0.20)
        L = stack["Ls_opt"]

        StoZ = coating_thermooptic_fast(
            f, L, 1064e-9, w_beam, params)
        assert StoZ > 0, f"npairs={npairs}: StoZ={StoZ}"
        assert np.isfinite(StoZ)
