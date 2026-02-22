"""Tests for OptimalBragg.materials and Material class."""
import numpy as np
import pytest
from OptimalBragg import Material
from OptimalBragg.materials import SiO2, TiTa2O5, Ta2O5, aSi_123, SiN_123, cSi_123, air, FusedSilica


class TestMaterial:
    def test_creates_from_dict(self):
        m = Material(SiO2)
        assert m.Name == "SiO2"
        assert m.Index == 1.45
        assert m.Y == 60e9

    def test_all_coating_materials_have_required_props(self):
        """Every coating material needs these for thermal noise calculations."""
        required = {'Name', 'Index', 'Y', 'Sigma', 'Phi', 'Alpha', 'Beta', 'CV', 'ThermalDiffusivity'}
        for name, mat in [('SiO2', SiO2), ('TiTa2O5', TiTa2O5), ('aSi_123', aSi_123), ('SiN_123', SiN_123)]:
            m = Material(mat)
            for prop in required:
                assert hasattr(m, prop), f"{name} missing {prop}"

    def test_substrate_materials_have_extra_props(self):
        """Substrate materials need additional properties for noise calculations."""
        extra = {'MassDensity', 'MassAlpha', 'MassCM', 'MassKappa', 'MirrorY', 'MirrorSigma'}
        for name, mat in [('FusedSilica', FusedSilica), ('cSi_123', cSi_123)]:
            m = Material(mat)
            for prop in extra:
                assert hasattr(m, prop), f"Substrate {name} missing {prop}"

    def test_air_minimal(self):
        m = Material(air)
        assert m.Index == 1.0


class TestQwStack:
    def test_basic_stack(self):
        from OptimalBragg import qw_stack
        stack = qw_stack(
            lam_ref=1064e-9,
            substrate=Material(SiO2),
            superstrate=Material(air),
            thin_films={"L": Material(SiO2), "H": Material(TiTa2O5)},
            pattern="LH" * 5,
        )
        assert 'Ls' in stack
        assert 'Ls_opt' in stack
        assert 'ns' in stack
        assert len(stack['Ls']) == 10
        assert len(stack['ns']) == 12  # sup + 10 layers + sub

    def test_qw_optical_thickness_is_quarter(self):
        from OptimalBragg import qw_stack
        stack = qw_stack(
            lam_ref=1064e-9,
            substrate=Material(SiO2),
            superstrate=Material(air),
            thin_films={"L": Material(SiO2), "H": Material(TiTa2O5)},
            pattern="LH" * 5,
        )
        # For QW stack: all Ls_opt should be 0.25
        assert np.allclose(stack['Ls_opt'], 0.25, atol=1e-15)

    def test_hwcap_optical_thickness(self):
        from OptimalBragg import qw_stack
        stack = qw_stack(
            lam_ref=1064e-9,
            substrate=Material(SiO2),
            superstrate=Material(air),
            thin_films={"L": Material(SiO2), "H": Material(TiTa2O5)},
            pattern="LH" * 3,
            hwcap="H",
        )
        # First layer (hwcap) should be 0.5, rest 0.25
        assert abs(stack['Ls_opt'][0] - 0.5) < 1e-15
        assert np.allclose(stack['Ls_opt'][1:], 0.25, atol=1e-15)

    def test_physical_optical_consistency(self):
        from OptimalBragg import qw_stack
        stack = qw_stack(
            lam_ref=1064e-9,
            substrate=Material(SiO2),
            superstrate=Material(air),
            thin_films={"L": Material(SiO2), "H": Material(TiTa2O5)},
            pattern="LH" * 5,
        )
        ns_layers = stack['ns'][1:-1]
        expected_phys = stack['lam_ref'] * stack['Ls_opt'] / ns_layers
        assert np.allclose(stack['Ls'], expected_phys)

    def test_material_property_arrays(self):
        from OptimalBragg import qw_stack
        stack = qw_stack(
            lam_ref=1064e-9,
            substrate=Material(SiO2),
            superstrate=Material(air),
            thin_films={"L": Material(SiO2), "H": Material(TiTa2O5)},
            pattern="LH" * 5,
        )
        n_layers = 10
        for key in ('alphas', 'Ys', 'sigmas', 'phis', 'ctes', 'betas', 'Cvs', 'thermaldiffs'):
            assert len(stack[key]) == n_layers, f"{key} has wrong length"
