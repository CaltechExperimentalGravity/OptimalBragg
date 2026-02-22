"""Tests for OptimalBragg.optimizer — coating optimization pipeline."""

import numpy as np
import pytest
from pathlib import Path

from OptimalBragg import Material, qw_stack
from OptimalBragg.materials import SiO2, TiTa2O5, FusedSilica, air


# ── Unit tests for _build_stack ────────────────────────────────────────

class TestBuildStack:
    def test_builds_correct_stack(self):
        from OptimalBragg.optimizer import _build_stack
        from OptimalBragg.io import load_materials_yaml

        mat_file = Path(__file__).parent.parent / "projects/aLIGO/materials.yml"
        if not mat_file.exists():
            pytest.skip("projects/aLIGO/materials.yml not found")

        materials = load_materials_yaml(str(mat_file))
        stack = _build_stack(materials, Npairs=6, optic="ETM")

        assert stack["lam_ref"] == pytest.approx(1064e-9)
        assert len(stack["Ls"]) == 12  # 6 bilayers = 12 layers
        assert len(stack["ns"]) == 14  # 12 layers + vac + sub
        assert stack["w_beam"] == pytest.approx(0.062)

    def test_pattern_correct(self):
        from OptimalBragg.optimizer import _build_stack
        from OptimalBragg.io import load_materials_yaml

        mat_file = Path(__file__).parent.parent / "projects/aLIGO/materials.yml"
        if not mat_file.exists():
            pytest.skip("projects/aLIGO/materials.yml not found")

        materials = load_materials_yaml(str(mat_file))
        stack = _build_stack(materials, Npairs=4)
        assert stack["pattern"] == "LH" * 4


# ── Integration test: short optimization ─────────────────────────────

class TestRunOptimization:
    @pytest.mark.slow
    def test_short_optimization(self, tmp_path):
        """Run a 2-bilayer optimization with minimal settings."""
        from OptimalBragg.optimizer import run_optimization

        # Write a minimal params file
        params = tmp_path / "test_params.yml"
        params.write_text("""\
costs:
    Trans1064:
        target: 0.5
        weight: 5
    Brownian:
        target: 20.0
        weight: 1

misc:
    pol: 'te'
    aoi: 0
    Npairs: 2
    Nfixed: 0
    Ncopies: 0
    Nparticles: 5
    atol: 0.1
    tol: 0.1
    init_method: 'halton'
    lambdaAUX: 0.5
    materials_file: 'materials.yml'
""")

        # Write a minimal materials file
        mat = tmp_path / "materials.yml"
        mat.write_text("""\
substrate:
  material: FusedSilica
  overrides:
    Temp: 295
    MassRadius: 0.17
    MassThickness: 0.20

thin_films:
  L:
    material: SiO2
  H:
    material: TiTa2O5

laser:
  wavelength: 1064e-9
  power: 125

optics:
  ETM:
    beam_radius: 0.062
""")

        result = run_optimization(str(params), save=False, optic="ETM")

        assert "L" in result
        assert "scalar_cost" in result
        assert "output" in result
        assert result["scalar_cost"] > 1.0  # multiplicative, always >= 1
        assert len(result["L"]) == 4  # 2 bilayers = 4 layers

    @pytest.mark.slow
    def test_saves_hdf5(self, tmp_path):
        """Verify HDF5 output is written."""
        from OptimalBragg.optimizer import run_optimization

        params = tmp_path / "ETM_params.yml"
        params.write_text("""\
costs:
    Trans1064:
        target: 0.5
        weight: 5

misc:
    pol: 'te'
    aoi: 0
    Npairs: 2
    Nfixed: 0
    Ncopies: 0
    Nparticles: 5
    atol: 0.1
    tol: 0.1
    init_method: 'halton'
    lambdaAUX: 0.5
    materials_file: 'materials.yml'
""")

        mat = tmp_path / "materials.yml"
        mat.write_text("""\
substrate:
  material: FusedSilica
  overrides:
    Temp: 295
    MassRadius: 0.17
    MassThickness: 0.20

thin_films:
  L:
    material: SiO2
  H:
    material: TiTa2O5

laser:
  wavelength: 1064e-9

optics:
  ETM:
    beam_radius: 0.062
""")

        result = run_optimization(str(params), save=True, optic="ETM")

        # Check HDF5 was created
        data_dir = tmp_path / "Data" / "ETM"
        assert data_dir.exists()
        hdf5_files = list(data_dir.glob("*.hdf5"))
        assert len(hdf5_files) == 1
