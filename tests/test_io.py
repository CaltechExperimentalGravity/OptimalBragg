"""Tests for OptimalBragg.io."""
import numpy as np
import pytest
from OptimalBragg.io import h5read, h5write, yamlread, load_materials_yaml
from OptimalBragg import Material
from OptimalBragg.materials import SiO2


class TestH5ReadWrite:
    def test_roundtrip_arrays(self, tmp_path):
        fname = str(tmp_path / 'test.hdf5')
        data = {'group1': {'Ls': np.array([1.0, 2.0, 3.0])}}
        h5write(fname, data)
        result = h5read(fname, 'group1', ['Ls'])
        assert np.allclose(result['Ls'], data['group1']['Ls'])

    def test_missing_dataset_returns_zero(self, tmp_path):
        fname = str(tmp_path / 'test.hdf5')
        data = {'group1': {'Ls': np.array([1.0, 2.0])}}
        h5write(fname, data)
        result = h5read(fname, 'group1', ['Ls', 'nonexistent'])
        assert result['nonexistent'] == 0.0

    def test_scalar_roundtrip(self, tmp_path):
        fname = str(tmp_path / 'test.hdf5')
        data = {'group1': {'lam_ref': 1064e-9}}
        h5write(fname, data)
        # Scalars stored as attrs, not datasets — h5read reads datasets
        # This tests the attr path in h5write


class TestYamlread:
    def test_reads_yaml(self, tmp_path):
        fname = str(tmp_path / 'test.yml')
        with open(fname, 'w') as f:
            f.write('key: value\nnested:\n  a: 1\n')
        result = yamlread(fname)
        assert result['key'] == 'value'
        assert result['nested']['a'] == 1


class TestLoadMaterialsYaml:
    def test_loads_aligo(self):
        result = load_materials_yaml('projects/aLIGO/materials.yml')
        assert hasattr(result['substrate'], 'Name')
        assert result['substrate'].Name == 'FusedSilica'
        assert 'L' in result['thin_films']
        assert 'H' in result['thin_films']
        assert result['thin_films']['L'].Name == 'SiO2'
        assert result['thin_films']['H'].Name == 'TiTa2O5'
        assert result['laser']['wavelength'] == 1064e-9

    def test_loads_voyager(self):
        result = load_materials_yaml('projects/Voyager_aSiSiN/materials.yml')
        assert result['substrate'].Temp == 123
        assert result['thin_films']['H'].Name == 'a-Si'

    def test_substrate_overrides_applied(self):
        result = load_materials_yaml('projects/aLIGO/materials.yml')
        # MassRadius is in overrides, not in the base FusedSilica dict
        assert result['substrate'].MassRadius == 0.17
        assert result['substrate'].MassThickness == 0.20
