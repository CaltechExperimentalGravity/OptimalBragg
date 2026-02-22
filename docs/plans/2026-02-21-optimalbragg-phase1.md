# Phase 1: OptimalBragg Package Skeleton + Materials Library

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Create the `OptimalBragg/` Python package with Material class, qw_stack(), I/O utilities, materials library, and materials YAML configs for all three projects.

**Architecture:** Single flat package `OptimalBragg/` with `__init__.py` (Material, qw_stack, public API), `materials.py` (library dicts), `io.py` (HDF5/YAML helpers). Projects get `materials.yml` files that reference the library. Stack dict is the central data structure (both `Ls` physical and `Ls_opt` optical).

**Tech Stack:** Python 3.10+, numpy, h5py, pyyaml

**Source references:**
- OptimalBragg original: `/tmp/OptimalBragg/OptimalBragg/`
- Our codebase: `/Users/rana/Desktop/Dropbox/GIT/40m/Coatings/`
- Material params: `Arms/aLIGO_SiO2Ta2O5.yaml`, `SiN_aSi/aSiSiN.yaml`, `/tmp/OptimalBragg/OptimalBragg/materials.py`

---

### Task 1: Create `OptimalBragg/materials.py`

**Files:**
- Create: `OptimalBragg/materials.py`

**Step 1: Write materials.py**

Copy from `/tmp/OptimalBragg/OptimalBragg/materials.py` verbatim, then add missing properties needed by our thermal noise code. Each material dict needs `Properties` + `References` keys. Add any properties present in our gwinc YAMLs but missing from OptimalBragg dicts:

For room-temp materials, cross-reference `Arms/aLIGO_SiO2Ta2O5.yaml`:
- `SiO2`: add `Temp: 295` (already has it as 300 — update to 295)
- `TiTa2O5`: add `rho`, `CP`, `Temp: 295`
- `Ta2O5`: already complete

For cryo materials, cross-reference `SiN_aSi/aSiSiN.yaml`:
- `aSi_123`: verify `CV: 372` (gwinc) vs `CV: 1.05e6` (OptimalBragg) — gwinc's is per-mass (J/kg/K), OptimalBragg's is per-volume (J/m³/K). Keep OptimalBragg convention (per-volume).
- `SiN_123`: verify `CV: 230` (gwinc) vs OptimalBragg's 230 — gwinc is per-mass. Need per-volume for our code.

Key: all CVs must be in J/m³/K (volumetric) since `thermoopticUtils.py` expects that. Our gwinc YAMLs use per-mass CVs multiplied by density internally.

**Step 2: Commit**

```bash
git add OptimalBragg/materials.py
git commit -m "feat(OptimalBragg): add materials library with references"
```

---

### Task 2: Create `OptimalBragg/io.py`

**Files:**
- Create: `OptimalBragg/io.py`

**Step 1: Write io.py**

Port from `/tmp/OptimalBragg/OptimalBragg/__init__.py` (h5read, h5write, yamlread) plus our `generic/coatingUtils.py` (importParams). Functions:

- `h5read(fname, group, targets)` — from OptimalBragg
- `h5write(fname, h5_dict)` — from OptimalBragg (handles Material objects)
- `yamlread(fname)` — from OptimalBragg (yaml.safe_load wrapper)
- `load_materials_yaml(fname)` — NEW: reads a project materials.yml, resolves material names against the library, applies overrides, returns substrate/thin_films/laser/optics dicts with Material objects

`load_materials_yaml` is the key new function. It replaces `gwinc.Struct.from_file()`:
```python
def load_materials_yaml(fname):
    """Load a project materials.yml and resolve against materials library.

    Returns dict with keys: substrate, thin_films, laser, optics
    """
    raw = yamlread(fname)
    # Look up each material name in OptimalBragg.materials module
    # Apply any overrides from the YAML
    # Return resolved dicts with Material objects
```

**Step 2: Commit**

```bash
git add OptimalBragg/io.py
git commit -m "feat(OptimalBragg): add HDF5/YAML I/O utilities"
```

---

### Task 3: Create `OptimalBragg/__init__.py` with Material class and qw_stack()

**Files:**
- Create: `OptimalBragg/__init__.py`

**Step 1: Write __init__.py**

Port Material class and qw_stack() from `/tmp/OptimalBragg/OptimalBragg/__init__.py`. Key changes from original:

1. `Material` class: same as OptimalBragg (`self.__dict__.update(props["Properties"])`)
2. `qw_stack()`: same signature as OptimalBragg, but also compute and store `Ls_opt` (optical thicknesses as fraction of lambda). The plan says both conventions in the stack dict.

```python
stack = {
    "lam_ref": lam_ref,
    "ns": np.array(n_stack),
    "Ls": np.array(L_stack),          # physical [m]
    "Ls_opt": np.array(L_opt_stack),   # optical [frac lambda]
    "alphas": ..., "Ys": ..., "sigmas": ..., "phis": ...,
    "ctes": ..., "betas": ..., "Cvs": ..., "thermaldiffs": ...,
    "sub": substrate, "sup": superstrate,
    "thin_films": thin_films, "pattern": pattern,
    "optimized": False,
}
```

Public API exports: `Material`, `qw_stack`, `h5read`, `h5write`, `yamlread`, `load_materials_yaml`

**Step 2: Commit**

```bash
git add OptimalBragg/__init__.py
git commit -m "feat(OptimalBragg): add Material class, qw_stack, public API"
```

---

### Task 4: Create project materials.yml files

**Files:**
- Create: `projects/aLIGO/materials.yml`
- Create: `projects/Voyager_aSiSiN/materials.yml`
- Create: `projects/Voyager_Ta2O5/materials.yml`

**Step 1: Write aLIGO materials.yml**

Translate `Arms/aLIGO_SiO2Ta2O5.yaml` (gwinc format) to the new format:

```yaml
substrate:
  material: SiO2
  overrides:
    Temp: 295
    MassRadius: 0.17
    MassThickness: 0.20
    MassDensity: 2200
    MirrorY: 72.7e9
    MirrorSigma: 0.167
    MassAlpha: 3.9e-7
    MassCM: 739
    MassKappa: 1.38
    MechanicalLossExponent: 0.77
    Alphas: 5.2e-12
    c2: 7.6e-12

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
  ITM:
    beam_radius: 0.055
```

**Step 2: Write Voyager_aSiSiN materials.yml**

Translate `SiN_aSi/aSiSiN.yaml`:

```yaml
substrate:
  material: cSi_123
  overrides:
    Temp: 123
    MassRadius: 0.215
    MassThickness: 0.55
    # ... etc from aSiSiN.yaml

thin_films:
  L:
    material: SiN_123
  H:
    material: aSi_123

laser:
  wavelength: 2050.15e-9
  power: 125

optics:
  ETM:
    beam_radius: 0.06
  ITM:
    beam_radius: 0.059
```

**Step 3: Write Voyager_Ta2O5 materials.yml**

Similar pattern using Ta2O5_123 + SiO2_123 materials.

**Step 4: Commit**

```bash
git add projects/
git commit -m "feat: add materials.yml configs for all three projects"
```

---

### Task 5: Write tests for Material, qw_stack, and materials loading

**Files:**
- Create: `tests/test_materials.py`
- Create: `tests/test_io.py`

**Step 1: Write test_materials.py**

```python
"""Tests for OptimalBragg.materials and Material class."""
import numpy as np
import pytest
from OptimalBragg import Material
from OptimalBragg.materials import SiO2, TiTa2O5, Ta2O5, aSi_123, SiN_123, cSi_123, air

class TestMaterial:
    def test_creates_from_dict(self):
        m = Material(SiO2)
        assert m.Name == "SiO2"
        assert m.Index == 1.45
        assert m.Y == 60e9

    def test_all_materials_have_required_props(self):
        required = {'Name', 'Index', 'Y', 'Sigma', 'Phi', 'Alpha', 'Beta', 'CV', 'ThermalDiffusivity'}
        for name, mat in [('SiO2', SiO2), ('TiTa2O5', TiTa2O5), ('aSi_123', aSi_123), ('SiN_123', SiN_123)]:
            m = Material(mat)
            for prop in required:
                assert hasattr(m, prop), f"{name} missing {prop}"

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

    def test_optical_physical_consistency(self):
        from OptimalBragg import qw_stack
        stack = qw_stack(
            lam_ref=1064e-9,
            substrate=Material(SiO2),
            superstrate=Material(air),
            thin_films={"L": Material(SiO2), "H": Material(TiTa2O5)},
            pattern="LH" * 5,
        )
        # For QW stack: Ls_opt should all be 0.25
        assert np.allclose(stack['Ls_opt'], 0.25)
        # Physical = lam_ref * Ls_opt / n
        ns_layers = stack['ns'][1:-1]
        expected_phys = stack['lam_ref'] * stack['Ls_opt'] / ns_layers
        assert np.allclose(stack['Ls'], expected_phys)
```

**Step 2: Write test_io.py**

```python
"""Tests for OptimalBragg.io."""
import numpy as np
import pytest
import tempfile, os
from OptimalBragg.io import h5read, h5write, yamlread

class TestH5ReadWrite:
    def test_roundtrip(self, tmp_path):
        fname = str(tmp_path / 'test.hdf5')
        data = {'Ls': np.array([1.0, 2.0, 3.0]), 'lam_ref': 1064e-9}
        h5write(fname, data)
        result = h5read(fname, '/', ['Ls'])
        assert np.allclose(result['Ls'], data['Ls'])

class TestYamlread:
    def test_reads_yaml(self, tmp_path):
        fname = str(tmp_path / 'test.yml')
        with open(fname, 'w') as f:
            f.write('key: value\nnested:\n  a: 1\n')
        result = yamlread(fname)
        assert result['key'] == 'value'
        assert result['nested']['a'] == 1
```

**Step 3: Run tests**

```bash
source ~/.zshrc && conda activate coatingDev && pytest tests/test_materials.py tests/test_io.py -v
```

Expected: all PASS

**Step 4: Commit**

```bash
git add tests/test_materials.py tests/test_io.py
git commit -m "test: add tests for Material, qw_stack, and I/O"
```

---

### Task 6: Update pyproject.toml

**Files:**
- Modify: `pyproject.toml`

**Step 1: Update package config**

Change `[tool.setuptools.packages.find]` to include OptimalBragg:
```toml
[tool.setuptools.packages.find]
include = ["OptimalBragg*", "generic*"]
```

Keep `generic*` for now (removed in Phase 8).

**Step 2: Commit**

```bash
git add pyproject.toml
git commit -m "build: include OptimalBragg in package find"
```

---

## Verification

After all tasks:
```bash
source ~/.zshrc && conda activate coatingDev && pytest tests/test_materials.py tests/test_io.py -v
```

All tests pass. The `OptimalBragg` package is importable with `Material`, `qw_stack`, `h5read`, `h5write`, `yamlread`, `load_materials_yaml`.
