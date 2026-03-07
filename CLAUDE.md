# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Behavioral Rules

These rules override all defaults. They exist because of real failures. Do not rationalize around them.

1. **When the user says "look at X" — read X immediately.** Do not regenerate, recalculate, or re-derive. Open the file they pointed to and read it. The user is telling you where the evidence is.

2. **Verify outputs before claiming success.** After any pipeline step (optimize, plot, MC, corner, report), examine the actual outputs: render plots to PNG and look at them, read the .rst report, check for NaN/missing traces/broken labels. Never say "done" based on exit code alone.

3. **Do not argue with the user about what they observe.** If the user says a plot is broken or a value is wrong, it IS broken. Read what they point you to. Do not run new calculations to "prove" it's fine.

4. **"This repo" means `/Users/rana/Desktop/Dropbox/GIT/40m/Coatings`.** Not `~/.claude/`, not `/tmp/`. Instruction files live here: `CLAUDE.md`, `CODEMAP.md`, `README.md`.

5. **Update THIS FILE when you learn a new rule.** Context compaction will erase conversation history. The only durable memory is in repo files. When a mistake reveals a new pattern, add it to the relevant section of this CLAUDE.md so future sessions inherit the lesson.

## Critical Units Convention

- `stack['Ls']` = physical layer thickness in **metres** (SI). Set by `qw_stack()` as `(lam_ref / 4) / n`.
- `stack['Ls_opt']` = optical thickness as **fraction of lambda** (dimensionless, ~0.25 for QW).
- `op2phys(L_opt, n)` returns `L_opt / n` — still a **fraction**, NOT metres. Multiply by `wavelength` to get metres.
- Always check units when writing code that touches `stack['Ls']`. A 1e6 error here silently corrupts thermal noise.

## Project Overview

Optical coating optimization for gravitational wave detector mirrors. The `OptimalBragg` package designs high-reflectivity dielectric stacks by globally optimizing a multiplicative cost function that balances transmissivity, thermal noise, manufacturing tolerance, and surface E-field constraints.

Four project configurations in `projects/`:
- **`projects/aLIGO/`** — aLIGO: SiO2/TiTa2O5 coatings at 1064 nm, 295 K, fused silica substrate
- **`projects/SFG/`** — 40m SFG cavity: SiO2/Ta2O5 coatings at 1064 nm, 295 K, fused silica substrate
- **`projects/Voyager_aSiSiN/`** — LIGO Voyager: aSi/SiN coatings at 2050 nm, 123 K, silicon substrate
- **`projects/Voyager_Ta2O5/`** — LIGO Voyager: Ta2O5/SiO2 coatings at 2050 nm, 123 K, silicon substrate

## Environment Setup

```bash
conda env create -f environment.yml
conda activate coatingDev
pip install -e .
```

Key dependencies: numpy, scipy, matplotlib, numba (>=0.56), emcee, arviz, h5py, pytest. The project uses Python >=3.10 with conda-forge.

## Common Workflows

**Run coating optimization** (CLI):
```bash
optimalbragg optimize projects/aLIGO/ETM_params.yml
optimalbragg optimize projects/Voyager_aSiSiN/ETM_params.yml
```
Output: `Data/ETM/ETM_Layers_YYMMDD_HHMMSS.hdf5`

**Or from a project directory:**
```bash
cd projects/aLIGO && python mkETM.py
```

**Sweep Npairs to find optimal layer count:**
```bash
optimalbragg sweep projects/aLIGO/ETM_params.yml
optimalbragg sweep projects/aLIGO/ETM_params.yml --min-pairs 16 --max-pairs 24
```

**Monte Carlo sensitivity analysis:**
```bash
optimalbragg mc Data/ETM/ETM_Layers_YYMMDD_HHMMSS.hdf5 5000
```

**Visualization:**
```bash
optimalbragg plot Data/ETM/ETM_Layers_YYMMDD_HHMMSS.hdf5
optimalbragg corner Data/ETM/ETM_MC.hdf5
```

**Run tests:**
```bash
pytest tests/ -v -k "not slow"           # Unit tests (~2 sec)
pytest tests/test_integration.py -m slow  # Integration test (~30 sec)
```

**Benchmarks:**
```bash
python benchmarks/bench_multidiel1.py     # Transfer matrix benchmark
python benchmarks/bench_thermooptic.py    # Thermooptic JIT vs numpy
```

## Architecture

### Core Package (`OptimalBragg/`)

- **`layers.py`** — Transfer matrix method (`multidiel1` with Numba JIT), `multilayer_diel()`, `field_zmag()`, `surfield()`, `calc_abs()`, `op2phys()`, Sellmeier dispersion, spectral reflectivity/transmission
- **`costs.py`** — Cost function components: `transmissionCost`, `brownianCost`, `thermoopticCost`, `sensitivityCost`, `surfEfieldCost`, `absorptionCost`. Master evaluator: `getMirrorCost(L, costs, stack, gam, verbose, misc)`. Multiplicative cost: `C = prod(1 + w_i * c_i)`
- **`noise.py`** — Thermal noise models: `coating_thermooptic_fast()` (Numba JIT), `coating_brownian()` (Hong et al.), `substrate_brownian()`, `substrate_thermoelastic()`, `brownian_proxy()`. All take `stack` dict.
- **`materials.py`** — Central materials library with references: `SiO2`, `TiTa2O5`, `Ta2O5`, `aSi_123`, `SiN_123`, `cSi_123`, `FusedSilica`, `air`
- **`optimizer.py`** — `run_optimization(params_yaml)` using `differential_evolution`
- **`plot.py`** — `plot_layers()`, `plot_spectral()`, `plot_noise()`, `plot_starfish()`, `plot_corner()`
- **`mc.py`** — `run_mc()` — Monte Carlo sensitivity via `emcee` (3D: high-n, low-n, thickness perturbations)
- **`report.py`** — Sphinx report generation: `generate_run_rst()`, `build_html()`
- **`io.py`** — `h5read()`, `h5write()`, `yamlread()`, `load_materials_yaml()`
- **`cli.py`** — CLI entry point: `optimalbragg {optimize,plot,mc,corner}`
- **`__init__.py`** — Public API: `Material` class, `qw_stack()`, `load_materials_yaml()`

### Configuration (YAML)

Each project directory has:
- **`materials.yml`** — Material properties referencing the central library, with optional overrides
- **`ETM_params.yml`, `ITM_params.yml`** — Cost function weights/targets, optimizer settings, `lambda2`/`lambda3` wavelength ratios

### Data Flow

`ETM_params.yml` + `materials.yml` → `optimalbragg optimize` → `.hdf5` → `optimalbragg plot` / `optimalbragg mc` → `.hdf5` → `optimalbragg corner`

### Monte Carlo

Uses `emcee` ensemble sampler (20 walkers, 3D parameter space). Perturbs: high-n index, low-n index, layer thickness — all as 0.5% Gaussian. Reads HDF5 optimizer output. Outputs to HDF5.

## Key Conventions

- All physics/optimization code lives in `OptimalBragg/` — project directories contain only YAML configs and thin wrapper scripts
- The `stack` dict is the central data structure — contains `ns`, `Ls`, `Ls_opt`, material property arrays, `sub`/`sup` Material objects
- Binary data files (`.mat`, `.hdf5`) are gitignored; `.mat` files use Git LFS
- Physical units: wavelengths in meters (SI), layer thicknesses as optical thickness (fraction of lambda_0), thermal noise in m^2/Hz (PSD), temperature in Kelvin
- Primary operating points: 2050 nm / 123 K (Voyager), 1064 nm / 295 K (aLIGO)
- Wavelength ratios are configured per-project via `lambda2`/`lambda3` in params YAML (not hardcoded)
- Corner plots use ArviZ `plot_pair()` (not `corner` package)
