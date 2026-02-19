# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Optical coating optimization for gravitational wave detector mirrors (LIGO Voyager, 40m prototype). Designs high-reflectivity dielectric stacks by globally optimizing a multi-objective cost function that balances transmissivity, thermal noise, manufacturing tolerance, and surface E-field constraints.

## Environment Setup

```bash
conda env create -f environment.yml
conda activate coatingDev
pip install -e .
```

Key dependencies: numpy, scipy, matplotlib, numba (>=0.56), emcee, corner, gwinc (>=0.6), lmfit, h5py, pytest. The project uses Python >=3.10 with conda-forge.

MATLAB engine for Python is needed for legacy MATLAB workflows (`$MATLABROOT/extern/engines/python`). Set `GWINCPATH` environment variable to point to MATLAB GWINC installation.

## Common Workflows

**Run coating optimization** (from the repo root):
```bash
cd SiN_aSi && python mkETM.py          # Optimize End Test Mass
cd SiN_aSi && python mkITM.py          # Optimize Input Test Mass
```
Output: `Data/ETM/ETM_Layers_YYMMDD_HHMMSS.hdf5`

**Monte Carlo sensitivity analysis:**
```bash
python doMC.py Data/ETM/ETM_Layers_YYMMDD_HHMMSS.hdf5 output_MC.hdf5 5000
```

**Visualization:**
```bash
python plot_ETM.py       # Design analysis dashboard
python cornerPlt.py      # MC corner plots from .hdf5
python plotlayers.py     # Layer structure + E-field
```

**Run tests:**
```bash
pytest tests/ -v -k "not slow"           # Unit tests (~1 sec)
pytest tests/test_integration.py -m slow  # Integration test (~30 sec)
```

**Benchmarks:**
```bash
python benchmarks/bench_multidiel1.py     # multidiel1 micro-benchmark
```

## Architecture

### Core Libraries (`generic/`)

Contains the shared physics and optimization code. All projects import directly from `generic/` — there are no project-local copies.

- **`coatingUtils.py`** — Transfer matrix method (`multidiel1` with Numba JIT), optical-to-physical thickness conversion, E-field depth profiles, Sellmeier dispersion, spectral reflectivity, absorption calculation, YAML parameter import
- **`optimUtils.py`** — Cost function components: transmission, Brownian noise proxy, thermooptic noise, layer sensitivity, surface E-field, absorption. Master evaluator: `getMirrorCost(L, costs, ifo, gam, verbose, misc)`. Uses consolidated `multidiel1` calls (2 per evaluation instead of 6).

### Configuration (YAML)

Each project directory has:
- **`ETM_params.yml`, `ITM_params.yml`** — Cost function weights/targets, optimizer settings (population size, tolerance, number of layer pairs), `lambdaAUX` wavelength ratio
- **GWINC structure file** (e.g. `aSiSiN.yaml`) — Material properties (elastic, thermal, optical), substrate, laser specs, mirror geometry

### Optimization

Python uses `scipy.optimize.differential_evolution` with multiprocessing (`workers=-1`). Layer thicknesses are optimized as optical thicknesses bounded to [0.05, 0.48]. Legacy MATLAB uses Particle Swarm Optimization (`runSwarm.m`).

### Monte Carlo

Uses `emcee` ensemble sampler (20 walkers, 4D parameter space). Perturbs: angle of incidence, high-n index, low-n index, layer thickness — all as 0.5% Gaussian. Reads HDF5 optimizer output. Outputs to HDF5.

### Data Flow

`ETM_params.yml` → `mkETM.py` → optimizer → `.hdf5` → `plot_ETM.py` / `doMC.py` → `.hdf5` → `cornerPlt.py`

## Key Conventions

- All projects import from `generic/` directly — no local copies
- Binary data files (`.mat`, `.hdf5`) are gitignored; `.mat` files use Git LFS
- Plotting uses matplotlib style `'gvELOG'` with `'bmh'` fallback
- Physical units: wavelengths in nm, layer thicknesses as optical thickness (fraction of lambda_0), thermal noise in m/sqrt(Hz), temperature in Kelvin
- Primary operating point: 2050 nm wavelength, 123 K (Voyager), with auxiliary wavelengths at 1550 nm
- AUX wavelength ratio is configured per-project via `lambdaAUX` in params YAML (not hardcoded)
