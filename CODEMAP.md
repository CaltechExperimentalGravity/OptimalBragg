# Code Map: Coating Optimization for Gravitational Wave Detectors

Comprehensive architecture reference for the optical coating optimization codebase.
For quick-start instructions, see [README.md](README.md).

---

## Architecture Overview

![OptimalBragg codemap](docs/codemap.svg)

---

## 1. Repository Layout

### Target structure (in progress)

```
Coatings/
|-- OptimalBragg/                    # Python package (core library)
|   |-- __init__.py                  #   Material class, qw_stack(), public API
|   |-- materials.py                 #   Central materials library with references
|   |-- io.py                        #   HDF5/YAML I/O, load_materials_yaml()
|   |-- layers.py                    #   Transfer matrix (Numba JIT), E-field, spectral
|   |-- costs.py                     #   Cost functions, multiplicative getMirrorCost
|   |-- noise.py                     #   Thermo-optic (JIT), Brownian (Hong), substrate
|   |-- optimizer.py                 #   run_optimization() entry point
|   |-- plot.py                      #   Layers, spectral, noise, corner, starfish
|   |-- mc.py                        #   emcee MCMC sensitivity analysis
|   `-- report.py                    #   Sphinx RST + HTML report generation
|
|-- projects/
|   |-- aLIGO/                       #   SiO2/TiTa2O5 at 1064 nm, 295 K
|   |   |-- materials.yml            #     Material properties (references library)
|   |   |-- ETM_params.yml           #     Cost weights, optimizer settings
|   |   |-- ITM_params.yml
|   |   |-- Data/                    #     HDF5 output (gitignored)
|   |   `-- Figures/                 #     Generated plots (gitignored)
|   |-- Voyager_aSiSiN/             #   aSi/SiN at 2050 nm, 123 K
|   |   |-- materials.yml
|   |   |-- ETM_params.yml
|   |   |-- ITM_params.yml
|   |   |-- Data/
|   |   `-- Figures/
|   `-- Voyager_Ta2O5/              #   Ta2O5/SiO2 at 2050 nm, 123 K
|       |-- materials.yml
|       |-- ETM_params.yml
|       |-- ITM_params.yml
|       |-- Data/
|       `-- Figures/
|
|-- generic/                         # Legacy shared libraries (being replaced by OptimalBragg)
|   |-- coatingUtils.py              #   Transfer matrix, E-field, Sellmeier, YAML import
|   |-- optimUtils.py                #   Cost functions, getMirrorCost evaluator
|   `-- thermoopticUtils.py          #   Numba JIT thermo-optic noise
|
|-- Arms/                            # Legacy project dir (migrating to projects/aLIGO/)
|-- SiN_aSi/                         # Legacy project dir (migrating to projects/Voyager_aSiSiN/)
|-- Ta2O5_Voyager/                   # Legacy project dir (migrating to projects/Voyager_Ta2O5/)
|
|-- tests/                           # Pytest test suite
|   |-- test_materials.py            #   Material class, qw_stack, materials library
|   |-- test_io.py                   #   HDF5/YAML I/O, load_materials_yaml
|   |-- conftest.py                  #   Shared fixtures
|   |-- test_coatingUtils.py         #   Transfer matrix physics tests
|   |-- test_optimUtils.py           #   Cost function tests
|   `-- test_integration.py          #   Full optimization integration test
|
|-- benchmarks/                      # Performance benchmarks
|-- docs/                            # Sphinx documentation + run reports
|   |-- codemap.dot                  #   Graphviz source for architecture diagram
|   |-- codemap.svg                  #   Rendered architecture diagram
|   `-- runs/                        #   Auto-generated run reports
|
|-- environment.yml                  # Conda environment specification
|-- pyproject.toml                   # Package config + pytest settings
|-- .gitlab-ci.yml                   # GitLab CI pipeline
|-- README.md                        # Project overview + install instructions
|-- CODEMAP.md                       # This file
`-- CLAUDE.md                        # AI assistant instructions
```

**Legacy directories** (`Arms/`, `SiN_aSi/`, `Ta2O5_Voyager/`, `generic/`) are being migrated to
`OptimalBragg/` + `projects/`. They remain functional during the transition.

---

## 2. Optimization Pipeline

### Full data flow for a single optimization run

```
                           CONFIGURATION
                    +-----------------------+
                    | ETM_params.yml        |     aSiSiN.yaml (or Ta2O5_ETM.yml)
                    |   costs:              |     Material properties, substrate,
                    |     Trans1064:         |     laser wavelength, mirror geometry
                    |       target / weight |         |
                    |     Brownian: ...     |         |
                    |   misc:              |         |
                    |     Npairs, tol, ...  |         |
                    +-----------+-----------+         |
                                |                     |
                                v                     v
                    +------------------------------- ---------+
                    |              mkETM.py                    |
                    |                                          |
                    |  1. importParams('ETM_params.yml')       |
                    |  2. ifo = gwinc.Struct.from_file(...)    |
                    |  3. gam = brownianProxy(ifo)             |
                    |  4. bounds = [(0.05, 0.48)] * (2N+1)     |
                    |                                          |
                    |  5. scipy.optimize.differential_evolution|
                    |     +----------------------------------+ |
                    |     | func = getMirrorCost             | |
                    |     | strategy = 'best1bin'            | |
                    |     | mutation = (0.05, 1.5)           | |
                    |     | popsize = Nparticles             | |
                    |     | workers = -1 (all cores)         | |
                    |     | maxiter = 2000                   | |
                    |     | polish = True (L-BFGS-B)         | |
                    |     +----------------------------------+ |
                    |                   |                      |
                    |                   | calls per candidate  |
                    |                   v                      |
                    |     +----------------------------------+ |
                    |     |       getMirrorCost(L, ...)       | |
                    |     |  For each cost with weight > 0:  | |
                    |     |    Trans1064  -> transmissionCost  | |
                    |     |    Trans532  -> transmissionCost  | |
                    |     |    TransOPLEV-> transmissionCost  | |
                    |     |    Brownian  -> brownianCost      | |
                    |     |    Thermooptic->thermoopticCost   | |
                    |     |    Lsens     -> sensitivityCost   | |
                    |     |    Esurf     -> surfEfieldCost    | |
                    |     |    Lstdev    -> stdevLCost        | |
                    |     |    Absorption-> (not implemented) | |
                    |     |                                   | |
                    |     |  scalar = sum(weight_i * cost_i)  | |
                    |     +----------------------------------+ |
                    |                                          |
                    |  6. Save to HDF5                         |
                    +------------------+-----------------------+
                                       |
                                       v
                      Data/ETM/ETM_Layers_YYMMDD_HHMMSS.hdf5
                      +------------------------------------+
                      | trajectory     (convergence curve) |
                      | vec_evol       (parameter history) |
                      | diffevo_output/                    |
                      |   L            (optical thickness) |
                      |   n            (refractive indices)|
                      |   scalarCost                       |
                      |   T1064, RPSL, T532, RAUX, ...     |
                      |   vectorCost/{Trans1064, ...}       |
                      +------------------------------------+
                                       |
                    +------------------+------------------+
                    |                                     |
                    v                                     v
          +------------------+                +-----------------------+
          | plot_ETM.py      |                | doMC.py               |
          | Design dashboard:|                | Monte Carlo analysis: |
          |  - Layer profile |                |  20 walkers, 4D,     |
          |  - Spectral R(l) |                |  0.5% Gaussian       |
          |  - Cost summary  |                |  perturbations       |
          +------------------+                +-----------+-----------+
                                                          |
                                                          v
                                              MC_output.hdf5
                                                          |
                                                          v
                                              +-----------------------+
                                              | cornerPlt.py          |
                                              | Corner plots:         |
                                              |  T_PSL vs T_AUX vs   |
                                              |  surface field, etc.  |
                                              +-----------------------+
```

### Key internal call: `getMirrorCost`

```
getMirrorCost(L, costs, ifo, gam, verbose, misc)
    |
    |-- Extend L with Ncopies and Nfixed bilayers
    |-- Build refractive index array n = [1, nL, nH, nL, nH, ..., nSub]
    |
    |-- For each cost term with weight > 0:
    |
    |   Trans1064 -----> transmissionCost(target, n, L, lamb=1.0)
    |                       \-> multidiel1(n, L, lamb, theta, pol)
    |
    |   Trans532 -----> transmissionCost(target, n, L, lamb=1550/2050)
    |                       \-> multidiel1(...)
    |
    |   TransOPLEV ---> transmissionCost(target, n, L, lamb=0.297)
    |                       \-> multidiel1(...)
    |
    |   Brownian -----> brownianCost(target, L, gam)
    |                       (proxy: zLow + gam * zHigh)
    |
    |   Thermooptic --> thermoopticCost(target, fTO, L, ifo)
    |                       \-> gwinc.noise.coatingthermal.coating_thermooptic(...)
    |
    |   Lsens --------> sensitivityCost(target, n, L) at PSL and AUX
    |                       \-> transmissionCost(target, n, 1.01*L)
    |                           (quadrature sum of both wavelengths)
    |
    |   Esurf --------> surfEfieldCost(target, n, L)
    |                       \-> multidiel1(...)
    |                       (formula: 50 * arcsinh(|1+r|^2))
    |
    |   Lstdev -------> stdevLCost(target, L)
    |                       (mean/std ratio penalty)
    |
    |-- scalar_cost = sum( weight_i * cost_i )
    |-- return scalar_cost  (or also verbose dict)
```

---

## 3. Core Library Function Map

### `coatingUtils.py` (generic/)

| Function | Signature | Purpose | Called by |
|----------|-----------|---------|----------|
| `multidiel1` | `(n, L, lamb, theta=0, pol='te')` | Transfer matrix reflectivity of dielectric stack | `transmissionCost`, `sensitivityCost`, `surfEfieldCost`, `specREFL`, `doMC.py` |
| `op2phys` | `(L, n)` | Convert optical thickness to physical thickness | `doMC.py`, plotting scripts |
| `lnprob` | `(x, mu, icov)` | Log-probability for Gaussian (MC sampling) | `doMC.py` (emcee sampler) |
| `surfaceField` | `(gamm, Ei=27.46)` | Surface E-field from amplitude reflectivity | `doMC.py` |
| `specREFL` | `(layers, dispFileName, lambda_0, lam, aoi, pol)` | Spectral reflectivity with dispersion | Plotting scripts |
| `fieldDepth` | `(L, n, lam=1064e-9, theta=0, pol='s', nPts=30)` | E-field squared vs depth (Arnon & Baumeister 1980) | `plotlayers.py`, `calcAbsorption` (indirectly) |
| `importParams` | `(paramFile)` | Load YAML config file | `mkETM.py`, `mkITM.py`, `getMirrorCost` (generic only) |
| `calcAbsorption` | `(Esq, L, nPts, alphaOdd, alphaEven)` | Integrated absorption from E-field profile [ppm] | Analysis scripts |
| `sellmeier` | `(B, C, lam=1064e-9)` | Refractive index from Sellmeier coefficients | Dispersion analysis scripts |

### `optimUtils.py` (generic/)

| Function | Signature | Purpose | Called by |
|----------|-----------|---------|----------|
| `transmissionCost` | `(target, n, L, lamb=1, theta=0, pol='te')` | Transmission cost: `\|target - T\|^2 / target^2` | `getMirrorCost` |
| `sensitivityCost` | `(target, n, L, lamb=1, theta=0, pol='te')` | Sensitivity to 1% thickness perturbation | `getMirrorCost` |
| `surfEfieldCost` | `(target, n, L, lamb=1, theta=0, pol='te')` | Surface E-field cost: `50 * arcsinh(\|1+r\|^2)` | `getMirrorCost` |
| `stdevLCost` | `(target, L)` | Thickness uniformity penalty (mean/std ratio) | `getMirrorCost` |
| `brownianProxy` | `(ifo)` | Pre-compute Brownian noise prefactor gamma | `mkETM.py`, `mkITM.py` |
| `brownianCost` | `(target, L, gam)` | Brownian noise proxy: `target * (zLow + gam * zHigh)` | `getMirrorCost` |
| `thermoopticCost` | `(target, fTarget, L, ifo)` | Thermo-optic noise via pygwinc | `getMirrorCost` |
| `getMirrorCost` | `(L, costs, ifo, gam, verbose=False, misc={})` | Master cost evaluator; weighted sum of all sub-costs | `devo()` in `mkETM.py`/`mkITM.py` |

**Note:** All projects import directly from `generic/` using the dict-style `getMirrorCost` signature. The `multidiel1` inner loop is Numba JIT-compiled for performance.

---

## 4. Cost Function Breakdown

All costs are computed inside `getMirrorCost`. Each cost term is independently weighted; set `weight: 0` to disable.

| Cost key | Physics | Formula / Method | Computed by |
|----------|---------|-----------------|-------------|
| **Trans1064** | Transmission at primary laser wavelength | `\|target - T\|^2 / target^2` where `T = 1 - \|r\|^2` at `lamb=1.0` | `transmissionCost` -> `multidiel1` |
| **Trans532** | Transmission at auxiliary wavelength | Same formula at `lamb = 1550/2050` (SiN_aSi) or `lamb = 1418.8/2128.2` (Ta2O5) | `transmissionCost` -> `multidiel1` |
| **TransOPLEV** | Transmission at optical lever wavelength | Same formula at `lamb = 0.297` (~608 nm) | `transmissionCost` -> `multidiel1` |
| **Brownian** | Coating Brownian thermal noise | Proxy: `target * (sum_low + gam * sum_high)` per E0900068 p.4 | `brownianCost` |
| **Thermooptic** | Thermo-optic noise (thermoelastic + thermorefractive) | `target * S_TO(f)` via pygwinc `coating_thermooptic` | `thermoopticCost` -> `gwinc` |
| **Lsens** | Sensitivity of T to 1% layer thickness error | Quadrature sum of PSL and AUX sensitivities: `sqrt(s_PSL^2 + s_AUX^2)` | `sensitivityCost` -> `transmissionCost` |
| **Esurf** | Surface electric field at HR face | `50 * arcsinh(\|1 + r\|^2)` | `surfEfieldCost` -> `multidiel1` |
| **Lstdev** | Layer thickness uniformity | `\|target - mean(L)/std(L)\|^2 / target^2` | `stdevLCost` |
| **Absorption** | Integrated coating absorption | **Not implemented** in optimizer (placeholder; `calcAbsorption` exists in `coatingUtils` for post-analysis) | -- |

### Scalar cost aggregation

```
scalar_cost = prod_i ( 1 + weight_i * cost_i )
```

Multiplicative cost function: each factor >= 1, so all objectives must be satisfied simultaneously.
The optimizer minimizes this scalar. After convergence, `polish=True` refines with L-BFGS-B.

---

## 5. Project Comparison

| | **SiN_aSi** | **Ta2O5_Voyager** | **aSi_Voyager** |
|---|---|---|---|
| **Status** | Active | Active | Dormant |
| **High-n material** | a-Si (n=3.65) | Ta2O5 (n=2.083) | a-Si |
| **Low-n material** | SiN (n=2.17) | SiO2 (n=1.435) | -- |
| **Primary wavelength** | 2050.15 nm | 2128.2 nm | -- |
| **Auxiliary wavelength** | 1550 nm | 1418.8 nm | -- |
| **Temperature** | 123 K | 123 K | -- |
| **ETM layer pairs** | 14 | 22 | 9 |
| **ITM layer pairs** | 19 | 19 | -- |
| **ETM popsize** | 200 | 10 | 20 |
| **ITM popsize** | 12 | 12 | -- |
| **ETM init method** | halton | latinhypercube | -- |
| **Config format** | New (dict-style) | New (dict-style) | Old (list-style) |
| **GWINC struct file** | aSiSiN.yaml | Ta2O5_ETM.yml / Ta2O5_ITM.yml | aSiModel.yaml |
| **Entry points** | mkETM.py, mkITM.py | mkETM.py, mkITM.py | mkMirror.py |
| **Output format** | HDF5 | HDF5 | .mat (legacy) |
| **Has doMC.py** | Yes | Yes | No |

---

## 6. Config Format Reference

### New dict-style format (`ETM_params.yml`, `ITM_params.yml`)

Used by all active projects. Loaded by `importParams()`, passed directly to `getMirrorCost`.

```yaml
# Cost function terms: each key maps to {target, weight}
# Set weight: 0 to disable a term
costs:
    Trans1064:                # Transmission at primary laser wavelength
        target: 5e-6         # Target value (5 ppm)
        weight: 5            # Relative weight in scalar cost
    Brownian:                # Brownian thermal noise proxy
        target: 0.1
        weight: 2
    Thermooptic:             # Thermo-optic noise
        target: 1e43         # Large target = small cost contribution
        weight: 2
    Lsens:                   # Layer thickness sensitivity
        target: 1e-7
        weight: 0            # Disabled
    Esurf:                   # Surface E-field
        target: 1e-9
        weight: 0
    Absorption:              # Integrated absorption (NOT IMPLEMENTED)
        target: 1e-4
        weight: 0
    Trans532:                # Transmission at auxiliary wavelength
        target: 1000e-6
        weight: 0
    TransOPLEV:              # Transmission at optical lever wavelength
        target: 0.05
        weight: 0
    Lstdev:                  # Layer thickness uniformity
        target: 0.5
        weight: 0

# Miscellaneous optimizer and physics parameters
misc:
    fTO: 100                 # Frequency for TO noise evaluation [Hz]
    pol: 'te'                # Polarization ('te' = s-pol, 'tm' = p-pol)
    aoi: 0                   # Angle of incidence [degrees]
    Npairs: 14               # Number of bilayer pairs to optimize
    Nfixed: 0                # Fixed bilayers appended after variable stack
    Ncopies: 0               # Copies of variable stack to append
    Nparticles: 200          # Population size for differential evolution
    atol: 3e-9               # Absolute tolerance for convergence
    tol: 1e-2                # Relative tolerance for convergence
    init_method: 'halton'    # Population initialization ('halton' or 'latinhypercube')
    gwincStructFile: 'aSiSiN.yaml'  # Path to GWINC material properties
```

### Old list-style format (`params.yml` in aSi_Voyager)

Legacy format. **Not used by active `mkETM.py`/`mkITM.py` scripts.**

```yaml
costs:    [Trans, coatBr, coatTO, sensL, surfE]
targets:  [5e-6, 1, 1e+42, 1e-2, 0]
weights:  [0.3, 1, 1e-9, 0.03, 1]
Npairs:   9
Nparticles: 20
fTO:      100
aoi:      0
pol:      'te'
gwincStructFile: aSiModel.m
```

---

## 7. Known Architecture Issues

1. **Absorption cost not implemented** -- The `Absorption` key is accepted in the config but does nothing in the optimizer. The underlying `calcAbsorption` function exists in `coatingUtils` for post-analysis only.

2. **`specREFL` hardcodes Ta2O5/SiO2 dispersion** -- The function loads dispersion data from a `.mat` file keyed to `'SiO2'` and `'Ta2O5'`, making it unusable for SiN/a-Si materials without modification.

### Resolved Issues (refactor/speed-and-robustness branch)

- ~~`generic_local/` duplication~~ -- Eliminated. All projects import from canonical `generic/`.
- ~~Old vs new `getMirrorCost` signature~~ -- Unified to dict-style API.
- ~~Hardcoded AUX wavelength ratios~~ -- Now configurable via `lambdaAUX` in params YAML.
- ~~`doMC.py` reads `.mat` files~~ -- Updated to read HDF5 output.
- ~~No `__init__.py`~~ -- `generic/__init__.py` added; proper package structure with `pyproject.toml`.
- ~~README describes only MATLAB workflows~~ -- Updated with Python pipeline docs.
