# Code Map: OptimalBragg — Coating Optimization for Gravitational Wave Detectors

Comprehensive architecture reference for the optical coating optimization codebase.
For quick-start instructions, see [README.md](README.md).

---

## Architecture Overview

![OptimalBragg codemap](docs/codemap.svg)

---

## 1. Repository Layout

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
|   |-- report.py                    #   Sphinx RST + HTML report generation
|   |-- cli.py                       #   CLI entry point (optimalbragg command)
|   `-- __main__.py                  #   python -m OptimalBragg support
|
|-- projects/
|   |-- aLIGO/                       #   SiO2/TiTa2O5 at 1064 nm, 295 K
|   |   |-- materials.yml            #     Material properties (references library)
|   |   |-- ETM_params.yml           #     Cost weights, optimizer settings
|   |   |-- ITM_params.yml
|   |   |-- mkETM.py / mkITM.py      #     Thin wrappers calling run_optimization()
|   |   |-- Data/                    #     HDF5 output (gitignored)
|   |   `-- Figures/                 #     Generated plots (gitignored)
|   |-- Voyager_aSiSiN/             #   aSi/SiN at 2050 nm, 123 K
|   `-- Voyager_Ta2O5/              #   Ta2O5/SiO2 at 2050 nm, 123 K
|
|-- tests/                           # Pytest test suite
|   |-- test_materials.py            #   Material class, qw_stack, materials library
|   |-- test_layers.py               #   Transfer matrix, E-field, spectral
|   |-- test_costs.py                #   Cost functions, getMirrorCost
|   |-- test_noise.py                #   All thermal noise models
|   |-- test_thermooptic.py          #   JIT vs numpy thermo-optic validation
|   |-- test_io.py                   #   HDF5/YAML I/O, load_materials_yaml
|   |-- test_optimizer.py            #   Stack building, optimization pipeline
|   |-- test_plot.py                 #   Plotting functions
|   |-- test_mc.py                   #   Monte Carlo perturbation generation
|   |-- test_report.py              #   Report quality assessment
|   |-- test_cli.py                  #   CLI help/argument parsing
|   |-- test_integration.py          #   Full optimization integration test
|   `-- conftest.py                  #   Shared fixtures
|
|-- benchmarks/                      # Performance benchmarks
|   |-- bench_multidiel1.py          #   Transfer matrix micro-benchmark
|   `-- bench_thermooptic.py         #   JIT vs numpy thermooptic comparison
|
|-- docs/                            # Sphinx documentation + run reports
|   |-- conf.py                      #   Sphinx config (pydata theme)
|   |-- index.rst                    #   Documentation root
|   |-- api/                         #   Autodoc API reference stubs
|   |-- runs/                        #   Auto-generated run reports
|   |-- codemap.dot                  #   Graphviz source for architecture diagram
|   `-- codemap.svg                  #   Rendered architecture diagram
|
|-- environment.yml                  # Conda environment specification
|-- pyproject.toml                   # Package config + pytest settings
|-- .github/workflows/ci.yml         # GitHub Actions CI + Pages
|-- README.md                        # Project overview + install instructions
|-- CODEMAP.md                       # This file
`-- CLAUDE.md                        # AI assistant instructions
```

**Legacy directories** (`Arms/`, `SiN_aSi/`, `Ta2O5_Voyager/`, `generic/`) contain the original code
before the OptimalBragg refactor. They are no longer used by the active pipeline.

---

## 2. Optimization Pipeline

### Full data flow

```
                         CONFIGURATION
              +------------------------+     +--------------------+
              | ETM_params.yml         |     | materials.yml      |
              |   costs:               |     |   substrate: SiO2  |
              |     Trans1:            |     |   thin_films:      |
              |       target / weight  |     |     L: SiO2        |
              |     Brownian: ...      |     |     H: TiTa2O5     |
              |   misc:               |     |   laser:           |
              |     Npairs, tol, ...   |     |     wavelength     |
              +----------+-------------+     +--------+-----------+
                         |                            |
                         v                            v
              +----------------------------------------------+
              |         run_optimization(params_yaml)         |
              |                                                |
              |  1. yamlread(params_yaml)                      |
              |  2. load_materials_yaml(materials_file)         |
              |  3. qw_stack(lam_ref, substrate, thin_films)    |
              |  4. gam = brownian_proxy(stack)                 |
              |  5. precompute_misc(costs, stack, misc)         |
              |                                                |
              |  6. scipy.optimize.differential_evolution       |
              |     +----------------------------------------+ |
              |     | func = getMirrorCost                   | |
              |     | strategy = 'best1bin'                  | |
              |     | workers = 1                            | |
              |     | maxiter = 2000                         | |
              |     | polish = True (L-BFGS-B)               | |
              |     +----------------------------------------+ |
              |                   |                            |
              |                   | calls per candidate        |
              |                   v                            |
              |     +----------------------------------------+ |
              |     |       getMirrorCost(L, ...)             | |
              |     |  For each cost with weight > 0:        | |
              |     |    Trans1     -> transmissionCost        | |
              |     |    Trans2     -> transmissionCost        | |
              |     |    Brownian   -> brownianCost            | |
              |     |    Thermooptic-> thermoopticCost         | |
              |     |    Absorption -> absorptionCost          | |
              |     |    Lsens      -> sensitivityCost         | |
              |     |    Esurf      -> surfEfieldCost          | |
              |     |    Lstdev     -> stdevLCost              | |
              |     |                                         | |
              |     |  scalar = prod(1 + weight_i * cost_i)   | |
              |     +----------------------------------------+ |
              |                                                |
              |  7. Save to HDF5                               |
              +---------------------+--------------------------+
                                    |
                                    v
                   Data/ETM/ETM_Layers_YYMMDD_HHMMSS.hdf5
                                    |
                 +------------------+------------------+
                 |                                     |
                 v                                     v
       +------------------+               +-----------------------+
       | optimalbragg plot |               | optimalbragg mc       |
       | Design dashboard: |               | Monte Carlo analysis: |
       |  - Layer profile  |               |  20 walkers, 3D,     |
       |  - Spectral R(l)  |               |  0.5% Gaussian       |
       |  - Noise budget   |               |  perturbations       |
       |  - Starfish       |               |  (n_H, n_L, dL)      |
       +------------------+               +-----------+-----------+
                                                       |
                                                       v
                                            MC_output.hdf5
                                                       |
                                                       v
                                           +-----------------------+
                                           | optimalbragg corner   |
                                           | ArviZ corner plots:   |
                                           |  T_PSL vs T_AUX vs   |
                                           |  surface field, etc.  |
                                           +-----------------------+
```

### Key internal call: `getMirrorCost`

```
getMirrorCost(L, costs, stack, gam, verbose, misc)
    |
    |-- Extend L with Ncopies and Nfixed bilayers
    |-- Build refractive index array n from stack dict
    |
    |-- Consolidated multidiel1 call (all wavelengths at once)
    |
    |-- For each cost term with weight > 0:
    |
    |   Trans1 -------> transmissionCost(target, n, L, lamb=1.0)
    |   Trans2 -------> transmissionCost(target, n, L, lamb=lambda2)
    |   Trans3 -------> transmissionCost(target, n, L, lamb=lambda3)
    |   Brownian ------> brownianCost(target, L, gam)
    |   Thermooptic ---> thermoopticCost(target, fTO, L, stack, ...)
    |   Absorption ----> absorptionCost(target, n, L, stack, ...)
    |   Lsens ---------> sensitivityCost(target, n, L) at PSL and AUX
    |   Esurf ---------> surfEfieldCost(target, n, L)
    |   Lstdev --------> stdevLCost(target, L)
    |
    |-- scalar_cost = prod( 1 + weight_i * cost_i )
    |-- return scalar_cost  (or also verbose dict)
```

---

## 3. Core Library Function Map

### `OptimalBragg/layers.py`

| Function | Purpose | Performance |
|----------|---------|-------------|
| `multidiel1(n, L, lamb, theta, pol)` | Transfer matrix reflectivity of dielectric stack | Numba JIT, ~7 us/call |
| `multilayer_diel(ns, Ls, lamb, aoi, pol)` | Wrapper converting physical to optical thicknesses | |
| `field_zmag(ns, Ls, lam, ...)` | E-field depth profile (Arnon & Baumeister) | |
| `surfield(rr, Ei, normalized)` | Surface E-field from reflectivity | |
| `calc_abs(Esq, Ls, alphas)` | Integrated absorption from E-field profile | |
| `op2phys(L, n)` | Optical to physical thickness conversion | |
| `field_zmag(L, n, ...)` | E-field squared vs depth | |
| `sellmeier(B, C, lam)` | Refractive index from Sellmeier coefficients | |
| `amp_refl(wavelengths, stack)` | Spectral amplitude reflectivity | |
| `refl(wavelengths, stack)` | Spectral power reflectivity | |
| `trans(wavelengths, stack)` | Spectral transmission | |

### `OptimalBragg/costs.py`

| Function | Purpose |
|----------|---------|
| `transmissionCost(target, n, L, lamb, theta, pol)` | `\|target - T\|^2 / target^2` |
| `brownianCost(target, L, gam)` | Proxy: `target * (zLow + gam * zHigh)` |
| `thermoopticCost(target, fTarget, L, stack, ...)` | Thermo-optic via JIT `coating_thermooptic_fast` |
| `absorptionCost(target, n, L, stack, ...)` | Integrated coating absorption |
| `sensitivityCost(target, n, L, ...)` | 1% thickness perturbation sensitivity |
| `surfEfieldCost(target, n, L, ...)` | `50 * arcsinh(\|1+r\|^2)` |
| `stdevLCost(target, L)` | Thickness uniformity (mean/std ratio) |
| `precompute_misc(costs, stack, misc)` | Cache invariant computations |
| `getMirrorCost(L, costs, stack, gam, verbose, misc)` | Master cost evaluator |

### `OptimalBragg/noise.py`

| Function | Physics | Performance |
|----------|---------|-------------|
| `coating_thermooptic_fast(f, L, lam, w, params)` | Thermo-optic PSD (single freq) | Numba JIT |
| `coating_thermooptic(ff, stack, w, r, d)` | Thermo-optic PSD (vectorized) | Pure numpy |
| `coating_brownian(f, stack, w)` | Hong et al. PRD 87, 082001 (2013) | |
| `substrate_brownian(f, stack, w)` | Substrate Brownian noise | |
| `substrate_thermoelastic(f, stack, w)` | Substrate thermoelastic noise | |
| `substrate_thermorefractive(f, stack, w)` | Substrate thermorefractive noise | |
| `brownian_proxy(stack)` | Pre-compute Brownian proxy gamma | |

### `OptimalBragg/materials.py`

Materials with full reference citations:

| Material | n | Use case |
|----------|---|----------|
| `SiO2` | 1.45 | Low-n coating, room temp |
| `TiTa2O5` | 2.06 | High-n coating (aLIGO), room temp |
| `Ta2O5` | 2.0 | High-n coating, room temp |
| `FusedSilica` | 1.45 | Substrate, room temp |
| `cSi_123` | 3.39 | Silicon substrate, 123 K |
| `aSi_123` | 3.65 | High-n coating (Voyager), 123 K |
| `SiN_123` | 2.17 | Low-n coating (Voyager), 123 K |
| `air` | 1.0 | Superstrate |

---

## 4. Cost Function Breakdown

All costs are computed inside `getMirrorCost`. Each cost term is independently weighted; set `weight: 0` to disable.

| Cost key | Physics | Formula / Method |
|----------|---------|-----------------|
| **Trans1** | Transmission at design wavelength | `\|target - T\|^2 / target^2` |
| **Trans2** | Transmission at second wavelength | Same formula at `lamb = lambda2` |
| **Trans3** | Transmission at third wavelength | Same formula at `lamb = lambda3` |
| **Brownian** | Coating Brownian thermal noise | Proxy: `target * (sum_low + gam * sum_high)` per E0900068 |
| **Thermooptic** | Thermo-optic noise | `target * S_TO(f)` via Numba JIT |
| **Absorption** | Integrated coating absorption | E-field profile integration |
| **Lsens** | Sensitivity of T to layer thickness error | Quadrature sum of PSL and AUX |
| **Esurf** | Surface electric field at HR face | `50 * arcsinh(\|1+r\|^2)` |
| **Lstdev** | Layer thickness uniformity | `\|target - mean(L)/std(L)\|^2 / target^2` |

### Scalar cost aggregation

```
scalar_cost = prod_i ( 1 + weight_i * cost_i )
```

Multiplicative cost function: each factor >= 1, so all objectives must be satisfied simultaneously.
The optimizer minimizes this scalar. After convergence, `polish=True` refines with L-BFGS-B.

---

## 5. Active Projects

| | **aLIGO** | **Voyager aSiSiN** | **Voyager Ta2O5** |
|---|---|---|---|
| **Directory** | `projects/aLIGO/` | `projects/Voyager_aSiSiN/` | `projects/Voyager_Ta2O5/` |
| **High-n material** | TiTa2O5 (n=2.06) | a-Si (n=3.65) | Ta2O5 (n=2.0) |
| **Low-n material** | SiO2 (n=1.45) | SiN (n=2.17) | SiO2 (n=1.435) |
| **Primary wavelength** | 1064 nm | 2050 nm | 2050 nm |
| **Temperature** | 295 K | 123 K | 123 K |
| **Substrate** | Fused silica | c-Si | c-Si |

---

## 6. Config Format Reference

### `ETM_params.yml` (cost function and optimizer settings)

```yaml
costs:
    Trans1:
        target: 5e-6
        weight: 15
    Brownian:
        target: 20.0
        weight: 2
    Thermooptic:
        target: 1.6e-42
        weight: 2
    # ... other terms with weight: 0 to disable

misc:
    Npairs: 18
    Nparticles: 500
    atol: 1e-10
    tol: 1e-3
    init_method: halton
    lambda2: 0.5
    materials_file: materials.yml
```

### `materials.yml` (material properties)

```yaml
substrate:
  material: FusedSilica        # from OptimalBragg.materials library
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
  ITM:
    beam_radius: 0.055
```

The `load_materials_yaml()` function resolves material names against the library,
applies overrides, and returns Material objects ready for `qw_stack()`.
