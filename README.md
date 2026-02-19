# Design of Mirror Coatings by global optimization of a cost function

Code to optimize HR mirror dielectric stack designs for gravitational wave detector test masses.
Designs are optimized via `scipy.optimize.differential_evolution`, balancing transmissivity,
thermal noise, manufacturing tolerance, surface E-field, and absorption constraints.

For a detailed architecture reference, see **[CODEMAP.md](CODEMAP.md)**.

## How to Install

1. Clone the repository:
   ```bash
   git clone <repo-url>
   cd Coatings
   ```
2. Create the conda environment and install the package:
   ```bash
   conda env create -f environment.yml
   conda activate coatingDev
   pip install -e .
   ```
   Key dependencies: numpy, scipy, matplotlib, numba (>=0.56), emcee, corner, gwinc (>=0.6), lmfit, h5py, pytest.

3. (Optional) For legacy MATLAB workflows, install the MATLAB engine for Python:
   ```bash
   cd $MATLABROOT/extern/engines/python
   python setup.py install
   ```
   Set `GWINCPATH` in your shell (e.g. `conda env config vars set GWINCPATH=<path_to_matgwinc>`).

## Quick Start

### Run a coating optimization

From a project directory (e.g. `SiN_aSi/` or `Ta2O5_Voyager/`):

```bash
python mkETM.py          # Optimize End Test Mass
python mkITM.py          # Optimize Input Test Mass
```

Output: `Data/ETM/ETM_Layers_YYMMDD_HHMMSS.hdf5`

### Monte Carlo sensitivity analysis

```bash
python doMC.py Data/ETM/ETM_Layers_YYMMDD_HHMMSS.hdf5 output_MC.hdf5 5000
```

Uses `emcee` (20 walkers, 4D parameter space) to perturb angle of incidence, refractive indices, and layer thicknesses by 0.5% Gaussian.

### Visualization

```bash
python plot_ETM.py       # Design analysis dashboard
python cornerPlt.py      # MC corner plots from .hdf5
python plotlayers.py     # Layer structure + E-field
```

## Optimal Objectives

For the mirror coatings, we have many constraints to satisfy:

1. Transmissivity at the primary laser wavelength (e.g. 5 ppm at 2050 nm).
2. Transmissivity at auxiliary wavelengths (1550 nm, optical lever).
3. Minimize Brownian thermal noise.
4. Minimize Thermo-Optic noise.
5. Minimize sensitivity of transmissivity to coating deposition errors.
6. Minimize E-field at HR surface.
7. Layer thickness uniformity.

Each objective is a weighted term in a scalar cost function. See [CODEMAP.md](CODEMAP.md) for the full cost function breakdown and architecture details.

## Active Projects

| Project | Materials | Primary wavelength | Temperature |
|---------|-----------|-------------------|-------------|
| `SiN_aSi/` | a-Si / SiN | 2050 nm | 123 K |
| `Ta2O5_Voyager/` | Ta2O5 / SiO2 | 2128 nm | 123 K |

## Running Tests

```bash
pytest tests/ -v -k "not slow"           # Unit tests (~1 sec)
pytest tests/test_integration.py -m slow  # Integration test (~30 sec)
python benchmarks/bench_multidiel1.py     # Performance benchmark
```

## Paper draft

A paper draft of this work lives at [this git repo](https://github.com/CaltechExperimentalGravity/OptimalCoatingDesign).
