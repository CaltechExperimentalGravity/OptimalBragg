Getting Started
===============

Installation
------------

1. Clone the repository::

    git clone git@github.com:CaltechExperimentalGravity/OptimalBragg.git
    cd OptimalBragg

2. Create the conda environment and install::

    conda env create -f environment.yml
    conda activate coatingDev
    pip install -e .

Running an Optimization
-----------------------

Using the CLI (from the repo root)::

    optimalbragg optimize projects/aLIGO/ETM_params.yml

Or from a project directory::

    cd projects/aLIGO
    python mkETM.py

Output is saved to ``Data/ETM/ETM_Layers_YYMMDD_HHMMSS.hdf5``.

The optimizer uses ``scipy.optimize.differential_evolution`` with
``workers=1`` (IPC overhead exceeds per-eval cost). A typical run with
18 bilayer pairs and population size 500 takes ~3 minutes.

Configuration
^^^^^^^^^^^^^

Each project directory has two YAML config files:

- **``ETM_params.yml``** (or ``ITM_params.yml``) configures cost function
  weights/targets, optimizer settings, and physics parameters.

- **``materials.yml``** defines material properties, referencing the central
  ``OptimalBragg.materials`` library with optional per-project overrides.

See ``projects/aLIGO/ETM_params.yml`` and ``projects/aLIGO/materials.yml``
for fully commented examples.

Reading Results
---------------

The HDF5 output contains:

- ``diffevo_output/L`` -- optimized optical thicknesses
- ``diffevo_output/n`` -- refractive index array (air + layers + substrate)
- ``diffevo_output/scalarCost`` -- final scalar cost
- ``diffevo_output/T1064``, ``RPSL`` -- transmission and reflectivity at PSL
- ``diffevo_output/vectorCost/*`` -- individual cost term values
- ``trajectory`` -- convergence curve

Visualization::

    optimalbragg plot Data/ETM/ETM_Layers_YYMMDD_HHMMSS.hdf5

Monte Carlo Analysis
--------------------

Perturb refractive indices and layer thicknesses to assess manufacturing
sensitivity::

    optimalbragg mc Data/ETM/ETM_Layers_YYMMDD_HHMMSS.hdf5 5000

Uses ``emcee`` with 20 walkers in a 3-dimensional parameter space
(high-n index, low-n index, thickness), each perturbed as 0.5% Gaussian.

Generate corner plots::

    optimalbragg corner Data/ETM/ETM_MC.hdf5

Running Tests
-------------

::

    pytest tests/ -v -k "not slow"           # Unit tests (~2 sec)
    pytest tests/test_integration.py -m slow  # Integration test
    python benchmarks/bench_multidiel1.py     # Transfer matrix benchmark
    python benchmarks/bench_thermooptic.py    # Thermooptic JIT vs numpy

Active Projects
---------------

==========================  ================  ==================  ===========
Project                     Materials         Primary wavelength  Temperature
==========================  ================  ==================  ===========
``projects/aLIGO/``         SiO2 / TiTa2O5   1064 nm             295 K
``projects/SFG/``           SiO2 / Ta2O5     1064 nm             295 K
``projects/Voyager_aSiSiN/``  a-Si / SiN     2050 nm             123 K
``projects/Voyager_Ta2O5/``   Ta2O5 / SiO2   2050 nm             123 K
==========================  ================  ==================  ===========
