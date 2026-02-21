Getting Started
===============

Installation
------------

1. Clone the repository::

    git clone git@git.ligo.org:40m/Coatings.git
    cd Coatings

2. Create the conda environment and install::

    conda env create -f environment.yml
    conda activate coatingDev
    pip install -e .

Running an Optimization
-----------------------

From the ``SiN_aSi/`` project directory::

    python mkETM.py    # Optimize End Test Mass
    python mkITM.py    # Optimize Input Test Mass

Output is saved to ``Data/ETM/ETM_Layers_YYMMDD_HHMMSS.hdf5``.

The optimizer uses ``scipy.optimize.differential_evolution`` with
``workers=1`` (IPC overhead exceeds per-eval cost). A typical run with
17 bilayer pairs and population size 500 takes ~3 minutes.

Configuration
^^^^^^^^^^^^^

Each project directory has a ``ETM_params.yml`` (or ``ITM_params.yml``)
that configures:

- **Cost function weights and targets** -- which objectives are active and
  how they trade off
- **Optimizer settings** -- population size, tolerances, initialization method
- **Physics parameters** -- angle of incidence, polarization, auxiliary
  wavelength ratio

See ``SiN_aSi/ETM_params.yml`` for a fully commented example.

Reading Results
---------------

The HDF5 output contains:

- ``diffevo_output/L`` -- optimized optical thicknesses
- ``diffevo_output/n`` -- refractive index array (air + layers + substrate)
- ``diffevo_output/scalarCost`` -- final scalar cost
- ``diffevo_output/T1064``, ``RPSL`` -- transmission and reflectivity at PSL
- ``diffevo_output/vectorCost/*`` -- individual cost term values
- ``trajectory`` -- convergence curve
- ``gwincStructFile`` -- path to the gwinc material properties file

Visualization::

    python plot_ETM.py      # Design analysis dashboard
    python plotlayers.py    # Layer structure + E-field profile

Monte Carlo Analysis
--------------------

Perturb material properties and thicknesses to assess manufacturing
sensitivity::

    python doMC.py Data/ETM/ETM_Layers_YYMMDD_HHMMSS.hdf5 output_MC.hdf5 5000

Uses ``emcee`` with 20 walkers in a 4-dimensional parameter space
(angle of incidence, high-n index, low-n index, thickness), each
perturbed as 0.5% Gaussian.

Running Tests
-------------

::

    pytest tests/ -v -k "not slow"          # Unit tests (~1 sec)
    pytest tests/test_integration.py -m slow  # Integration test
    python benchmarks/bench_multidiel1.py     # Performance benchmark

Active Projects
---------------

==================  ================  ==================  ===========
Project             Materials         Primary wavelength  Temperature
==================  ================  ==================  ===========
``Arms/``           SiO2 / Ti:Ta2O5   1064 nm             295 K
``SiN_aSi/``        a-Si / SiN        2050 nm             123 K
``Ta2O5_Voyager/``  Ta2O5 / SiO2      2128 nm             123 K
==================  ================  ==================  ===========

The ``Arms/`` project targets aLIGO arm cavity mirrors (ETM and ITM)
with SiO2/Ti:Ta2O5 coatings on fused silica substrates. It includes a
532 nm auxiliary wavelength constraint for green locking, with
``lambdaAUX: 0.5`` in the params YAML.

Running aLIGO optimizations::

    cd Arms
    python mkETM.py    # End Test Mass
    python mkITM.py    # Input Test Mass
    python plot_ETM.py # Generate plots + Sphinx run report
