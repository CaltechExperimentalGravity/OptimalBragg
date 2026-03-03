Cost Function
=============

The optimizer minimizes a scalar cost function that is a multiplicative
product of individual objective terms. Each term is configured in the
project's ``ETM_params.yml`` with a ``target`` and ``weight``. Set
``weight: 0`` to disable a term.

.. math::

   C_\text{total} = \prod_i (1 + w_i \cdot c_i)

The master evaluator is :func:`~OptimalBragg.costs.getMirrorCost`.

Cost Terms
----------

Trans1064 -- Primary Laser Transmission
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Targets a specific power transmission at the primary laser wavelength
(e.g., 5 ppm at 2050 nm for an ETM):

.. math::

   c_\text{Trans1064} = \left|\frac{T_\text{target} - T}{T_\text{target}}\right|^2

where :math:`T = 1 - |r|^2` from :func:`~OptimalBragg.layers.multidiel1`.

Trans532 -- Auxiliary Wavelength Transmission
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Same formula as Trans1064, evaluated at the auxiliary wavelength ratio
``lambdaAUX`` (configured per project in the YAML ``misc`` section).

TransOPLEV -- Optical Lever Transmission
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Same formula, evaluated at ``lamb = 0.297`` (~608 nm for a 2050 nm design).

Brownian -- Coating Brownian Thermal Noise
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A fast proxy for Brownian thermal noise (LIGO-E0900068 p. 4):

.. math::

   c_\text{Brownian} = t \cdot (z_\text{low} + \gamma \cdot z_\text{high})

where :math:`z_\text{low}` and :math:`z_\text{high}` are the total
optical thicknesses of low-index and high-index layers, :math:`\gamma`
encodes material loss angles and Young's moduli, and :math:`t` is the
target scaling factor.

The pre-factor :math:`\gamma` is computed once before optimization by
:func:`~OptimalBragg.noise.brownian_proxy`.

Thermooptic -- Thermo-Optic Noise
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Evaluates the full thermo-optic noise PSD at a single frequency using
pygwinc:

.. math::

   c_\text{TO} = t \cdot S_\text{TO}(f_\text{target})

This is the most expensive cost term per evaluation, as it calls into
the gwinc thermal noise calculation.

Lsens -- Layer Thickness Sensitivity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Measures how sensitive the transmission is to a 1% perturbation in all
layer thicknesses:

.. math::

   c_\text{Lsens} = \sqrt{s_\text{PSL}^2 + s_\text{AUX}^2}

where each :math:`s` is the transmission cost evaluated at
:math:`1.01 \times L`.

Esurf -- Surface Electric Field
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Penalizes high electric fields at the HR surface (which can cause
laser damage):

.. math::

   c_\text{Esurf} = 50 \cdot \text{arcsinh}(|1 + r|^2)

Lstdev -- Layer Thickness Uniformity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Penalizes deviations from a target mean-to-std ratio of layer
thicknesses:

.. math::

   c_\text{Lstdev} = \left|\frac{t - \mu_L / \sigma_L}{t}\right|^2

Absorption
^^^^^^^^^^

**Not implemented** in the optimizer. The underlying
:func:`~OptimalBragg.layers.calc_abs` function exists for
post-optimization analysis.

Optimizer
---------

The optimizer is ``scipy.optimize.differential_evolution`` with:

- Strategy: ``best1bin``
- Mutation: ``(0.05, 1.5)``
- Population: configurable via ``Nparticles``
- Workers: ``1`` (IPC overhead exceeds per-eval cost)
- Polish: ``True`` (L-BFGS-B refinement after convergence)

Layer thicknesses are bounded to ``[0.05, 0.48]`` in units of the design
wavelength.
