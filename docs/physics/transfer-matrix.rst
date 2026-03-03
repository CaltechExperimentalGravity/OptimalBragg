Transfer Matrix Method
======================

The core physics computation is the transfer matrix method (TMM) for
multilayer dielectric stacks, implemented in
:func:`~OptimalBragg.layers.multidiel1`.

Theory
------

A dielectric stack consists of :math:`M` layers with refractive indices
:math:`n_1, n_2, \ldots, n_M` sandwiched between an incident medium
(:math:`n_a`, typically air) and a substrate (:math:`n_s`).

At each interface between media :math:`i` and :math:`i+1`, the Fresnel
reflection coefficient is:

.. math::

   r_i = \frac{n_i^T - n_{i+1}^T}{n_i^T + n_{i+1}^T}

where :math:`n^T` is the effective index accounting for polarization:

- **TE (s-pol):** :math:`n^T = n \cos\theta`
- **TM (p-pol):** :math:`n^T = n / \cos\theta`

The total reflection coefficient is computed by recursive application
starting from the substrate:

.. math::

   \Gamma_i = \frac{r_i + \Gamma_{i+1} e^{-2j\delta_i}}{1 + r_i \Gamma_{i+1} e^{-2j\delta_i}}

where :math:`\delta_i = 2\pi L_i / \lambda` is the phase accumulated
in layer :math:`i`, with :math:`L_i` being the optical thickness
normalized to the design wavelength.

Implementation
--------------

The inner loop of this recursion is the computational bottleneck -- it
is called millions of times during a single optimization run. We use
Numba JIT compilation (``@numba.njit(cache=True)``) to achieve ~7 us
per call for a 14-bilayer stack, compared to ~200-500 us in pure Python.

The function :func:`~OptimalBragg.layers.multidiel1` accepts an array
of wavelengths, enabling a single call to evaluate reflectivity at
multiple wavelengths simultaneously. This is exploited in
:func:`~OptimalBragg.costs.getMirrorCost` to consolidate what would
otherwise be 6 separate calls into just 2.

Quarter-Wave Stacks
-------------------

A standard quarter-wave (QW) stack has all layers with optical thickness
:math:`L = 0.25` (one quarter of the design wavelength). The power
reflectivity of a QW stack with :math:`N` bilayers is:

.. math::

   R = \left(\frac{1 - (n_H/n_L)^{2N} n_s}{1 + (n_H/n_L)^{2N} n_s}\right)^2

The optimizer searches for non-QW designs that trade reflectivity for
lower thermal noise, better auxiliary wavelength transmission, or other
objectives.
