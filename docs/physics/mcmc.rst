Monte Carlo Sensitivity Analysis
=================================

After finding an optimal coating design, we assess its robustness to
manufacturing uncertainties using a Monte Carlo (MC) analysis with the
``emcee`` ensemble sampler.

Perturbation Model
------------------

The MC perturbs four parameters simultaneously:

.. list-table::
   :header-rows: 1

   * - Parameter
     - Perturbation
     - Physical meaning
   * - Angle of incidence
     - 0.5% Gaussian
     - Alignment tolerance
   * - High-n refractive index
     - 0.5% Gaussian
     - Deposition variability
   * - Low-n refractive index
     - 0.5% Gaussian
     - Deposition variability
   * - Layer thickness (all)
     - 0.5% Gaussian
     - Thickness control

Each parameter is drawn from a 4-D multivariate Gaussian centered at
zero with covariance :math:`\Sigma = (0.005)^2 I_4`. The perturbation
is applied multiplicatively: :math:`x_\text{perturbed} = x_0 (1 + \epsilon)`.

emcee Configuration
-------------------

- **Sampler:** ``emcee.EnsembleSampler`` (affine-invariant ensemble)
- **Walkers:** 64
- **Burn-in:** 1000 steps (discarded via ``sampler.reset()``)
- **Production:** 5000 steps
- **Target distribution:** Multivariate Gaussian
  (``lnprob`` in :mod:`generic.coatingUtils`)

The sampler generates correlated draws from the perturbation
distribution. A subset of ``nSamples`` points from the flattened chain
is used for the MC evaluation.

Per-Sample Evaluation
---------------------

For each perturbed parameter set, the code:

1. Scales the refractive index arrays (low-n and high-n independently)
2. Scales all layer thicknesses uniformly
3. Recomputes :math:`T_\text{1064}` and :math:`T_\text{AUX}` via
   :func:`~generic.coatingUtils.multidiel1`
4. Computes the surface E-field via
   :func:`~generic.coatingUtils.surfaceField`
5. Evaluates the full thermal noise spectrum (Brownian + thermo-optic)
   via ``pygwinc``

Output
------

The MC output HDF5 file contains:

- ``MCout`` — shape ``(5, nSamples)``:
  rows are :math:`T_\text{1064}` [ppm], :math:`T_\text{AUX}` [%],
  thermo-optic noise at 100 Hz [:math:`\times 10^{-21}`],
  Brownian noise at 100 Hz [:math:`\times 10^{-21}`],
  surface E-field [V/m]
- ``TOnoise`` — full TO noise spectra, shape ``(nSamples+1, len(ff))``
- ``Brnoise`` — full Brownian spectra, shape ``(nSamples+1, len(ff))``

Interpreting Corner Plots
-------------------------

The ``cornerPlt.py`` script generates corner plots (pair plots) from
the MC output using the ``corner`` package. Each panel shows:

- **Diagonal:** Marginalized 1-D histogram for one observable
- **Off-diagonal:** 2-D joint distribution between two observables

Key things to look for:

- **T_1064 spread:** If the 1-sigma width exceeds the target
  (e.g., 5 ppm), the design is sensitive to manufacturing errors.
- **T_AUX vs T_1064 correlation:** Strong correlation means the two
  wavelengths respond similarly to perturbations (typical for
  quarter-wave-like designs).
- **Thermal noise stability:** A tight distribution means the noise
  performance is robust.

Usage
-----

::

    python doMC.py Data/ETM/ETM_Layers_YYMMDD_HHMMSS.hdf5 output_MC.hdf5 5000
    python cornerPlt.py output_MC.hdf5
