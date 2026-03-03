Thermal Noise
=============

Coating thermal noise is the dominant noise source in the detection band
of ground-based gravitational wave detectors. This code computes three
components: Brownian, thermo-elastic, and thermo-refractive (the latter
two together form "thermo-optic" noise).

Coating Brownian Noise
----------------------

Brownian noise arises from mechanical dissipation in the coating
materials, quantified by their loss angles :math:`\phi`. The power
spectral density of displacement noise is (Harry et al., CQG 19, 897,
2002):

.. math::

   S_\text{Br}(f) = \frac{2 k_B T}{\pi^{3/2} f}
     \frac{1 - \sigma_s^2}{Y_s w}
     \left( d_L \frac{Y_s}{Y_L} \phi_L + d_H \frac{Y_s}{Y_H} \phi_H \right)

where :math:`d_L, d_H` are the total physical thicknesses of low-n and
high-n layers, :math:`Y` are Young's moduli, :math:`\phi` are loss
angles, and :math:`w` is the beam radius.

The optimizer uses a **fast proxy** (LIGO-E0900068 p. 4) that avoids
recomputing all material constants on every evaluation:

.. math::

   c_\text{Brownian} = t \cdot (z_L + \gamma \, z_H)

where :math:`z_L, z_H` are the total *optical* thicknesses and
:math:`\gamma` encodes the material property ratios (pre-computed once
by :func:`~OptimalBragg.noise.brownian_proxy`).

Thermo-Optic Noise
-------------------

Thermo-optic noise combines two mechanisms driven by statistical
temperature fluctuations in the coating (Evans et al., PRD 78, 102003,
2008; LIGO-T080101):

**Thermo-elastic (TE):** Thermal expansion changes the coating
thickness, shifting the reflected phase:

.. math::

   \frac{d\phi}{dT}\bigg|_\text{TE} = \frac{4\pi}{\lambda}
     \sum_i \bar\alpha_i \, d_i

where :math:`\bar\alpha_i` is the effective expansion coefficient
(corrected for Poisson's ratio and substrate constraint).

**Thermo-refractive (TR):** Temperature changes the refractive indices
via :math:`dn/dT` (the :math:`\beta` coefficient):

.. math::

   \frac{d\phi}{dT}\bigg|_\text{TR} =
     \sum_i \frac{d\phi}{d\delta_i} (\beta_i + \bar\sigma_i n_i) \, d_i

The total thermo-optic noise PSD is:

.. math::

   S_\text{TO}(f) = S_\text{surf}(f) \cdot g_\text{TO}(\xi)
     \cdot \left(\frac{d\phi_\text{TE}}{dT} + \frac{d\phi_\text{TR}}{dT}\right)^2

where :math:`S_\text{surf}` is the surface temperature fluctuation
spectrum and :math:`g_\text{TO}(\xi)` is the finite coating thickness
correction (Braginsky & Vyatchanin, PLA 312, 244, 2003).

Implementation
--------------

The canonical implementation lives in ``pygwinc``
(:func:`gwinc.noise.coatingthermal.coating_thermooptic`). For the
optimizer hot path, we use a JIT-compiled version in
:mod:`OptimalBragg.noise` that inlines all sub-calculations into a
single Numba-compiled function. This gives identical results with ~10x
speedup.

The JIT function :func:`~OptimalBragg.noise.coating_thermooptic_fast`
takes a pre-extracted ``stack_params`` tuple of scalars, created once
per optimization.

Finite Mirror Corrections
^^^^^^^^^^^^^^^^^^^^^^^^^

The code includes the Bondu-Hello-Vinet finite mirror size correction,
which accounts for the mirror being a finite cylinder rather than an
infinite half-space. This uses a Bessel function expansion (300 terms)
to compute the correction factor :math:`C_\text{fsm}`.

Material Properties
-------------------

Material properties are stored in the gwinc YAML structure files:

- **aLIGO** (``aLIGO_SiO2Ta2O5.yaml``): SiO2 (low-n) / Ti:Ta2O5
  (high-n) at 294 K. Ti-doped tantala has lower mechanical loss than
  pure Ta2O5 (:math:`\phi \approx 2.3 \times 10^{-4}` vs
  :math:`4 \times 10^{-4}`).
- **Voyager** (``aSiSiN.yaml``): SiN (low-n) / a-Si (high-n) at 123 K
  on silicon substrates. Cryogenic operation reduces thermal noise by
  :math:`\sqrt{T}`.

Key properties per material: refractive index :math:`n`, Young's
modulus :math:`Y`, Poisson ratio :math:`\sigma`, loss angle :math:`\phi`,
thermal expansion :math:`\alpha`, thermo-refractive coefficient
:math:`\beta`, heat capacity :math:`C_V`, thermal diffusivity
:math:`\kappa`.
