Material Properties
===================

Central library of thin-film and substrate material properties used by the
optimizer, noise models, and Monte Carlo analysis.  Source code:
:mod:`OptimalBragg.materials`.

.. note::

   Values marked **[?]** have no literature reference and need citation.
   The ``ThermalDiffusivity`` key is a legacy misnomer --- it stores
   *thermal conductivity* :math:`\kappa` in W/m/K.

Room-Temperature Materials (295 K, 1064 nm)
--------------------------------------------

Thin-film coatings
^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 22 12 22 22 22

   * - Property
     - Unit
     - SiO2
     - TiO2:Ta2O5
     - Ta2O5
   * - Refractive index (*n*)
     - ---
     - 1.45 [a0]_
     - 2.06 [b0]_
     - 2.06 [c0]_ [c2]_
   * - Young's modulus (*Y*)
     - GPa
     - 60 [a0]_
     - 140 [b0]_
     - 140 [c0]_
   * - Poisson ratio (:math:`\sigma`)
     - ---
     - 0.17 [a0]_
     - 0.23 [b0]_
     - 0.23 [c0]_
   * - Heat capacity (*C_V*)
     - 10\ :sup:`6` J/m\ :sup:`3`/K
     - 1.641 [a0]_
     - 1.728 [b0]_
     - 2.096 [c0]_
   * - CTE (:math:`\alpha`)
     - 10\ :sup:`-6` /K
     - 0.51 [a0]_
     - 3.6 [b0]_
     - 3.6 [c0]_
   * - Thermal conductivity (:math:`\kappa`)
     - W/m/K
     - 2.0 [a0]_
     - 1.67 [b0]_
     - 1.67 [c0]_
   * - Thermo-refractive (:math:`\beta`)
     - 10\ :sup:`-6` /K
     - 8.0 [a0]_
     - 14.0 [b0]_
     - 2.3 [c0]_
   * - Loss angle (:math:`\phi`)
     - 10\ :sup:`-4` rad
     - 0.5 [a0]_
     - 2.0 [b0]_
     - 3.8 [c0]_
   * - Absorption
     - 1/m
     - 0 [a3]_
     - 100 **[?]**
     - 40 [c1]_
   * - Density (:math:`\rho`)
     - kg/m\ :sup:`3`
     - 2203 [a1]_
     - ---
     - ---

Fused silica substrate
^^^^^^^^^^^^^^^^^^^^^^

Source: aLIGO gwinc YAML [fs0]_.

.. list-table::
   :header-rows: 1
   :widths: 35 15 30

   * - Property
     - Unit
     - Value
   * - Refractive index (*n*)
     - ---
     - 1.45
   * - Young's modulus (*Y*)
     - GPa
     - 72.7
   * - Poisson ratio (:math:`\sigma`)
     - ---
     - 0.167
   * - Heat capacity (*C_V*)
     - 10\ :sup:`6` J/m\ :sup:`3`/K
     - 1.626
   * - CTE (:math:`\alpha`)
     - 10\ :sup:`-7` /K
     - 3.9
   * - Thermal conductivity (:math:`\kappa`)
     - W/m/K
     - 1.38
   * - Thermo-refractive (:math:`\beta`)
     - 10\ :sup:`-6` /K
     - 8.0
   * - Loss angle (:math:`\phi`)
     - 10\ :sup:`-5` rad
     - 5.0
   * - Density (:math:`\rho`)
     - kg/m\ :sup:`3`
     - 2200
   * - Specific heat (*C_P*)
     - J/kg/K
     - 739
   * - Mech. loss exponent
     - ---
     - 0.77 **[?]**
   * - Structural loss *c2*
     - ---
     - 7.6e-12 **[?]**
   * - Surface loss (*alpha_s*)
     - ---
     - 5.2e-12 **[?]**

Cryogenic Materials (123 K, 2050 nm)
-------------------------------------

Thin-film coatings
^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 10 14 14 14 14 14

   * - Property
     - Unit
     - a-Si
     - SiN
     - SiO2
     - Ta2O5
     - c-Si (sub.)
   * - Refractive index (*n*)
     - ---
     - 3.65 [d6]_
     - 2.17 [e5]_
     - 1.435 [f2]_
     - 2.083 [g3]_
     - 3.50 [h0]_
   * - Young's modulus (*Y*)
     - GPa
     - 147 [d1]_
     - 270 [e0]_ [e1]_
     - 72 **[?]**
     - 136 **[?]**
     - 155.8 [h0]_
   * - Poisson ratio (:math:`\sigma`)
     - ---
     - 0.23 [d2]_
     - 0.25 [e0]_
     - 0.17 **[?]**
     - 0.22 **[?]**
     - 0.27 [h0]_
   * - Heat capacity (*C_V*)
     - 10\ :sup:`6` J/m\ :sup:`3`/K
     - 1.050 [d3]_
     - 0.729 [e2]_ [e3]_
     - 0.744 **[?]**
     - 1.355 **[?]**
     - 0.699 [h0]_
   * - CTE (:math:`\alpha`)
     - 10\ :sup:`-6` /K
     - ~0 **[?]**
     - 2.6 [e4]_
     - 0.0145 **[?]**
     - 0.09 **[?]**
     - ~0 [h0]_
   * - Thermal conductivity (:math:`\kappa`)
     - W/m/K
     - 1.03 [d3]_
     - 0.27 [e2]_
     - 1.05 [f1]_
     - 1.03 **[?]**
     - 700 [h0]_
   * - Thermo-refractive (:math:`\beta`)
     - 10\ :sup:`-4` /K
     - 1.4 [d4]_
     - 0.4 [e6]_
     - 0.042 **[?]**
     - 0.004 [g1]_
     - 1.0 [h1]_
   * - Loss angle (:math:`\phi`)
     - rad
     - 2e-5 [d5]_ [d8]_
     - 8e-5 [e0]_ [e1]_
     - 2e-4 **[?]**
     - 5e-4 [g2]_
     - 3e-13 [h2]_
   * - Absorption
     - 1/m
     - 20 [d7]_
     - 546 [e0]_
     - 245 [f4]_
     - 35 [g4]_ [g5]_
     - 0

Crystalline silicon substrate (c-Si, 123 K)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 35 15 30

   * - Property
     - Unit
     - Value
   * - Density (:math:`\rho`)
     - kg/m\ :sup:`3`
     - 2329 [h0]_
   * - Specific heat (*C_P*)
     - J/kg/K
     - 300 [h0]_
   * - Mech. loss exponent
     - ---
     - 1 **[?]**
   * - Structural loss *c2*
     - ---
     - 3e-13 [h2]_
   * - Surface loss (*alpha_s*)
     - ---
     - 5.2e-12 **[?]**

Project Overrides
-----------------

Per-project ``materials.yml`` files can override central library values.
Currently only mirror geometry parameters are overridden (not material physics).

.. list-table::
   :header-rows: 1
   :widths: 25 20 20 15

   * - Project
     - Material
     - Property
     - Value
   * - aLIGO
     - FusedSilica
     - MassRadius
     - 0.17 m
   * - aLIGO
     - FusedSilica
     - MassThickness
     - 0.20 m
   * - Voyager (aSi/SiN)
     - c-Si
     - MassRadius
     - 0.215 m
   * - Voyager (aSi/SiN)
     - c-Si
     - MassThickness
     - 0.55 m
   * - Voyager (Ta2O5)
     - c-Si
     - MassRadius
     - 0.215 m
   * - Voyager (Ta2O5)
     - c-Si
     - MassThickness
     - 0.55 m

Placeholder Materials
---------------------

The following materials are declared in the library but have no properties
populated yet:

- ``aSi`` --- amorphous silicon (room temperature)
- ``SiN`` --- silicon nitride (room temperature)
- ``GaAs`` --- gallium arsenide
- ``AlGaAs`` --- aluminum gallium arsenide
- ``H2O_123`` --- water at 123 K

References
----------

**SiO2 (295 K)**

.. [a0] Evans *et al.*, "Thermo-optic noise in coated mirrors for
   high-precision optical measurements," `arXiv:0912.0107
   <https://arxiv.org/pdf/0912.0107.pdf>`_
.. [a1] NIST Ceramics Data Portal, `SiO2 Elasticity
   <https://srdata.nist.gov/CeramicDataPortal/Elasticity/SiO2>`_
.. [a3] Steinlechner *et al.*, Class. Quantum Grav. **37** 095004 (2020),
   `doi:10.1088/1361-6382/ab77e9
   <https://iopscience.iop.org/article/10.1088/1361-6382/ab77e9>`_

**TiO2:Ta2O5 (295 K)**

.. [b0] Evans *et al.*, `arXiv:0912.0107
   <https://arxiv.org/pdf/0912.0107.pdf>`_ (same as [a0]_)

**Ta2O5 (295 K)**

.. [c0] Evans *et al.*, `arXiv:0912.0107
   <https://arxiv.org/pdf/0912.0107.pdf>`_ (same as [a0]_)
.. [c1] Yang *et al.*, OIC 2019 FA.6,
   `doi:10.1364/OIC.2019.FA.6
   <https://opg.optica.org/abstract.cfm?uri=OIC-2019-FA.6>`_
.. [c2] Cetinorgu-Goldenberg *et al.*, J. Appl. Phys. **113** 143511 (2013),
   `doi:10.1063/1.4819325 <https://doi.org/10.1063/1.4819325>`_

**Fused silica substrate**

.. [fs0] aLIGO gwinc YAML (``aLIGO_SiO2Ta2O5.yaml``)

**a-Si (123 K)**

.. [d1] Li, PhD thesis, Univ. Glasgow (2012), sec. 5.5.5,
   `link <https://theses.gla.ac.uk/3671/>`_
.. [d2] Vajente *et al.*, PRD **103** 042001 (2021),
   `doi:10.1103/PhysRevD.103.042001
   <https://link.aps.org/doi/10.1103/PhysRevD.103.042001>`_
.. [d3] Zink *et al.*, PRL **96** 055902 (2006),
   `doi:10.1103/PhysRevLett.96.055902
   <https://link.aps.org/doi/10.1103/PhysRevLett.96.055902>`_
.. [d4] Inci & Yoshino, J. Appl. Phys. **89** 6956 (2001),
   `doi:10.1063/1.1383056 <http://dx.doi.org/10.1063/1.1383056>`_
.. [d5] Vajente *et al.*, PRD **103** 042001 (2021),
   `doi:10.1103/PhysRevD.103.042001
   <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.103.042001>`_
.. [d6] Birney *et al.*, PRL **120** 263602 (2018),
   `doi:10.1103/PhysRevLett.120.263602
   <https://link.aps.org/doi/10.1103/PhysRevLett.120.263602>`_
.. [d7] Personal comm. from M. Ruiz to RXA (March 2023)
.. [d8] Personal comm. from M. Ruiz to PS & RXA (October 2023)

**SiN (123 K)**

.. [e0] Prasai *et al.*, PRD **98** 102001 (2018),
   `doi:10.1103/PhysRevD.98.102001
   <https://link.aps.org/doi/10.1103/PhysRevD.98.102001>`_
.. [e1] Martin *et al.*, PRD **96** 022007 (2017),
   `doi:10.1103/PhysRevD.96.022007
   <https://link.aps.org/doi/10.1103/PhysRevD.96.022007>`_
.. [e2] ECS J. Solid State Sci. Technol. **6** P691 (2017)
.. [e3] Okhotin *et al.*, J. Heat Transfer **130** 082403 (2008),
   DOI: 10.1115/1.2945904
.. [e4] Inci & Yoshino, Appl. Opt. **51** 7229 (2012),
   `doi:10.1364/AO.51.007229 <https://doi.org/10.1364/AO.51.007229>`_
.. [e5] Zink *et al.*, PRL **96** 055902 (2006),
   `doi:10.1103/PhysRevLett.96.055902
   <https://link.aps.org/doi/10.1103/PhysRevLett.96.055902>`_
.. [e6] IEEE Photon. J. **8** 2500109 (2016),
   DOI: 10.1109/JPHOT.2016.2561622

**SiO2 (123 K)**

.. [f1] Luo *et al.*, Proc. ITHERM 2002,
   `doi:10.1109/ITHERM.2002.1012450
   <http://dx.doi.org/10.1109/ITHERM.2002.1012450>`_
.. [f2] Malitson, J. Opt. Soc. Am. **55** 1205 (1965),
   `IEEE <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=317500>`_
.. [f4] `refractiveindex.info <https://refractiveindex.org>`_

**Ta2O5 (123 K)**

.. [g1] Inci & Yoshino, J. Appl. Phys. **89** 6956 (2001),
   `doi:10.1063/1.1383056 <http://dx.doi.org/10.1063/1.1383056>`_
.. [g2] Granata *et al.* (2019), `arXiv:1903.06094
   <https://arxiv.org/pdf/1903.06094.pdf>`_
.. [g3] Cetinorgu-Goldenberg *et al.*, J. Appl. Phys. **113** 143511 (2013),
   `doi:10.1063/1.4819325 <https://doi.org/10.1063/1.4819325>`_
.. [g4] Yang *et al.*, OIC 2019 FA.6,
   `doi:10.1364/OIC.2019.FA.6
   <https://opg.optica.org/abstract.cfm?uri=OIC-2019-FA.6>`_
.. [g5] `refractiveindex.info <https://refractiveindex.org>`_

**Crystalline silicon substrate (123 K)**

.. [h0] Ioffe Institute NSM, `Silicon
   <http://www.ioffe.ru/SVA/NSM/Semicond/Si/index.html>`_
.. [h1] Braginsky *et al.* (2006), `arXiv:physics/0606168
   <http://arxiv.org/abs/physics/0606168>`_
.. [h2] Braginsky & Mitrofanov, Phys. Lett. A **82** 435 (1981),
   `doi:10.1016/0375-9601(81)90635-6
   <https://doi.org/10.1016/0375-9601(81)90635-6>`_
