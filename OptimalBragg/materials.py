# Compilation of thin film and substrate material properties.
#
# Property key definitions:
# -----------------------------------------------------------------
# Name                   : Material name (string)
# Index                  : Refractive index at operating wavelength
# Y                      : Young's modulus [N/m^2]
# Sigma                  : Poisson ratio (dimensionless)
# Phi                    : Mechanical loss angle [rad]
# Alpha                  : Coefficient of thermal expansion [1/K]
# Beta                   : Thermo-refractive coefficient (dn/dT) [1/K]
# CV                     : Volumetric heat capacity [J/m^3/K]
# ThermalDiffusivity     : Thermal conductivity [W/m/K]  (name is legacy misnomer)
# Absorption             : Optical absorption [1/m]
# Temp                   : Reference temperature [K]
#
# Substrate-specific additional keys:
# rho / MassDensity      : Mass density [kg/m^3]
# MassAlpha              : Substrate CTE [1/K]
# MassCM                 : Specific heat [J/kg/K]
# MassKappa              : Substrate thermal conductivity [W/m/K]
# MirrorY                : Substrate Young's modulus [N/m^2]
# MirrorSigma            : Substrate Poisson ratio
# MechanicalLossExponent : Exponent in c2 * f^exponent
# c2                     : Prefactor for structural loss
# Alphas                 : Surface loss limit
# -----------------------------------------------------------------
#
# DISCLAIMER: There is no such thing as a comprehensive list of thin film
# properties for all wavelengths, temperatures, thicknesses, etc. The
# numbers below are a starting point for coating design — not definitive.
#
# TODO:
#     o Fill room-temp dicts for aSi, SiN, GaAs, AlGaAs
#     o Use Sellmeier coefficients for all thin films to cover any wavelength
#       >> Mine data from refractiveindex.org if not available?


# =====================================================================
#  Room temperature coatings @ 1064 nm, 295 K
# =====================================================================

air = {"Properties": {"Index": 1.00, "Name": "air"}}

SiO2 = {
    "Properties": {
        "Name": "SiO2",                   # See [0]
        "Y": 60e9,                         # See [0]
        "Sigma": 0.17,                     # See [0]
        "CV": 1.6412e6,                    # See [0], volumetric [J/m^3/K]
        "rho": 2203,                       # See [1]
        "CP": 744.98,                      # CV/rho [J/kg/K]
        "Alpha": 0.51e-6,                  # See [0]
        "MechanicalLossExponent": 1,       # c2 * f^exponent
        "c2": 3e-13,
        "ThermalDiffusivity": 2.0,         # Thermal conductivity [W/m/K], See [0]
        "Kappa": 0.5,                      # 1/ThermalDiffusivity (legacy)
        "Beta": 8e-6,                      # See [0]
        "Phi": 0.5e-4,                     # See [0]
        "Index": 1.45,                     # See [0]
        "Absorption": 0,                   # See [3]; negligible at 1064 nm
        "Temp": 295,                       # aLIGO operating temperature
    },
    "References": {
        0: "https://arxiv.org/pdf/0912.0107.pdf",
        1: "https://srdata.nist.gov/CeramicDataPortal/Elasticity/SiO2",
        2: "https://doi.org/10.1364/AO.51.006789",
        3: "https://iopscience.iop.org/article/10.1088/1361-6382/ab77e9",
    },
}

TiTa2O5 = {
    "Properties": {
        "Name": "TiTa2O5",                # See [0]
        "Y": 140e9,                        # See [0]
        "Sigma": 0.23,                     # See [0]
        "CV": 1.7283e6,                    # See [0], volumetric [J/m^3/K]
        "Alpha": 3.6e-6,                   # See [0]
        "ThermalDiffusivity": 1.67,        # Thermal conductivity [W/m/K], See [0]
        "Beta": 14e-6,                     # See [0]
        "Phi": 2e-4,                       # See [0]
        "Index": 2.06,                     # See [0]
        "Absorption": 100,                 # [1/m]
        "Temp": 295,                       # aLIGO operating temperature
    },
    "References": {
        0: "https://arxiv.org/pdf/0912.0107.pdf",
    },
}

Ta2O5 = {
    "Properties": {
        "Name": "Ta2O5",                   # See [0]
        "Y": 140e9,                        # See [0]
        "Sigma": 0.23,                     # See [0]
        "CV": 2.0961e6,                    # See [0], volumetric [J/m^3/K]
        "Alpha": 3.6e-6,                   # See [0]
        "ThermalDiffusivity": 1.67,        # Thermal conductivity [W/m/K], See [0]
        "Beta": 2.3e-6,                    # See [0]
        "Phi": 3.8e-4,                     # See [0]
        "Index": 2.06,                     # See [0], [2]
        "Absorption": 20e-6 / 500e-9,      # ~40 [1/m], See [1]
        "Temp": 295,
    },
    "References": {
        0: "https://arxiv.org/pdf/0912.0107.pdf",
        1: "https://opg.optica.org/abstract.cfm?uri=OIC-2019-FA.6",
        2: "https://doi.org/10.1063/1.4819325",
    },
}


# =====================================================================
#  Fused silica substrate — room temperature, 295 K
#  Source: aLIGO gwinc YAML (aLIGO_SiO2Ta2O5.yaml)
# =====================================================================

FusedSilica = {
    "Properties": {
        "Name": "FusedSilica",
        "Index": 1.45,                     # Substrate refractive index
        "Y": 72.7e9,                       # = MirrorY
        "Sigma": 0.167,                    # = MirrorSigma
        "Phi": 5e-5,                       # Bulk fused silica loss at ~100 Hz
        "Alpha": 3.9e-7,                   # = MassAlpha, CTE [1/K]
        "Beta": 8e-6,                      # dn/dT for fused silica at 1064 nm
        "CV": 1.6258e6,                    # MassCM * MassDensity = 739 * 2200
        "ThermalDiffusivity": 1.38,        # = MassKappa, thermal conductivity [W/m/K]
        "Absorption": 0,                   # Substrate bulk absorption handled separately
        "Temp": 295,
        # Aliases for OptimalBragg noise functions (sub.CP, sub.Kappa)
        "CP": 739,                         # = MassCM, specific heat [J/kg/K]
        "Kappa": 1.38,                     # = MassKappa, thermal conductivity [W/m/K]
        # Substrate-specific properties (from gwinc YAML)
        "rho": 2200,                       # = MassDensity [kg/m^3]
        "MassDensity": 2200,
        "MassAlpha": 3.9e-7,               # CTE [1/K]
        "MassCM": 739,                     # Specific heat [J/kg/K]
        "MassKappa": 1.38,                 # Thermal conductivity [W/m/K]
        "MirrorY": 72.7e9,                 # Young's modulus [N/m^2]
        "MirrorSigma": 0.167,              # Poisson ratio
        "MechanicalLossExponent": 0.77,
        "c2": 7.6e-12,
        "Alphas": 5.2e-12,                 # Surface loss limit
    },
    "References": {
        0: "aLIGO gwinc YAML (aLIGO_SiO2Ta2O5.yaml)",
        1: "https://arxiv.org/pdf/0912.0107.pdf",
    },
}


# =====================================================================
#  Cryogenic coatings @ 2050 nm, 123 K  (LIGO Voyager)
# =====================================================================

cSi_123 = {
    "Properties": {
        "Name": "c-Si",
        "Y": 155.8e9,                      # See [0]
        "Sigma": 0.27,                     # See [0]
        "CP": 300,                         # See [0], specific heat [J/kg/K]
        "rho": 2329,                       # See [0], [kg/m^3]
        "CV": 0.6987e6,                    # See [0], volumetric = rho * CP
        "Alpha": 1e-9,                     # See [0], nearly zero at 123 K
        "Index": 3.5,                      # See [0], 3.38*(1 + 4e-5*T)
        "Alphas": 5.2e-12,                 # Surface loss limit
        "MechanicalLossExponent": 1,       # c2 * f^exponent
        "Kappa": 700,                      # See [0], thermal conductivity [W/m/K]
        "ThermalDiffusivity": 0.00143,     # See [0], 1/Kappa
        "Beta": 1e-4,                      # See [1]
        "Phi": 3e-13,                      # At 1 Hz, See [2]
        "c2": 3e-13,                       # See [2]
        "Temp": 123,                       # Kelvin
        "Absorption": 0,                   # Substrate; bulk absorption handled separately
        # Substrate-specific properties (from Voyager gwinc YAML)
        "MassDensity": 2329,
        "MassAlpha": 1e-9,                 # CTE [1/K]
        "MassCM": 300,                     # Specific heat [J/kg/K]
        "MassKappa": 700,                  # Thermal conductivity [W/m/K]
        "MirrorY": 155.8e9,                # Young's modulus [N/m^2]
        "MirrorSigma": 0.27,               # Poisson ratio
    },
    "References": {
        0: "http://www.ioffe.ru/SVA/NSM/Semicond/Si/index.html",
        1: "http://arxiv.org/abs/physics/0606168",
        2: "https://doi.org/10.1016/0375-9601(81)90635-6",
    },
}

aSi_123 = {
    "Properties": {
        "Name": "a-Si",                    # See [0]
        "Y": 147e9,                        # See [1]
        "Sigma": 0.23,                     # See [2]
        "CV": 1.05e6,                      # See [3], volumetric [J/m^3/K]
        "Alpha": 1e-9,                     # Near zero at cryo
        "ThermalDiffusivity": 1.03,        # Thermal conductivity [W/m/K], See [3]
        "Beta": 1.4e-4,                    # See [4]
        "Phi": 2e-5,                       # See [5], See [8] (H:aSi)
        "Index": 3.65,                     # See [6]
        "Absorption": 20,                  # 27 See [8], 540 See [6], 20 See [7]
        "Temp": 123,
    },
    "References": {
        0: "https://en.wikipedia.org/wiki/Amorphous_silicon",
        1: "https://theses.gla.ac.uk/3671/, 5.5.5",
        2: "https://link.aps.org/doi/10.1103/PhysRevD.103.042001",
        3: "https://link.aps.org/doi/10.1103/PhysRevLett.96.055902",
        4: "http://dx.doi.org/10.1063/1.1383056",
        5: "https://journals.aps.org/prd/abstract/10.1103/PhysRevD.103.042001",
        6: "https://link.aps.org/doi/10.1103/PhysRevLett.120.263602",
        7: "Personal comm from Manel Ruiz to RXA 3/2023",
        8: "Personal comm from Manel Ruiz to PS & RXA 10/2023",
    },
}

SiN_123 = {
    "Properties": {
        "Name": "SiN",
        "Y": 270.0e9,                      # See [0, 1]
        "Sigma": 0.25,                     # See [0]
        "CV": 0.729e6,                     # Volumetric [J/m^3/K] = 230 J/kg/K * 3170 kg/m^3
                                            # See [2, 3]
        "Alpha": 2.6e-6,                   # See [4]
        "ThermalDiffusivity": 0.27,        # Thermal conductivity [W/m/K], See [2]
        "Phi": 0.8e-4,                     # See [0, 1]
        "Index": 2.17,                     # See [5]
        "Beta": 4e-5,                      # See [6]
        "Absorption": 546,                 # See [0]; 27 or 6 also reported
        "Temp": 123,
    },
    "References": {
        0: "https://link.aps.org/doi/10.1103/PhysRevD.98.102001",
        1: "https://link.aps.org/doi/10.1103/PhysRevD.96.022007",
        2: "2017 ECS J. Solid State Sci. Technol. 6 P691",
        3: "DOI: 10.1115/1.2945904",
        4: "https://doi.org/10.1364/AO.51.007229",
        5: "https://link.aps.org/doi/10.1103/PhysRevLett.96.055902",
        6: "10.1109/JPHOT.2016.2561622",
    },
}

SiO2_123 = {
    "Properties": {
        "Name": "silica",
        "Y": 72e9,
        "Sigma": 0.17,
        "CV": 0.744e6,                     # Volumetric [J/m^3/K]
        "Alpha": 0.0145e-6,
        "Beta": 4.2e-6,
        "ThermalDiffusivity": 1.05,        # Thermal conductivity [W/m/K], See [1]
        "Phi": 2e-4,
        "Index": 1.43545,                  # Calculated using [2]
        "Absorption": 245,                 # 295 See [3], 245 See [4]
        "Temp": 123,
    },
    "References": {
        0: "https://wiki.ligo.org/OPT/SilicaCoatingProp",
        1: "http://dx.doi.org/10.1109/ITHERM.2002.1012450",
        2: "http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=317500",
        3: "https://doi.org/10.1364/AO.51.006789",
        4: "https://refractiveindex.org",
    },
}

Ta2O5_123 = {
    "Properties": {
        "Name": "tantala",                 # See [0]
        "Y": 136e9,
        "Sigma": 0.22,
        "CP": 165.68,                      # [J/kg/K]
        "rho": 8180,                       # [kg/m^3]
        "CV": 1.355e6,                     # Volumetric = rho * CP [J/m^3/K]
        "Alpha": 0.09e-6,
        "Beta": 0.4e-6,                    # See [1]
        "ThermalDiffusivity": 1.03,        # Thermal conductivity [W/m/K]
        "Phi": 5e-4,                       # See [2]
        "Index": 2.083,                    # See [3]
        "Absorption": 35,                  # 27553 See [3]; optimistic 35 using [4, 5]
        "Temp": 123,
    },
    "References": {
        0: "https://doi.org/10.1364/AO.48.004536",
        1: "http://dx.doi.org/10.1063/1.1383056",
        2: "https://arxiv.org/pdf/1903.06094.pdf",
        3: "https://doi.org/10.1063/1.4819325",
        4: "https://opg.optica.org/abstract.cfm?uri=OIC-2019-FA.6",
        5: "https://refractiveindex.org",
    },
}
