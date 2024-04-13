# Compilation of thin film material properties including:
# Y                   : Young's modulus in [N / m**2]
# Sigma               : Poisson ratio in [kg / m**3]
# CP                  : Constant pressure specific heat in [J / kg / K]
# rho                 : Volume density in [kg / m**3]
# CV                  : Constant volume specific heat in [J / kg / K]
# Alpha               : Coefficient of thermal expansion in [1 / K]
# ThermalDiffusivity  : (inverse thermal conductivity) in [K m / W]
# Phi                 : Mechanical loss angle in [rad]
# Index               : Refractive index
# Absorption          : Absorption coefficient in [1 / m]
# Beta                : Thermorefractive coefficient in [1 / K]
# Kappa               : Thermal conductivity in [W / m / K]

""" 
TODO: 
    o Fill room-temp dicts using all/either of:
    >> https://arxiv.org/pdf/0912.0107.pdf
    >> https://wiki.ligo.org/OPT/CoatingProperties
    >> https://git.ligo.org/voyager/mariner40/-/wikis/optics/dichroic-coatings

WISHLIST:
    o Use sellmeier coefficients for all thin films to cover any wavelength
    >> Mine data from refractiveindex.org if not available?
"""
aSi = {}
SiN = {}
GaAs = {}
AlGaAs = {}

# Room temperature @ 1.06 um

air = {"Properties": {"Index": 1.00, "Name": "air"}}

SiO2 = {
    "Properties": {
        "Name": "SiO2",  # See [0]
        "Y": 60e9,  # See [0]
        "Sigma": 0.17,  # See [0]
        "CV": 1.6412e6,  # See [0]
        "rho": 2203,  # See [1]
        "CP": 744.98,  # CV/rho
        "Alpha": 0.51e-6,  # See [0]
        "MechanicalLossExponent": 1,  # c2 * f^exponent
        "c2": 3e-13,
        "ThermalDiffusivity": 2.0,  # See [0]
        "Kappa": 0.5,  # 1/thermaldiffusivity
        "Beta": 8e-6,  # See [0]
        "Phi": 0.5e-4,  # See [0]
        "Index": 1.45,  # - 0.8j,  # See [0]
        "Absorption": 0,  #  See[3], or 162.14 See [2]
        "Temp": 300,
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
        "Name": "TiTa2O5",  # See [0]
        "Y": 140e9,  # See [0]
        "Sigma": 0.23,  # See [0]
        "CV": 1.7283e6,  # See [0]
        "Alpha": 3.6e-6,  # See [0]
        "ThermalDiffusivity": 1.67,  # See [0]
        "Beta": 14e-6,  # See [0]
        "Phi": 2e-4,  # See [0]
        "Index": 2.06,  # See [0]
        "Absorption": 100,  # ??
    },
    "References": {0: "https://arxiv.org/pdf/0912.0107.pdf"},
}
Ta2O5 = {
    "Properties": {
        "Name": "Ta2O5",  # See [0]
        "Y": 140e9,  # See [0]
        "Sigma": 0.23,  # See [0]
        "CV": 2.0961e6,  # See [0]
        "Alpha": 3.6e-6,  # See [0]
        "ThermalDiffusivity": 1.67,  # See [0]
        "Beta": 2.3e-6,  # See [0]
        "Phi": 3.8e-4,  # See [0]
        "Index": 2.06,  # - 0.0017775j  # See [0], [2]
        "Absorption": 20e-6 / 500e-9,  # See [1]
    },
    "References": {
        0: "https://arxiv.org/pdf/0912.0107.pdf",
        1: "https://opg.optica.org/abstract.cfm?uri=OIC-2019-FA.6",
        2: "https://doi.org/10.1063/1.4819325",
    },
}

# Cryo-Si (123 K) @ 2 um
cSi_123 = {
    "Properties": {
        "Name": "c-Si",
        "Y": 155.8e9,  # See [0]
        "Sigma": 0.27,  # See [0]
        "CP": 300,  # See [0]
        "rho": 2329,  # See [0]
        "CV": 0.6987e6,  # See [0]
        "Alpha": 1e-9,  # See [0]
        "Index": 3.5,  # See [0] 3.38 * (1 + 4e-5 * T)
        "Alphas": 5.2e-12,  # Surface loss limit??
        "MechanicalLossExponent": 1,  # c2 * f^exponent
        "Kappa": 700,  # See [0]
        "ThermalDiffusivity": 0.00143,  # See [0]
        "Beta": 1e-4,  # See [1]
        "Phi": 3e-13,  # At 1 Hz, See [2]
        "c2": 3e-13,  # See [2]
        "Temp": 123,  # Kelvin
    },
    "References": {
        0: "http://www.ioffe.ru/SVA/NSM/Semicond/Si/index.html",
        1: "http://arxiv.org/abs/physics/0606168",
        2: "https://doi.org/10.1016/0375-9601(81)90635-6",
    },
}

aSi_123 = {
    "Properties": {
        "Name": "a-Si",  # See [0]
        "Y": 147e9,  # See [1]
        "Sigma": 0.23,  # See [2]
        "CV": 1.05e6,  # See [3]
        "Alpha": 1e-9,
        "ThermalDiffusivity": 1.03,  # See [3]
        "Beta": 1.4e-4,  # See [4]
        "Phi": 2e-5,  # See [5], See [8] (H:aSi)
        "Index": 3.65,  # See [6]
        "Absorption": 20,  # 27 See [8] (H:aSi), 540.35 see [6], also 10e-6/0.5e-6 see [7]
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
        "Y": 270.0e9,  # 103.7e9, See [0,1]
        "Sigma": 0.25,  # See [0]
        "CV": 230,  # See [2, 3]
        "Alpha": 2.6e-6,  # See [4]
        "ThermalDiffusivity": 0.27,  # See [2]
        "Phi": 0.8e-4,  # See [0, 1]
        "Index": 2.17,  # See [5]
        "Beta": 4e-5,  # See [6]
        "Absorption": 546,  # See [0], 27? 6?
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
        "CV": 0.744e6,
        "Alpha": 0.0145e-6,
        "Beta": 4.2e-6,
        "ThermalDiffusivity": 1.05,  # See [1]
        "Phi": 2e-4,
        "Index": 1.43545,  # Calculated using [2]
        "Absorption": 245,  # 295 See [3], 245 see [4]
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
        "Name": "tantala",  # See [0]
        "Y": 136e9,
        "Sigma": 0.22,
        "CP": 165.68,
        "rho": 8180,
        "CV": 1.355e6,
        "Alpha": 0.09e-6,
        "Beta": 0.4e-6,  # See [1]
        "ThermalDiffusivity": 1.03,
        "Phi": 5e-4,  # See [2]
        "Index": 2.083,  # See [3]
        "Absorption": 35,  # 27553 See [3], or optimistic 35 using [4-5]
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


H2O_123 = {}


def n_SiO2(λ):
    """Reference https://slac.stanford.edu/grp/arb/tn/arbvol5/ARDB436.pdf"""

    # Input wavelength λ needs to be in microns
    def sellmeier(
        A=1,
        Bs=[
            0,
        ],
        Cs=[
            0,
        ],
    ):
        n_squared = A
        for B, C in zip(Bs, Cs):
            n_squared += B / (λ**2 - C)
        return n_squared

    pass
    # return sqrt(
    # sellmeier(1.3107237, [0.7935797, 0.9237144], [1.0959659e-2, 1080])
    # )
