"""OptimalBragg — optical coating design and optimization.

Provides Material handling, quarter-wave stack construction,
and I/O utilities for HDF5 and YAML formats.
"""

import numpy as np

from OptimalBragg.io import h5read, h5write, yamlread, load_materials_yaml


class Material:
    """Thin-film or substrate material with properties as attributes.

    Parameters
    ----------
    props : dict
        Must contain a ``"Properties"`` key whose value is a dict of
        material properties (Index, Y, Sigma, Phi, etc.).  Each
        key-value pair is set as an instance attribute.
    """

    def __init__(self, props):
        self.__dict__.update(props["Properties"])


def qw_stack(lam_ref, substrate, superstrate, thin_films, pattern, hwcap=""):
    """Build a stack attribute dict from a quarter-wave design at a single wavelength.

    Parameters
    ----------
    lam_ref : float
        Reference wavelength [m].
    substrate : Material
        Substrate material (e.g. fused silica, c-Si).
    superstrate : Material
        Superstrate material (e.g. air).
    thin_films : dict
        Mapping of pattern characters to Material objects,
        e.g. ``{"H": high_index_mat, "L": low_index_mat}``.
    pattern : str
        Layer structure from superstrate side, e.g. ``"HLHLHLHL"``.
    hwcap : str, optional
        Half-wave cap layers at the superstrate interface.  Each
        character maps to a key in *thin_films*.  Default is no cap.

    Returns
    -------
    dict
        Stack attributes including refractive indices, physical and
        optical thicknesses, material properties, and metadata.
    """
    # Build refractive-index and physical-thickness lists
    n_stack = [superstrate.Index]
    L_stack = []

    # Half-wave cap layers (optical thickness = 0.5)
    if hwcap:
        n_stack += [thin_films[X].Index for X in hwcap]
        L_stack = [(lam_ref / 2) / nj for nj in n_stack[1:]]

    # Quarter-wave layers (optical thickness = 0.25)
    n_stack += [thin_films[X].Index for X in pattern]
    L_stack += [(lam_ref / 4) / nj for nj in n_stack[len(hwcap) + 1:]]

    # Complete pattern including any cap
    pattern = hwcap + pattern

    # Gather material properties for each layer
    a_stack = [thin_films[X].Absorption for X in pattern]
    Y_stack = [thin_films[X].Y for X in pattern]
    sigma_stack = [thin_films[X].Sigma for X in pattern]
    phi_stack = [thin_films[X].Phi for X in pattern]
    cte_stack = [thin_films[X].Alpha for X in pattern]
    beta_stack = [thin_films[X].Beta for X in pattern]
    Cv_stack = [thin_films[X].CV for X in pattern]
    kd_stack = [thin_films[X].ThermalDiffusivity for X in pattern]

    # Substrate index
    n_stack.append(substrate.Index)

    # Optical thicknesses: n * d / lambda
    L_opt_stack = [Lj * nj / lam_ref for Lj, nj in zip(L_stack, n_stack[1:])]

    stack = {
        "lam_ref": lam_ref,
        "ns": np.array(n_stack),
        "Ls": np.array(L_stack),
        "Ls_opt": np.array(L_opt_stack),
        "alphas": np.array(a_stack),
        "ctes": np.array(cte_stack),
        "betas": np.array(beta_stack),
        "Ys": np.array(Y_stack),
        "sigmas": np.array(sigma_stack),
        "phis": np.array(phi_stack),
        "Cvs": np.array(Cv_stack),
        "thermaldiffs": np.array(kd_stack),
        "sub": substrate,
        "sup": superstrate,
        "thin_films": thin_films,
        "pattern": pattern,
        "optimized": False,
    }
    return stack
