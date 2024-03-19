import h5py
import yaml
import numpy as np

from .layers import *
from .optimizer import *
from .plot import *


class Material:
    def __init__(self, props):
        self.__dict__.update(props["Properties"])


def h5read(fname, group, targets):
    """Helper function to read hdf5 files with
    coating designs, optimization, results, etc.

    Args:
        fname (str): Path to h5file
        group (str): Group
        targets (list): Datasets to read

    Returns:
        dict: Output datasets
    """
    data = {}
    with h5py.File(fname, "r") as f:
        for target in targets:
            try:
                data[target] = np.array(f[group][target])
            except:
                data[target] = 0.0
    return data


def h5write(fname, h5_dict):
    """Helper function to write hdf5 files with
    coating designs, optimization, results, etc.

    Args:
        fname (str): Path to h5file
        targets (list): Dictionary of depth <= 3
    """
    with h5py.File(fname, "w") as f:
        for k, v in h5_dict.items():
            if isinstance(v, dict):
                f.create_group(k)
                for kk, vv in v.items():
                    if isinstance(vv, dict):
                        f[k].create_group(kk)
                        for kkk, vvv in vv.items():
                            if (
                                isinstance(vvv, int)
                                or isinstance(vvv, float)
                                or isinstance(vvv, str)
                            ):
                                f[k][kk].attrs.create(kkk, vvv)
                            elif isinstance(vvv, Material):
                                f[k][kk].create_group(vvv.Name)
                                f[k][kk][vvv.Name].attrs.create("key", kkk)
                                for mk, mv in vvv.__dict__.items():
                                    f[k][kk][vvv.Name].attrs.create(mk, mv)
                            else:
                                f[k][kk].create_dataset(kkk, data=vvv)
                    elif isinstance(vv, Material):
                        f[k].create_group(vv.Name)
                        f[k][vv.Name].attrs.create("key", kk)
                        for mk, mv in vv.__dict__.items():
                            f[k][vv.Name].attrs.create(mk, mv)
                    else:
                        if (
                            isinstance(vv, int)
                            or isinstance(vv, float)
                            or isinstance(vv, str)
                        ):
                            f[k].attrs.create(kk, vv)
                        else:
                            f[k].create_dataset(kk, data=vv)
            elif isinstance(v, Material):
                f.create_group(v.Name)
                f[v.Name].attrs.create("key", k)
                for mk, mv in v.__dict__.items():
                    f[v.Name].attrs.create(mk, mv)
            else:
                if (
                    isinstance(v, int)
                    or isinstance(v, float)
                    or isinstance(v, str)
                ):
                    f.attrs.create(k, v)
                else:
                    f.create_dataset(k, data=v)


def yamlread(fname):
    """Helper function to read yaml file with parameters
    for the coating design, optimization, results, etc.

    Args:
        fname (str): Path to yaml file

    Returns:
        params (dict): Output parameters

    """
    with open(fname, "r") as f:
        params = yaml.safe_load(f)
    return params


def qw_stack(lam_ref, substrate, superstrate, thin_films, pattern, hwcap=""):
    n_stack = [superstrate.Index]
    L_stack = []
    if hwcap:
        n_stack += [thin_films[X].Index for X in hwcap]
        L_stack = [(lam_ref / 2) / nj for nj in n_stack[1:]]
    # Rest of quarter-wave stack
    n_stack += [thin_films[X].Index for X in pattern]
    L_stack += [(lam_ref / 4) / nj for nj in n_stack[len(hwcap) + 1 :]]

    # Complete pattern if needed
    pattern = hwcap + pattern
    a_stack = [thin_films[X].Absorption for X in pattern]
    Y_stack = [thin_films[X].Y for X in pattern]
    sigma_stack = [thin_films[X].Sigma for X in pattern]
    phi_stack = [thin_films[X].Phi for X in pattern]
    cte_stack = [thin_films[X].Alpha for X in pattern]
    beta_stack = [thin_films[X].Beta for X in pattern]
    Cv_stack = [thin_films[X].CV for X in pattern]
    kd_stack = [thin_films[X].ThermalDiffusivity for X in pattern]

    # Substrate
    n_stack.append(substrate.Index)

    stack = {
        "lam_ref": lam_ref,
        "ns": np.array(n_stack),
        "Ls": np.array(L_stack),
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
