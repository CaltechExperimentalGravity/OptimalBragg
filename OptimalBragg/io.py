"""HDF5 and YAML I/O utilities for OptimalBragg.

Provides reading/writing of coating designs, optimization results,
and material configuration files.
"""

import h5py
import numpy as np
import yaml


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
        h5_dict (dict): Dictionary of depth <= 3.
            Material objects are serialized as HDF5 groups
            with their attributes stored as HDF5 attrs.
    """
    from OptimalBragg import Material

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
        dict: Output parameters
    """
    with open(fname, "r") as f:
        params = yaml.safe_load(f)
    return params


def load_materials_yaml(fname):
    """Read a project materials.yml file and resolve material names
    against the OptimalBragg materials library.

    The materials.yml format::

        substrate:
          material: SiO2          # name in OptimalBragg.materials
          overrides:
            Temp: 295
            MassRadius: 0.17
        thin_films:
          L:
            material: SiO2
          H:
            material: TiTa2O5
        laser:
          wavelength: 1064e-9
          power: 125
        optics:
          ETM:
            beam_radius: 0.062
          ITM:
            beam_radius: 0.055

    Args:
        fname (str): Path to materials.yml file

    Returns:
        dict: Contains keys ``substrate`` (Material), ``thin_films``
            (dict of ``{"L": Material, "H": Material}``), ``laser``
            (dict), and ``optics`` (dict).
    """
    import OptimalBragg.materials as matlib
    from OptimalBragg import Material

    raw = yamlread(fname)

    # Resolve substrate
    sub_name = raw["substrate"]["material"]
    sub_dict = getattr(matlib, sub_name)
    # Deep copy and apply overrides
    sub_props = dict(sub_dict["Properties"])
    sub_props.update(raw["substrate"].get("overrides", {}))
    substrate = Material({"Properties": sub_props})

    # Resolve thin films
    thin_films = {}
    for key, spec in raw["thin_films"].items():
        if isinstance(spec, str):
            mat_dict = getattr(matlib, spec)
            thin_films[key] = Material(mat_dict)
        else:
            mat_name = spec["material"]
            mat_dict = getattr(matlib, mat_name)
            mat_props = dict(mat_dict["Properties"])
            mat_props.update({k: v for k, v in spec.items() if k != "material"})
            thin_films[key] = Material({"Properties": mat_props})

    # Ensure numeric values in laser config are floats (YAML safe_load
    # treats values like '1064e-9' as strings when there's no decimal point)
    laser = raw.get("laser", {})
    for k, v in laser.items():
        if isinstance(v, str):
            try:
                laser[k] = float(v)
            except ValueError:
                pass

    result = {
        "substrate": substrate,
        "thin_films": thin_films,
        "laser": laser,
        "optics": raw.get("optics", {}),
    }
    return result
