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


def _h5_write_value(group, key, value):
    """Write a single value into an HDF5 group.

    Scalars (int, float, str) are written as scalar datasets so they
    can be read back with ``group[key]``.  Material objects are stored
    as sub-groups with their properties as attrs.
    """
    from OptimalBragg import Material

    if isinstance(value, Material):
        g = group.create_group(value.Name)
        g.attrs.create("key", key)
        for mk, mv in value.__dict__.items():
            g.attrs.create(mk, mv)
    elif isinstance(value, dict):
        g = group.create_group(key)
        for kk, vv in value.items():
            _h5_write_value(g, kk, vv)
    else:
        group.create_dataset(key, data=value)


def h5write(fname, h5_dict):
    """Helper function to write hdf5 files with
    coating designs, optimization, results, etc.

    Args:
        fname (str): Path to h5file
        h5_dict (dict): Dictionary of arbitrary depth.
            Material objects are serialized as HDF5 groups
            with their attributes stored as HDF5 attrs.
    """
    with h5py.File(fname, "w") as f:
        for k, v in h5_dict.items():
            _h5_write_value(f, k, v)


def _coerce_numerics(obj):
    """Recursively convert string values that look like numbers to floats.

    YAML safe_load treats values like ``5e-6`` or ``1064e-9`` (no decimal
    point) as strings.  This function walks a nested dict/list and converts
    any such strings to ``float``.
    """
    if isinstance(obj, dict):
        return {k: _coerce_numerics(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_coerce_numerics(v) for v in obj]
    if isinstance(obj, str):
        try:
            return float(obj)
        except ValueError:
            return obj
    return obj


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
    return _coerce_numerics(params)


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

    result = {
        "substrate": substrate,
        "thin_films": thin_films,
        "laser": raw.get("laser", {}),
        "optics": raw.get("optics", {}),
    }
    return result
