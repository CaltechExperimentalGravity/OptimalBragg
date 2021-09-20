#! /usr/bin/python

import sys, glob, os
import h5py
from generic_local.coatingUtils import *

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

try:
    plt.style.use("pacostyle")
else:
    plt.style.use("default")

from matplotlib.ticker import FormatStrFormatter

if len(sys.argv) == 1:
    fname = max(glob.iglob("Data/ETM/*Layers*.hdf5"), key=os.path.getctime)
    # fname = fname[5:] # rm 'data' from the name
else:
    # For example fname = 'ETM_Layers_190519_1459.hdf5'
    fname = str(sys.argv[1])


def h5read(targets):
    data = {}
    with h5py.File(fname, "r+") as f:
        for target in targets:
            try:
                data[target] = np.array(f["diffevo_output/" + target])
            except:
                data[target] = 0.0
    return data


def read_dispersion_data(filepath):
    return np.genfromtxt(filepath, delimiter=",", skip_header=1)


def interpolate_dispersions():
    """Interpolate dispersion curves using spline

    Returns:
        ndarray: Complex dispersion
    """

    # Read fused silica dispersion
    fs_short_wave = read_dispersion_data("generic_local/SiO2_n_k_short.csv")
    fs_long_wave = read_dispersion_data("generic_local/SiO2_n_k_long.csv")
    fused_silica_n_k = np.concatenate((fs_short_wave, fs_long_wave))

    # Read tantala dispersion
    tantala_n_k = read_dispersion_data("generic_local/Ta2O5_n_k.csv")

    # Read silicon dispersion
    si_short_wave = read_dispersion_data("generic_local/Silicon_n_k_short.csv")
    si_long_wave = read_dispersion_data("generic_local/Silicon_n_k_long.csv")
    silicon_n_k = np.concatenate((si_short_wave, si_long_wave))

    # Use cubic spline through scipy interpolate
    fs_n_interpolator = interp1d(
        fused_silica_n_k[:, 0],
        fused_silica_n_k[:, 1],
        kind="cubic",
        fill_value=(fused_silica_n_k[0, 1], fused_silica_n_k[-1, 1]),
        bounds_error=False,
    )
    fs_k_interpolator = interp1d(
        fused_silica_n_k[:, 0],
        fused_silica_n_k[:, 2],
        kind="cubic",
        fill_value=(fused_silica_n_k[0, 2], fused_silica_n_k[-1, 2]),
        bounds_error=False,
    )

    ta_n_interpolator = interp1d(
        tantala_n_k[:, 0],
        tantala_n_k[:, 1],
        kind="cubic",
        fill_value=(tantala_n_k[0, 1], tantala_n_k[-1, 1]),
        bounds_error=False,
    )
    ta_k_interpolator = interp1d(
        tantala_n_k[:, 0],
        tantala_n_k[:, 2],
        kind="cubic",
        fill_value=(tantala_n_k[0, 2], tantala_n_k[-1, 2]),
        bounds_error=False,
    )

    si_n_interpolator = interp1d(
        silicon_n_k[:, 0],
        silicon_n_k[:, 1],
        kind="cubic",
        fill_value=(silicon_n_k[0, 1], silicon_n_k[-1, 1]),
        bounds_error=False,
    )
    si_k_interpolator = interp1d(
        silicon_n_k[:, 0],
        silicon_n_k[:, 2],
        kind="cubic",
        fill_value=(silicon_n_k[0, 2], silicon_n_k[-1, 2]),
        bounds_error=False,
    )

    return (
        fs_n_interpolator,
        fs_k_interpolator,
        ta_n_interpolator,
        ta_k_interpolator,
        si_n_interpolator,
        si_k_interpolator,
    )


def emissivity(reflectivity, transmissivity):
    """Compute emissivity assuming it's equal to absorptivity

    Args:
        reflectivity (float): wavelength dependent reflectivity
        transmissivity (float): wavelength dependent transmissivity

    Returns:
        emissivity (float): emissivity (between 0 and 1)
    """
    return 1 - reflectivity - transmissivity


if __name__ == "__main__":
    data = h5read(targets=["n", "L"])
    n_data = data["n"]
    L_data = data["L"]
    print(fname)
    # Read interpolated dispersions of thin films
    fs_n, fs_k, ta_n, ta_k, si_n, si_k = interpolate_dispersions()
    wavelengths = np.linspace(1 * um, 50 * um, 2 ** 12) / um

    # Evaluate the reflectivity of the stack at
    # each wavelength without extinction coefficients
    emissivity_0 = np.zeros_like(wavelengths)
    emissivity_1 = np.zeros_like(wavelengths)

    for j, wavelength in enumerate(wavelengths):
        # Evaluate dispersion
        n_low = fs_n(wavelength) - 1j * fs_k(wavelength)
        n_high = ta_n(wavelength) - 1j * ta_k(wavelength)
        n_bulk = si_n(wavelength) - 1j * si_k(wavelength)

        # Construct stack index array using n_data array
        n_mask = n_data.copy().astype(complex)
        n_mask[n_mask < 2] = n_low
        n_mask[n_mask > 2] = n_high
        n_mask[0] = 1.0
        n_mask[-1] = n_bulk

        L_corr = L_data * n_mask[1:-1] / n_data[1:-1]
        dz = np.real(L_corr / n_mask[1:-1])
        # print(dz * n_mask[1:-1].imag)
        attenuation = np.exp(
            2 * dz * n_mask[1:-1].imag * 2 * np.pi / wavelength / 2.1282
        )
        total_att = np.prod(attenuation)

        # Compute reflectivity, and impedance in the presence of nonzero ext. coeff.
        rj, zj = multidiel1(n_mask, L_corr, wavelength / 2.1282)
        Tj = 1 - np.abs(rj) ** 2

        # print(Tj, total_att)

        # From the fraction that doesn't get reflected (i.e. transmitted)
        # compute the absorption (i.e. emission) by the total attenuation
        emissivity_1[j] = total_att * Tj

    # plt.figure()
    # plt.plot(wavelengths, emissivity_1)
    # plt.ylim(-0.25, 1.15)
    # plt.show()

    plt.figure()
    plt.plot(wavelengths, fs_k(wavelengths), label="SiO2")
    plt.plot(wavelengths, ta_k(wavelengths), label="Ta2O5")
    plt.show()
