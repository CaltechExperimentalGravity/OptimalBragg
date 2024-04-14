import time
import h5py
import numpy as np
import matplotlib.pyplot as plt
from physunits import um, nm, ppm, Hz

from OptimalBragg.materials import *
from OptimalBragg import qw_stack, h5write, Material
from OptimalBragg.layers import *
from OptimalBragg.plot import plot_layers, plot_spectral
from OptimalBragg.optimizer import diff_evo

lam_ref = 1550 * nm

# Initialize QW stack but override with user defined pre-designed stack
stack = qw_stack(
    lam_ref=lam_ref,
    substrate=Material(SiO2),
    superstrate=Material(air),
    thin_films={"L": Material(SiO2), "H": Material(Ta2O5)},
    pattern="LH" * 4,
    hwcap="H",
)
# stack["ns"] = np.array([1.0, 2.1, 1.45, 2.1, 1.45, 2.1, 1.45])
# stack["Ls"] = np.array([0.8548, 268.4, 204.2, 90.18, 61.42]) * nm
T_ref = trans(lam_ref, stack)
stack["T_ref"] = T_ref
print(Rf"R < {(1 - T_ref)*100:.8f} % at {lam_ref/um:.2f} um.")

# Optimization over multiple wavelength AR and absorption
lam_m = 1545 * nm
lam_p = 1564 * nm
T_p = trans(lam_p, stack)
T_m = trans(lam_m, stack)

# Reference (initial) stack
stack["init"] = {"ns": stack["ns"], "Ls": stack["Ls"], "T_ref": T_ref}
multi_target = {
    "R": {
        "target": {
            lam_ref: 10 * ppm,
            lam_p: 10 * ppm,
            lam_m: 10 * ppm,
        },
        "weight": {lam_ref: 1, lam_p: 1, lam_m: 1},
    },
    "abs": {"target": 25 * ppm, "weight": 1e-2},
}
optimization_result = diff_evo(stack, multi_target)
stack["optimized"] = True

# Update thicknesses and other optimized attributes
stack["Ls"] = optimization_result["Ls"]
T_ref = trans(lam_ref, stack)
_, Enorm = field_zmag(
    stack["ns"], stack["Ls"], n_pts=2**8, lam=stack["lam_ref"]
)
intAbs = calc_abs(Enorm, stack["Ls"], stack["alphas"])
stack["Absorption"] = intAbs
stack["T_ref"] = T_ref

# Results
plot_layers(stack)
plt.show()

wavelengths = np.linspace(0.95 * lam_m, 1.05 * lam_p, 2**12)
plot_spectral(wavelengths, stack, markers={"R": [lam_p, lam_m, lam_ref]})
plt.show()

# Save to hdf5
time_tag = time.strftime("%Y%m%d-%H%M%S")
h5write(
    Rf"./AR1550_R_{(1-T_ref)/ppm:.0f}_A_{intAbs/ppm:.0f}_ppm_{time_tag}.h5",
    stack,
)
