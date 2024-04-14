# Design and global optimization of dielectric coatings

## Description
Thin film coating design, optimization, and analysis tools based on [G. Venugopalan et. al.](https://doi.org/10.1364/OE.513807)

## Install
From within your favorite python environment (e.g. conda) run:

```bash 
python -m pip install OptimalBragg
```

## Examples

### Quarter-wave high-reflectivity (HR) coating
``` python
import numpy as np
import matplotlib.pyplot as plt
from physunits import um, nm, ppm

from OptimalBragg.materials import *
from OptimalBragg import qw_stack, Material
from OptimalBragg.layers import *
from OptimalBragg.plot import plot_layers, plot_spectral

lam_ref = 1064 * nm
silica = Material(SiO2)
tantala = Material(Ta2O5)
Nlayers = 11

# Makes arbitrary quarter-wave stack
stack = qw_stack(
    lam_ref,
    substrate=silica,
    superstrate=Material(air),
    thin_films={"A": silica, "B": tantala},
    pattern="BA" * Nlayers,
)

# Results
T_ref = trans(lam_ref, stack)
print(Rf"T = {T_ref/ppm:.1f} ppm at {lam_ref/um:.2f} um.")

# Show layer structure and spectral refl/trans
plot_layers(stack)

rel_lambdas = np.linspace(0.75 * lam_ref, 1.25 * lam_ref, 2**10)

# Let's pretend I have custom dispersion data for these thin films
lam_disp = np.array([0.532, 0.633, 0.780, 0.852, 1.064, 1.083, 1.550]) * um
nSiO2 = np.interp(
    rel_lambdas,
    lam_disp,
    [1.4607, 1.4570, 1.4537, 1.4525, 1.4496, 1.4494, 1.444],
)
nTa2O5 = np.interp(
    rel_lambdas,
    lam_disp,
    [2.24, 2.1979, 2.1628, 2.1515, 2.1297, 2.1282, 2.1046],
)

plot_spectral(rel_lambdas, stack)
plot_spectral(rel_lambdas, stack, dispersion={"A": nSiO2, "B": nTa2O5})
plt.show()
```

### Optimize an existing anti-reflective (AR) coating
```python
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