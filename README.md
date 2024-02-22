# Design and global optimization of dielectric coatings

## Description

## Install
From within your favorite python environment (e.g. conda) run:

```bash 
python -m pip install coatings
```

## Examples

### Quarter-wave high-reflectivity (HR) coating
``` python
import numpy as np
import matplotlib.pyplot as plt
from physunits import um, nm, ppm

from materials import *
from utils import qw_stack, Material
from utils.layers import stack_R
from utils.plot import plot_layers, plot_spectral

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
T_ref = 1 - stack_R(lam_ref, stack)
print(Rf"T = {T_ref/ppm:.1f} ppm at {lam_ref/um:.2f} um.")

# Show layer structure and spectral refl/trans
plot_layers(stack)

rel_lambdas = np.linspace(0.85 * lam_ref, 1.15 * lam_ref, 2**8)
plot_spectral(rel_lambdas, stack)
plt.show()
```

### Optimize an existing anti-reflective (AR) coating
```python
import time
import h5py
import numpy as np
import matplotlib.pyplot as plt
from physunits import um, nm, ppm, Hz

from materials import *
from utils import qw_stack, h5write, Material
from utils.layers import stack_R, field_zmag, calc_abs
from utils.plot import plot_layers, plot_spectral
from utils.optimizer import diff_evo

lam_ref = 1550 * nm

# Initialize QW stack but override with user defined pre-designed stack
stack = qw_stack(
    lam_ref=lam_ref,
    substrate=Material(SiO2),
    superstrate=Material(air),
    thin_films={"L": Material(SiO2), "H": Material(Ta2O5)},
    pattern="LH" * 2,
    hwcap="H",
)
stack["ns"] = np.array([1.0, 2.1, 1.45, 2.1, 1.45, 2.1, 1.45])
stack["Ls"] = np.array([0.8548, 268.4, 204.2, 90.18, 61.42]) * nm
T_ref = 1 - stack_R(lam_ref, stack)
stack["T_ref"] = T_ref
print(Rf"R < {(1 - T_ref)*100:.4f} % at {lam_ref/um:.2f} um.")

# Optimization over multiple wavelength AR and absorption
lam_m = 1545 * nm
lam_p = 1564 * nm
T_p = 1 - stack_R(lam_p, stack)
T_m = 1 - stack_R(lam_m, stack)

# Reference (initial) stack
stack["init"] = {"ns": stack["ns"], "Ls": stack["Ls"], "T_ref": T_ref}
multi_target = {
    "T": {
        "target": {
            lam_ref: 1 - 200 * ppm,
            lam_p: 1 - 500 * ppm,
            lam_m: 1 - 500 * ppm,
        },
        "weight": {lam_ref: 25, lam_p: 15, lam_m: 15},
    },
    "Absorption": {"target": 25 * ppm, "weight": 1e-5},
}
optimization_result = diff_evo(stack, multi_target)
stack["optimized"] = True

# Update thicknesses and other optimized attributes
stack["Ls"] = optimization_result["Ls"]
T_ref = 1 - stack_R(lam_ref, stack)
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
plot_spectral(wavelengths, stack, markers=[lam_p, lam_m])
plt.show()

# Save to hdf5
time_tag = time.strftime("%Y%m%d-%H%M%S")
h5write(
    Rf"./AR1550_R_{(1-T_ref)/ppm:.0f}_A_{intAbs/ppm:.0f}_ppm_{time_tag}.h5",
    stack,
)
