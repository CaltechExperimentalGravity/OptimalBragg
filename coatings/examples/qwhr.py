import numpy as np
import matplotlib.pyplot as plt
from physunits import um, nm, ppm

from coatings.materials import *
from coatings import qw_stack, Material
from coatings.layers import stack_R
from coatings.plot import plot_layers, plot_spectral

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
