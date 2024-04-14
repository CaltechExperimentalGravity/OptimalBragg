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
