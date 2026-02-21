import time
import h5py
import numpy as np
import matplotlib.pyplot as plt
from physunits import cm, inch, um, nm, ppm, Hz

from OptimalBragg.materials import *
from OptimalBragg import qw_stack, h5write, Material
from OptimalBragg.layers import *
from OptimalBragg.noise import coating_noise
from OptimalBragg.plot import plot_layers, plot_spectral, plot_noises
from OptimalBragg.optimizer import diff_evo

lam_psl = 1064 * nm
lam_als = 0.5*lam_psl
optic = 'ITM'


if __name__ == '__main__':
    # Initialize QW stack but override with user defined pre-designed stack
    stack = qw_stack(
        lam_ref=lam_psl,
        substrate=Material(SiO2),
        superstrate=Material(air),
        thin_films={"L": Material(SiO2), "H": Material(Ta2O5)},
        pattern="LH" * 15,
        hwcap="H",
    )
    
    T_ref = trans(lam_psl, stack)
    stack["T_ref"] = T_ref
    print(Rf"T < {(T_ref)*100:.8f} % at {lam_psl/um:.2f} um.")

    # Optimization over multiple wavelength AR and absorption
    T_p = trans(lam_psl, stack)
    T_a = trans(lam_als, stack)

    ffreqs = np.logspace(0, 4, 2**10)
    mirror_mass = 0.25
    mirror_diam = 3*inch
    mirror_thic = 1*inch
    if optic == 'ETM':
        w_TM_psl = 1.4*cm 
        target_T_psl = 13.7*ppm
    elif optic == 'ITM':
        w_TM_psl = 1.0*cm
        target_T_psl = 13840*ppm
    w_TM_aux = w_TM_psl / np.sqrt(2)
    power_psl = 1e3
    power_aux = 0.01

    TO_pars = {"freq":100, "w_beam":w_TM_psl}
    Br_pars = {"freq":100, "w_beam":w_TM_psl, "power":power_psl, "m_mirror":mirror_mass}
    
    # Calc noise
    coat_Sbr, coat_Ste, coat_Str, coat_Sto = coating_noise(ffreqs, stack, w_beam=w_TM_psl, power=power_psl, r_mirror=0.5*mirror_diam, d_mirror=mirror_thic, m_mirror=mirror_mass)
    noise_labels = {"Sbr":coat_Sbr, "Sto":coat_Sto}
    plot_noises(ffreqs, noise_labels, plot_total=True)
    plt.show()

    # Reference (initial) stack
    stack["init"] = {"ns": stack["ns"], "Ls": stack["Ls"], "T_ref": T_ref}
    multi_target = {
        "T": {
            "target": {
                lam_psl: target_T_psl,
                lam_als: 100 * ppm,
            },
            "weight": {lam_psl: 1, lam_als: 1},
        },
        # "Sbr": {"target":1e-45,
        #         "weight":1e-1},
        # "Sto": {"target":1e-45,
        #         "weight":1e-1},
        # "abs": {"target": 50 * ppm, "weight": 1e-2},
    }
    
    optimization_result = diff_evo(stack, multi_target, to_pars=TO_pars,
    br_pars=Br_pars)
    stack["optimized"] = True

    # Update thicknesses and other optimized attributes
    stack["Ls"] = optimization_result["Ls"]
    T_ref = trans(lam_psl, stack)
    T_als = trans(lam_als, stack)

    _, Enorm = field_zmag(
        stack["ns"], stack["Ls"], n_pts=2**8, lam=stack["lam_ref"]
    )
    intAbs = calc_abs(Enorm, stack["Ls"], stack["alphas"])
    stack["Absorption"] = intAbs
    stack["T_ref"] = T_ref

    # Results
    plot_layers(stack)
    plt.show()

    wavelengths = np.linspace(0.95 * lam_als, 1.15 * lam_psl, 2**12)
    plot_spectral(wavelengths, stack, markers={"T": [lam_als, lam_psl]})
    plt.show()

    # Calc noise
    coat_Sbr, coat_Ste, coat_Str, coat_Sto = coating_noise(ffreqs, stack, w_beam=w_TM_psl, power=power_psl, r_mirror=0.5*mirror_diam, d_mirror=mirror_thic, m_mirror=mirror_mass)
    noise_labels = {"Sbr":coat_Sbr, "Sto":coat_Sto}
    plot_noises(ffreqs, noise_labels, plot_total=True)
    plt.show()
    
    # Save to hdf5
    time_tag = time.strftime("%Y%m%d-%H%M%S")
    h5write(
        Rf"./{optic}_Tpsl_{(T_ref)/ppm:.0f}ppm_Taux_{(T_als)/ppm:.0f}_ppm_Apsl_{intAbs/ppm:.0f}ppm_{time_tag}.h5",
        stack,
    )