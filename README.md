# Design of Mirror Coatings by global optimization of a cost function

For the mirror coatings, we have many constraints to satisfy:

1. Transmissivity of XX ppm at a wavelength of 1064 nm.
1. Similar at other wavelengths
1. Minimize Brownian thermal noise
1. Minimize Thermo-Optic noise
1. Minimize sensitivity of transmissivity to coating deposition errors
1. Minimize E-field at HR surface
1. Layers cannot be less than 1 nm thick or more than 1 micron.

# Description of functions
The `generic` directory contains the most up-to-date versions of the various scripts.
 * `multidiel1.m` --- used to calculate coating reflectivity.
 * `getMirrorCost.m` --- Generic cost function template which will be used by global optimizer, e.g. PSO.
 * `runSwarm.m` --- Generic calling function that sets up the optimization problem and calls the PSO.
 * `doSens.m` --- Perturbs a given parameter by +/- 1% and computes sensitivity of coating reflectivity to this perturbation.
 * `thermalNoiseFuncs` --- Contains a collection of GWINC functions that compute thermal noise psds.
 * `calcEField.m` --- Computes the E-field squared (normalized to surface E-field units) as a function of penetration depth.
 * `op2phys.m` --- Convert optical thickness to physical thickness (in units of lambda0)
 * `theta2.m` --- Snel's law angle calculator
 * `makeSpecREFLplot.m` --- as the name suggests
 * `plotLayers.m` --- as the name suggests
 * `pythonAddOns` --- a set of functions for making matplotlib plots from the output of PSO optimization. Also, for MC Hammering. See below.
All other functions in `generic` are adaptations of the above list to specific design cases.


# Monte-Carlo analysis of sensitivity of coating design to assumed model parameters
Once a coating design has been generated, we'd like to see how sensitive it is to 
 * Manufacturing tolerances
 * Assumed values of model parameters, e.g. refractive indices, dispersion, angle of incidence etc
Tools to do this sort of analysis are available in `generic/pythonAddOns`. A description of the relevant files:
1. `coatingUtils.py` has some functions for calculating reflectivity etc. for a given coating design.
2. `doMC.py` takes in the output file from the MATLAB PSO optimization, and perturbs the following:
	* Angle of incidence
	* Refractive index of high and low index layers (systematically, by the same fractional amount for each type of layer)
	* Thickness of layers (systematically, by the same fractional amount)
3. In the above step, the perturbed variables are assumed to be i.i.d. Gaussian random variables, with a width of 0.5%
4. The `emcee` package is used for sampling from the 4D multivariate distribution described above.
5. We are interested in the effect of the perturbations on:
	* Reflectivity at wavelengths of interest
	* Thermal noise properties (you'll need the python version of `gwinc` to evaluate these)
	* Surface E-field 
	* Absorption (to be added)
6. `cornerPlt.py` takes in the output file from `doMC.py` and generates a visualization of the MC simulation.
