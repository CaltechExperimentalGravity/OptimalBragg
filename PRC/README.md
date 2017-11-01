# PR3/SR3 coating Design

#### notes

Set T_1064 = 10 ppm and T_532 = 99.9%, seems to work now.


### GV Feb 22 2017 ###
Made some changes to the code to optimize R for 1064nm p-pol and T for 532nm p and s-pol on the HR side, and take into account non-normal incidence
Also calculated AR side coating, with only R @ 532nm for p and s-pol used to compute the error function
Current best params:
HR Side
T_1064_p = 50.7 ppm
T_532_p = 0.99
T_532_s = 0.99

AR side
R_532_p = 51.5ppm
R_532_s = 50ppm
nothing special for 1064nm...

### GV April 16 2017 ###
- Cleaned up some scripts
- separate functions for calculating derivatives, making plot 
- now using dispersion data from Ramin for calculations as well as plotting
- changed constraints on layer thickness, first layer now required to be nearly half-wave thick @1064nm, all layers thicker than ~50nm

#GV July 04 2017
- Updating to version of code used to generate E1700016-v9 [https://dcc.ligo.org/E1700016-v9]
- Major changes are 
	a. Fixed bug in code that wasn't taking into account Dispersion properly
	b. Tweaked sensitivity function slightly
- MATLAB functions are:
	1. runPR3.m            --- This allows multiple calls to the function doPR3.m, with different number of layer pairs.
	2. doPR3.m             --- Main function in which T / R targets are to be defined. Most user config has to be done in here.
	3. getMirrorCost.m     --- Calculates the cost function. User may need to tweak weights to various terms in here for better performance, or add terms to cost function.
	4. doSens.m            --- Calculates numerical derivatives w.r.t. various params, for reducing sensitivity of designed coating to said params.
	5. theta2.m            --- Function that calculates angle of incidence in second media for Snel's Law.
	6. makeSpecREFLplot.m  --- Function that makes the spectran R/T plots for both polarizations for a given coating design, takes dispersion into account.
	7. op2phys.m           --- Function that converts optical thicknesses from mutidiel1.m to physical thicknesses. Required for taking into account dispersion properly.
	8. plotLayers.m        --- Makes a bar-histogram style plot of the coating design.
	9. doMC.m              --- Calculates R/T for a given coating design with perturbations to various parameters. Ensemble size, errors etc can be tweaked inside this function.
				   I have set this up right now to export the results for plotting with Python.
	10. plotTO120.m        --- makes plot of thermal noise for a given coating design (not tweaked since porting over from Rana's code)
	11. SilicaTantala300.m --- GWINC style IFO parameter file. This probably can be gotten rid of for speed of code execution.
