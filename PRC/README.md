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
