% Code for optimizing dielectric coating design using MATLAB particle swarm
% Developed for optimizing 40m PR3/SR3 coatings, generalized to this version
% https://git.ligo.org/40m/Coatings
% Adapted from RXA code by GV, September 2017
clc
% clear
close all
addpath('../');

%Start a parpool with 40 workers
delete(gcp('nocreate'));
parpool(40);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    USER CONFIG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Some flags
plotFlag = 1;
saveFlag = 1;

%%Load dispersion data...
load dispersion_revised.mat;
na = 1.0000;           %Vacuum
NUMTOOLS.n1_IR = interp1(SiO2(:,1),SiO2(:,2),1064,'pchip');
NUMTOOLS.n2_IR = interp1(Ta2O5(:,1),Ta2O5(:,2),1064,'pchip');
NUMTOOLS.nb_IR = 1.449641;   %From Ramin

NUMTOOLS.n1_green = interp1(SiO2(:,1),SiO2(:,2),532,'pchip');
NUMTOOLS.n2_green = interp1(Ta2O5(:,1),Ta2O5(:,2),532,'pchip');
NUMTOOLS.nb_green = 1.4607; %Estimate

NUMTOOLS.n1_ratio = NUMTOOLS.n1_green / NUMTOOLS.n1_IR;
NUMTOOLS.n2_ratio = NUMTOOLS.n2_green / NUMTOOLS.n2_IR;

% Initial guess vector of layer thicknesses in units of lambda
no_of_stacks = 20;  
x0 = [];
for kk = 1:no_of_stacks
    x0 = [x0 2/8 2/8];
end
 
clear costOut                                                                 
clear ifo    

NUMTOOLS.lambda = 1064e-9;            % Wavelength at which R and T are to be evaluated
NUMTOOLS.opt_name       = 'PR3';
NUMTOOLS.coatingType    = 'HR';
NUMTOOLS.T_1            = 50e-6;      % Target transmission at 1064 nm --- 50ppm for HR, 99% for AR
NUMTOOLS.T_2s           = 0.995;      % Target transmission at 532 nm, s-pol --- 99.5% for HR, 99.8% for AR
NUMTOOLS.T_2p           = 0.995;      % Target transmission at 532 nm, p-pol ---- same as above
NUMTOOLS.aoi            = 41.1;       % Angle of incidence (degrees)  --- 41.1deg for HR, 24.8deg for AR
NUMTOOLS.aoi_green      = 41.1;       % Angle of incidence (degrees)  ---24.746 deg for green
NUMTOOLS.wBeam        = 5e-3;          % Beam size on optic being optimized
NUMTOOLS.f_optimize   = 100;           % Frequency at which to evaluate noise being optimized
NUMTOOLS.noise_weight = 1e42;          % to equate Brownian and TO noise
NUMTOOLS.include_sens   = 1;           %Include derivatives in sensitivity function (1 or 0)
NUMTOOLS.include_thermal   = 1;        %Include TN term in sensitivity function (1 for HR or 0 for AR)


ifo = SilicaTantala300;                 % Load a standard gwinc style parameter file with various material properties defined             
ifo.Laser.Wavelength = NUMTOOLS.lambda;
% load this lookup table for speedup
load besselzeros.mat
ifo.Constants.BesselZeros = besselzeros;
ifo.Optics.ETM.BeamRadius = NUMTOOLS.wBeam; %using the gwinc infrastructure for calculating TO noise for ease...
ifo.Optics.ITM.BeamRadius = NUMTOOLS.wBeam;
NUMTOOLS.ifo = ifo;

%set up the weighting to convert MO cost to scalar cost
% Indices of the yy term
% 1 = Transmission at 1064, p-pol
% 2 = Transmission at 532, p-pol
% 3 = Transmission at 532, s-pol
% 4 = HR surface field
% 5 = Brownian Noise
% 6 = sensitivity to change in layer thickness @1064nm p- +/-1%...
% 7 = sensitivity to change in layer thickness @532nm p-pol +/-1%...
% 8 = sensitivity to change in layer thickness @532nm s-pol +/-1%...
% 9 = sensitivity to change in n1 @1064nm p- +/-1%...
% 10 = sensitivity to change in n1 @532nm p-pol +/-1%...
% 11 = sensitivity to change in n1 @532nm s-pol +/-1%...
% 12 = sensitivity to change in n2 @1064nm p- +/-1%...
% 13 = sensitivity to change in n2 @532nm p-pol +/-1%...
% 14 = sensitivity to change in n2 @532nm s-pol +/-1%...
% 15 = sensitivity to change in aoi @1064nm p- +/-1%...
% 16 = sensitivity to change in aoi @532nm p-pol +/-1%...
% 17 = sensitivity to change in aoi @532nm s-pol +/-1%...

weights = [5555. 1111. 1111. 100. 0.0001 1. 10. 10. 1. 10. 10. 1. 10. 10. 1. 10. 10.];
NUMTOOLS.wt = weights;

% setting the bounds for the variables to be searched over
LB = 0.050 * ones(size(x0));     % lower bounds on the layer thicknesses, lambda/50
UB = 0.51 * ones(size(x0));     % upper bound, lambda/2          
nvars = length(UB);

if strcmp(NUMTOOLS.coatingType,'HR')
    LB(1) = 0.48;  %Constrain topost layer to be nearly halfwave thick...
    UB(1) = 0.52;
	NUMTOOLS.firstLayer = 'SiO2';     %Depending whether coating is on the HR side or AR side
else
	NUMTOOLS.firstLayer = 'Ta2O5';
end

% set particle swarm optimization options
hybridopts = optimoptions('fmincon',...
                  'Display','iter',...
                  'MaxIter', 911,...
                  'TolFun', 1e-4,...
                  'MaxFunEvals', 11111);

options = optimoptions('particleswarm',...
               'SwarmSize', nvars*155,...     % 
               'UseParallel', 1,...
               'MaxIter', 911,...
               'SelfAdjustment',   1.49,...
               'SocialAdjustment', 1.49,...
               'TolFun', 1e-4,...
               'Display', 'iter',...
               'HybridFcn',{@fmincon, hybridopts});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   END USER CONFIG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Run particle swarm optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = x0(:);                       
%% Do the optimization
% RUNS the Particle Swarm ========
tic
[xout, fval, exitflag] =...
    particleswarm(@(x) getMirrorCost([x], NUMTOOLS, 0),...
                                    nvars, LB, UB, options);
toc

% Check for thin layers
if find(xout < 0.05)
    disp('Bad Layer Thickness: invalid results')
end         
%% SAVE layer structure, IFOmodel, and noise
costOut = getMirrorCost(xout, NUMTOOLS, 1);
tnowstr = datestr(now, 'yymmdd_HHMM');
savename = ['Data/' NUMTOOLS.opt_name '_' NUMTOOLS.coatingType '_' num2str(no_of_stacks) '_layers_' tnowstr];
save(savename, 'costOut', 'ifo');

%%Make a spectral reflectivity plot
makeSpecREFLplot(savename,NUMTOOLS.firstLayer);

%% printing plot to file
title(' ')
orient landscape
set(gcf,'PaperPositionMode','auto')

fname = ['Figures/' NUMTOOLS.opt_name '_R_' tnowstr];

print('-depsc','-r600', fname)
[a,b] = system(['../makePDF.sh ' fname '.eps']);
if a ~= 0
    disp('PDF Generation Error')
else
    [a,b] = system(['rm ' fname '.eps']);
end

%% plot layer structure and thermo-optic noise
% plotTO120
% exit;
