% Code for optimizing dielectric coating design using MATLAB particle swarm
% Adapted to optimize aLIGO ETM coating requirements - E0900068-v5
% https://git.ligo.org/40m/Coatings
clc
% clear
close all
addpath('../');
addpath('thermalNoiseFuncs/');

%Start a parpool with 40 workers
delete(gcp('nocreate'));
parpool(40);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    USER CONFIG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

NUMTOOLS.lambda         = 1064e-9;  % Wavelength at which R and T are to be evaluated
NUMTOOLS.opt_name       = 'aLIGO_ETM';
NUMTOOLS.coatingType    = 'HR';
NUMTOOLS.T_1            = 5e-6;     % Target transmission at 1064 nm --- 5ppm for HR
NUMTOOLS.T_2p           = 0.02;     % Target transmission at 532 nm, p-pol ---- same as above
NUMTOOLS.surfaceField   = 0.01;     % Target Surface E Field, [V/m]
NUMTOOLS.aoi            = 0;       % Angle of incidence (degrees)  --- 41.1deg for HR, 24.8deg for AR
NUMTOOLS.aoi_green      = 0;       % Angle of incidence (degrees)  ---24.746 deg for green
NUMTOOLS.wBeam        = 6e-2;       % Beam size on optic being optimized
NUMTOOLS.f_optimize   = 100;        % Frequency at which to evaluate noise being optimized
NUMTOOLS.noise_weight = 1e21;       % to equate Brownian and TO noise
NUMTOOLS.include_sens   = 0;        % Include derivatives in sensitivity function (1 or 0)
NUMTOOLS.include_brownian  = 1;     % Include coating brownian term in sensitivity function
NUMTOOLS.include_TO   = 0;          % Include thermo-optic term in sensitivity function


ifo = SilicaTantala300;                 % Load a standard gwinc style parameter file with various material properties defined             
ifo.Laser.Wavelength = NUMTOOLS.lambda;
ifo.Materials.Coating.Phihighn = loss_highn;
% load this lookup table for speedup
load besselzeros.mat
ifo.Constants.BesselZeros = besselzeros;
ifo.Optics.ETM.BeamRadius = NUMTOOLS.wBeam; %using the gwinc infrastructure for calculating TO noise for ease...
ifo.Optics.ITM.BeamRadius = NUMTOOLS.wBeam;
%set some of the material params in this structure to match that used to optimize coating, i.e. values from Ramin
ifo.Materials.Coating.Indexhighn = NUMTOOLS.n2_IR;
ifo.Materials.Coating.Indexlown  = NUMTOOLS.n1_IR;
ifo.Materials.Substrate.RefractiveIndex = 1.449641; % Corning datasheet
NUMTOOLS.ifo = ifo;

%set up the weighting to convert MO cost to scalar cost
% Indices of the yy term
% 1 = Transmission at 1064, p-pol
% 2 = Transmission at 532, p-pol
% 3 = HR surface field
% 4 = Brownian Noise
% 5 = TO Noise
% 6 = sensitivity to change in layer thickness @1064nm p- +/-1%...
% 7 = sensitivity to change in layer thickness @532nm p-pol +/-1%...
% 8 = sensitivity to change in n1 @1064nm p- +/-1%...
% 9 = sensitivity to change in n1 @532nm p-pol +/-1%...
% 10 = sensitivity to change in n2 @1064nm p- +/-1%...
% 11 = sensitivity to change in n2 @532nm p-pol +/-1%...
% 12 = sensitivity of surface field term to perturbations in n1, n2 and L

weights = [10000 1000 1 10 1e22 5 1 5 1 5 1 0.01];
NUMTOOLS.wt = weights;

% setting the bounds for the variables to be searched over
LB = 0.030 * ones(size(x0));     % lower bounds on the layer thicknesses, lambda/50
UB = 0.51 * ones(size(x0));     % upper bound, lambda/2          
nvars = length(UB);

if strcmp(NUMTOOLS.coatingType,'HR')
%    LB(1) = 0.45;  %Constrain topost layer to be nearly halfwave thick...
%    UB(1) = 0.55;
	NUMTOOLS.firstLayer = 'SiO2';     %Depending whether coating is on the HR side or AR side
else
	NUMTOOLS.firstLayer = 'Ta2O5';
end

% set particle swarm optimization options
hybridopts = optimoptions('fmincon',...
                  'Display','iter',...
                  'MaxIter', 911,...
                  'TolFun', 1e-4,...
                  'MaxFunEvals', 51111);

options = optimoptions('particleswarm',...
               'SwarmSize', nvars*55,...
               'UseParallel', 1,...
               'MaxIter', 1111,...
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
    particleswarm(@(x) getCost_aLIGO_ETM([x], NUMTOOLS, 0),...
                                    nvars, LB, UB, options);
toc

% Check for thin layers
if find(xout < 0.05)
    disp('Bad Layer Thickness: invalid results')
end         
%% SAVE layer structure, IFOmodel, and noise
costOut = getCost_aLIGO_ETM(xout, NUMTOOLS, 1);
tnowstr = datestr(now, 'yymmdd_HHMM');
savename = ['Data/' NUMTOOLS.opt_name '_' NUMTOOLS.coatingType '_' num2str(no_of_stacks) '_layers_' tnowstr];
save(savename, 'costOut', 'ifo');
