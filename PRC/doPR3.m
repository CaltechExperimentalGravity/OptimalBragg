% Code for optimizing the PR3/SR3 coating for the 40m
% https://dcc.ligo.org/E1700016
% Adapted from RXA code by GV, April 2017
clc
% clear
close all
addpath('../');

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


plotFlag = 1;
saveFlag = 1;

% Initial guess vector of layer thicknesses in units of lambda
%no_of_stacks = 20;  
% no_of_stacks = 6;  
x0 = [];
for kk = 1:no_of_stacks
    x0 = [x0 2/8 2/8];
end
 
%load('Data/170623_HR/PR3_HR_20_layers_170622_1451.mat');             
% load('Data/170621_AR/PR3_layers_170621_1911.mat');             
% x0 = TNout.L;                                                               
clear TNout                                                                 
clear ifo    


NUMTOOLS.lambda = 1064e-9;            % Wavelength at which R and T are to be evaluated
NUMTOOLS.opt_name       = 'PR3';
NUMTOOLS.coatingType    = 'AR';
NUMTOOLS.T_1            = 0.99;      % Target transmission at 1064 nm --- 50ppm for HR, 99% for AR
NUMTOOLS.T_2s           = 0.998;      % Target transmission at 532 nm, s-pol --- 99.5% for HR, 99.8% for AR
NUMTOOLS.T_2p           = 0.998;      % Target transmission at 532 nm, p-pol ---- same as above
NUMTOOLS.aoi            = 24.967;       % Angle of incidence (degrees)  --- 41.1deg for HR, 24.8deg for AR
NUMTOOLS.aoi_green      = 24.746;       % Angle of incidence (degrees)  ---24.746 deg for green
NUMTOOLS.firstLayer     = 'Ta2O5';     %Depending whether coating is on the HR side or AR side
NUMTOOLS.wBeam        = 5e-3;          % Beam size on optic being optimized
NUMTOOLS.f_optimize   = 100;           % Frequency at which to evaluate noise being optimized
NUMTOOLS.noise_weight = 1e42;          % to equate Brownian and TO noise
NUMTOOLS.include_sens   = 1;           %Include derivatives in sensitivity function (1 or 0)
NUMTOOLS.include_thermal   = 0;        %Include TN term in sensitivity function (1 for HR or 0 for AR)


ifo = SilicaTantala300;                 % Load a standard gwinc style parameter file with various material properties defined             
ifo.Laser.Wavelength = NUMTOOLS.lambda;
% load this lookup table for speedup
load ../besselzeros.mat
ifo.Constants.BesselZeros = besselzeros;
ifo.Optics.ETM.BeamRadius = NUMTOOLS.wBeam; %using the gwinc infrastructure for calculating TO noise for ease...
ifo.Optics.ITM.BeamRadius = NUMTOOLS.wBeam;
NUMTOOLS.ifo = ifo;

% setting the bounds for the variables to be searched over
LB = 0.050 * ones(size(x0));     % lower bounds on the layer thicknesses, lambda/50
UB = 0.51 * ones(size(x0));     % upper bound, lambda/2          
nvars = length(UB);

if strcmp(NUMTOOLS.coatingType,'HR')
    LB(1) = 0.48;  %Constrain topost layer to be nearly halfwave thick...
    UB(1) = 0.52;
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
TNout = getMirrorCost(xout, NUMTOOLS, 1);
tnowstr = datestr(now, 'yymmdd_HHMM');
savename = ['Data/' NUMTOOLS.opt_name '_' NUMTOOLS.coatingType '_' num2str(no_of_stacks) '_layers_' tnowstr];
save(savename, 'TNout', 'ifo');

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
