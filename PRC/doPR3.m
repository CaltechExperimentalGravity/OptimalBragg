% Code for optimizing the PR3/SR3 coating for the 40m
% https://dcc.ligo.org/E1700016
% Adapted from RXA code by GV, April 2017
clc
clear
close all
addpath('../');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    USER CONFIG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
plotFlag = 1;
saveFlag = 1;

% Initial guess vector of layer thicknesses in units of lambda
no_of_stacks = 3;  
x0 = [];
for kk = 1:no_of_stacks
    x0 = [x0 2/8 2/8];
end

NUMTOOLS.lambda = 1064e-9;            % Wavelength at which R and T are to be evaluated
NUMTOOLS.opt_name       = 'PR3';
NUMTOOLS.coatingType    = 'HR';
NUMTOOLS.T_1            = 50e-6;      % Target transmission at 1064 nm
NUMTOOLS.T_2s           = 0.99;      % Target transmission at 532 nm, s-pol
NUMTOOLS.T_2p           = 0.99;      % Target transmission at 532 nm, p-pol
NUMTOOLS.aoi            = 41.1;       % Angle of incidence (degrees)
NUMTOOLS.firstLayer     = 'SiO2';     %Depending whether coating is on the HR side or AR side
NUMTOOLS.wBeam        = 5e-3;          % Beam size on optic being optimized
NUMTOOLS.f_optimize   = 100;           % Frequency at which to evaluate noise being optimized
NUMTOOLS.noise_weight = 1e42;          % to equate Brownian and TO noise
NUMTOOLS.include_sens   = 1;           %Include derivatives in sensitivity function (1 or 0)
NUMTOOLS.include_thermal   = 1;        %Include TN term in sensitivity function (1 or 0)


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
UB = 0.510 * ones(size(x0));     % upper bound, lambda/2          
nvars = length(UB);

if strcmp(NUMTOOLS.coatingType,'HR')
    LB(1) = 0.48;  %Constrain topost layer to be nearly halfwave thick...
    UB(1) = 0.52;
end

% set particle swarm optimization options
hybridopts = optimoptions('fmincon',...
                  'Display','iter',...
                  'MaxIter', 911,...
                  'TolFun', 1e-3,...
                  'MaxFunEvals', 11111);

options = optimoptions('particleswarm',...
               'SwarmSize', nvars*155,...     % 
               'UseParallel', 1,...
               'MaxIter', 911,...
               'SelfAdjustment',   1.49,...
               'SocialAdjustment', 1.49,...
               'TolFun', 1e-3,...
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
savename = ['Data/' NUMTOOLS.opt_name '_layers_' tnowstr];
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
