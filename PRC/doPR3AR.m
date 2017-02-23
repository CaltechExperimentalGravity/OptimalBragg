% Code for optimizing the PR3/SR3 coating for the 40m
% https://dcc.ligo.org/E1700016
% Adapted from RA code by GV, Feb 2017

addpath('../');

% Initial guess vector of layer thicknesses in units of lambda
no_of_stacks = 19;  
x0 = [];
for kk = 1:no_of_stacks
    x0 = [x0 2/8 2/8];
end
% last layer should be SiO2, since the substrate is SiO2

%Make the layer structure a column...
x0 = x0(:);                       

%Laser wavelength
NUMTOOLS.lambda = 1064e-9;            % Wavelength at which R and T are to be evaluated

% NUMTOOLS.func         = 'getMirrorCost_AR';

NUMTOOLS.opt_name       = 'PR3AR';
NUMTOOLS.T_1            = 1. - 50e-6;      % Target transmission at 1064 nm
NUMTOOLS.T_2s           = 1. - 50e-6;      % Target transmission at 532 nm, s-pol
NUMTOOLS.T_2p           = 1. - 50e-6;      % Target transmission at 532 nm, p-pol
NUMTOOLS.aoi            = 24.9;       % Angle of incidence (degrees), assuming 2 degree wedge


ifo = SilicaTantala300;                 % Load a standard gwinc style parameter file with various material properties defined             
ifo.Laser.Wavelength = NUMTOOLS.lambda;

% load this lookup table for speedup
load ../besselzeros.mat
ifo.Constants.BesselZeros = besselzeros; 
NUMTOOLS.ifo = ifo;
%NUMTOOLS.ifo = precompIFO(aSiModel, 'gwinc');

NUMTOOLS.wBeam        = 5e-3;          % Beam size on optic being optimized
NUMTOOLS.f_optimize   = 100;           % Frequency at which to evaluate noise being optimized
NUMTOOLS.noise_weight = 1e42;          % to equate Brownian and TO noise
ifo.Optics.ETM.BeamRadius = NUMTOOLS.wBeam; %using the gwinc infrastructure for calculating TO noise for ease...
ifo.Optics.ITM.BeamRadius = NUMTOOLS.wBeam;

%% Do the optimization
% setting the bounds for the variables to be searched over
LB = 0.010 * ones(size(x0));     % lower bounds on the layer thicknesses, lambda/100
UB = 0.510 * ones(size(x0));     % upper bound, lambda/2          
nvars = length(UB);

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


% RUNS the Particle Swarm ========
tic
[xout, fval, exitflag] =...
    particleswarm(@(x) getMirrorCost_AR([x], NUMTOOLS, 0),...
                                    nvars, LB, UB, options);
toc

% Check for thin layers
if find(xout < 0.001)
    disp('Bad Layer Thickness: invalid results')
end         
%% Run it one final time with the flag option turned on
TNout = getMirrorCost_AR(xout, NUMTOOLS, 1);
tnowstr = datestr(now, 'yymmdd_HHMM');

%SbrZ = getCoatBrownian(100, ifo, ifo.Optics.ETM.BeamRadius, xout);

%disp(['Brownian displacement noise = '...
%    num2str(sqrt(SbrZ)*1e21) ' zm/rHz @ 100 Hz'])

%% SAVE layer structure, IFOmodel, and noise
savename = ['Data/' NUMTOOLS.opt_name '_layers_' tnowstr];
save(savename, 'TNout', 'ifo');

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
%plotTO120
% exit;
