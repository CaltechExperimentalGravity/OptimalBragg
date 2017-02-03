% Code for optimizing the PR3/SR3 coating for the 40m
% https://dcc.ligo.org/E1700016

addpath('../');
%addpath(genpath('../../gwincDev'));  % add GWINCdev path to get TO coating noise

% Initial guess vector of layer thicknesses in units of lambda
no_of_stacks = 19;  % 7 pairs for 7 ppm

x0 = [];
for kk = 1:no_of_stacks
    x0 = [x0 2/8 2/8];
end
% last layer should be SiO2, since the substrate is Si
% -- this makes the total number of layers ODD
%x0 = [x0 1/4];
x0 = x0(:);                        % Make it a column

NUMTOOLS.lambda = 1064e-9;

NUMTOOLS.func         = 'getMirrorCost';

NUMTOOLS.opt_name     = 'PR3';
NUMTOOLS.T_1            = 1e-6;  % less than 50 ppm at 1064 nm
NUMTOOLS.T_2            = 1 - 1e-6;      % more than 99%    at  532 nm

ifo = SilicaTantala300;
ifo.Laser.Wavelength = NUMTOOLS.lambda;

% load this lookup table for speedup
load ../besselzeros.mat
ifo.Constants.BesselZeros = besselzeros; 
NUMTOOLS.ifo = ifo;
%NUMTOOLS.ifo = precompIFO(aSiModel, 'gwinc');

NUMTOOLS.wBeam        = 5e-3;
%NUMTOOLS.wBeam = 291e-6;     % for little Reference Cavities
NUMTOOLS.f_optimize   = 100;
NUMTOOLS.noise_weight = 1e42;  % to equate Brownian and TO noise
ifo.Optics.ETM.BeamRadius = NUMTOOLS.wBeam;
ifo.Optics.ITM.BeamRadius = NUMTOOLS.wBeam;

%% Do the optimization
% setting the bounds for the variables to be searched over
LB = 0.010 * ones(size(x0));     % bounds on the layer thicknesses
UB = 0.510 * ones(size(x0));     %              "
nvars = length(UB);

% set particle swarm optimization options
hybridopts = optimoptions('fmincon',...
                  'Display','iter',...
                  'MaxIter', 611,...
                  'TolFun', 1e-3,...
                  'MaxFunEvals', 2111);

options = optimoptions('particleswarm',...
               'SwarmSize', nvars*55,...     % 
               'UseParallel', 1,...
               'MaxIter', 211,...
               'SelfAdjustment',   1.49,...
               'SocialAdjustment', 1.49,...
               'TolFun', 1e-2,...
               'Display', 'iter',...
               'HybridFcn',{@fmincon, hybridopts});


% RUNS the Particle Swarm ========
tic
[xout, fval, exitflag] =...
    particleswarm(@(x) getMirrorCost([x], NUMTOOLS, 0),...
                                    nvars, LB, UB, options);
toc

% Check for thin layers
if find(xout < 0.001)
    disp('Bad Layer Thickness: invalid results')
end         
%% Run it one final time with the flag option turned on
TNout = getMirrorCost(xout, NUMTOOLS, 1);
tnowstr = datestr(now, 'yymmdd_HHMM');

%SbrZ = getCoatBrownian(100, ifo, ifo.Optics.ETM.BeamRadius, xout);

%disp(['Brownian displacement noise = '...
%    num2str(sqrt(SbrZ)*1e21) ' zm/rHz @ 100 Hz'])

%% SAVE layer structure, IFOmodel, and noise
savename = ['Data/' NUMTOOLS.opt_name '_layers_' tnowstr];
save(savename, 'TNout', 'ifo');

%% printing
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
