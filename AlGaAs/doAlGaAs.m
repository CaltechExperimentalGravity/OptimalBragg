% Example code for coating optimization
% 

addpath(genpath('../'));
addpath(genpath('../generic/'));
%addpath(genpath('../../gwincDev'));  % add GWINCdev path to get TO coating noise

% Initial guess vector of layer thicknesses in units of lambda
no_of_stacks = 40;  

x0 = [];

x0 = [x0; ones(2*no_of_stacks,1)];

% add a layer of GaAs at the bottom
%x0 = [x0; 1/4];
x0 = x0(:);                        % Make it a column

NUMTOOLS.lambda = 1064e-9;

NUMTOOLS.func         = 'getMirrorCost';

NUMTOOLS.opt_name     = 'ETM';
NUMTOOLS.T            = 5e-6;

ifo = AlGaAsModel;
ifo.Laser.Wavelength = NUMTOOLS.lambda;
ifo.Materials.Substrate.Temp = 293;

% load this lookup table for speedup
load ../besselzeros.mat
ifo.Constants.BesselZeros = besselzeros;
NUMTOOLS.ifo = ifo;

NUMTOOLS.wBeam        = 0.06;
%NUMTOOLS.wBeam = 291e-6;     % for little Reference Cavities
NUMTOOLS.f_optimize   = 100;
ifo.Optics.ETM.BeamRadius = NUMTOOLS.wBeam;
ifo.Optics.ITM.BeamRadius = NUMTOOLS.wBeam;

%% Do the optimization
% setting the bounds for the variables to be searched over
LB = 0.002 * ones(size(x0));     % bounds on the layer thicknesses
UB = 0.500 * ones(size(x0));     %              "
nvars = length(UB);

% set particle swarm optimization options
hybridopts = optimoptions('fmincon',...
                  'Display','iter',...
                  'MaxIter', 611,...
                  'TolFun', 1e-3,...
                  'MaxFunEvals', 2111);

options = optimoptions('particleswarm',...
               'SwarmSize', nvars*154,...     % 
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
    particleswarm(@(x) getMirrorCost(x, NUMTOOLS, 0),...
                                    nvars, LB, UB, options);
toc

% Check for thin layers
if find(xout < 0.001)
    disp('Bad Layer Thickness: invalid results')
end

%% SAVE layer structure, IFOmodel, and noise
savename = ['Data/' NUMTOOLS.opt_name '_layers_' tnowstr];
save(savename, 'TNout', 'ifo');

%% Run it one final time with the flag option turned on
TNout = getMirrorCost(xout, NUMTOOLS, 1);
tnowstr = datestr(now, 'yymmdd_HHMM');

SbrZ = getCoatBrownian(100, ifo, ifo.Optics.ETM.BeamRadius, xout);

disp(['Brownian displacement noise = '...
    num2str(sqrt(SbrZ)*1e21) ' zm/rHz @ 100 Hz'])


%% printing
title(' ')
orient landscape
set(gcf,'PaperPositionMode','auto')

fname = ['Figures/' NUMTOOLS.opt_name '_R_' tnowstr];

print('-depsc','-r600', fname)
[a,b] = system(['Figures/makePDF.sh ' fname '.eps']);
if a ~= 0
    disp('PDF Generation Error')
else
    [a,b] = system(['rm ' fname '.eps']);
end

%% plot layer structure and thermo-optic noise
plotTOnoise
