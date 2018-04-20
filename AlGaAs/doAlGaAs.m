% PSO for AlGaAs coating design for ETM/ITM
%

try % set path if its not there already
    op2phys(4, 2);
catch
    addpath(genpath('../'));
    addpath(genpath('../generic/'));
    %addpath(genpath('../../gwincDev'));  % add GWINCdev path to get TO coating noise
end

% Initial guess vector of layer thicknesses in units of lambda
no_of_pairs = 44;
x0 = [];
x0 = [x0; 0.25*ones(2*no_of_pairs,1)];

% add a layer of GaAs at the bottom
x0 = [x0; 1/4];
x0 = x0(:);                        % Make it a column

NUMTOOLS.lambda    = 1064e-9;
NUMTOOLS.func      = 'getMirrorCost';
NUMTOOLS.opt_name  = 'ETM';
NUMTOOLS.T         = 5e-6;  % desired power transmission

ifo = AlGaAsModel;
ifo.Laser.Wavelength = NUMTOOLS.lambda;
ifo.Materials.Substrate.Temp = 293;
ifo.Materials.Coating.Type = 'AlGaAs';  % high index (GaAs) cap layer

% load this lookup table for speedup
load ../besselzeros.mat
ifo.Constants.BesselZeros = besselzeros;
NUMTOOLS.ifo = ifo;

NUMTOOLS.wBeam        = 0.065;
%NUMTOOLS.wBeam = 291e-6;     % for little Reference Cavities
NUMTOOLS.f_optimize   = 100;  % minimize the noise at this
                              % frequency only
ifo.Optics.ETM.BeamRadius = NUMTOOLS.wBeam;
ifo.Optics.ITM.BeamRadius = NUMTOOLS.wBeam;

%% Do the optimization
% setting the bounds for the variables to be searched over
minThick = 0.002;
maxThick = 0.600;
LB = minThick * ones(size(x0));     % bounds on the layer thicknesses
UB = maxThick * ones(size(x0));     %              "
nvars = length(UB);

nswarms = 19;
x0 = 0.25 + 0.1*(rand(nswarms, nvars) *(maxThick-minThick) + minThick);

% set particle swarm optimization options
hybridopts = optimoptions('fmincon',...
                  'Display','iter',...
                  'MaxIter', 611,...
                  'TolFun', 1e-2,...
                  'MaxFunEvals', 21111);

options = optimoptions('particleswarm',...
                       'InitialSwarmMatrix', x0,...
                       'SwarmSize', nvars*nswarms,...
                       'UseParallel', 1,...
                       'MaxIter', 211,...
                       'SelfAdjustment',   0.49,...
                       'SocialAdjustment', 0.49,...
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

%% Run it one final time with the flag option turned on
TNout = getMirrorCost(xout, NUMTOOLS, 1);
tnowstr = datestr(now, 'yymmdd_HHMM');

%% SAVE layer structure, IFOmodel, and noise
savename = ['Data/' NUMTOOLS.opt_name '_layers_' tnowstr];
save(savename, 'TNout', 'ifo');



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
