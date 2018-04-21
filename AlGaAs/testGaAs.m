% PSO for AlGaAs coating design for ETM/ITM
%

try % set path if its not there already
    op2phys(4, 2);
catch
    addpath(genpath('../'));
    addpath(genpath('../generic/'));
end

% Initial guess vector of layer thicknesses in units of lambda
%no_of_pairs = 9;
%x0 = [];
%x0 = [x0; 0.25*ones(2*no_of_pairs,1)];
x0 =    [0.1896 0.1121 0.4995 0.1000 0.4598];
x0 = [x0 0.1695 0.2760 0.2145 0.2510 0.2388];
x0 = [x0 0.2403 0.2508 0.2368 0.2553 0.2375];
x0 = [x0 0.2571 0.2391 0.2564 0.2414 0.2550];
x0 = [x0 0.2437 0.2533 0.2459 0.2515 0.2480];
x0 = [x0 0.2498 0.2498 0.2482 0.2514 0.2469];
x0 = [x0 0.2528 0.2457 0.2539 0.2447 0.2549];
x0 = [x0 0.2439 0.2556 0.2433 0.2562 0.2427];
x0 = [x0 0.2566 0.2423 0.2571 0.2420 0.2574];
x0 = [x0 0.2417 0.2577 0.2414 0.2579 0.2412];
x0 = [x0 0.2581 0.2409 0.2585 0.2405 0.2587];
x0 = [x0 0.2401 0.2556];
% add a layer of GaAs at the bottom
%x0 = [x0; 1/4];
xout = x0(:);                        % Make it a column

NUMTOOLS.lambda    = 1064e-9;
NUMTOOLS.func      = 'getMirrorCost';
NUMTOOLS.opt_name  = 'ETM';
NUMTOOLS.T         = 100e-6;  % desired power transmission

ifo = AlGaAsModel;
ifo.Laser.Wavelength = NUMTOOLS.lambda;
ifo.Materials.Substrate.Temp = 305;
ifo.Materials.Coating.Type = 'AlGaAs';  % high index (GaAs) cap layer

% load this lookup table for speedup
load ../besselzeros.mat
ifo.Constants.BesselZeros = besselzeros;
NUMTOOLS.ifo = ifo;

%NUMTOOLS.wBeam     = 0.065;
NUMTOOLS.wBeam      = 215e-6;     % for little Reference Cavities
NUMTOOLS.f_optimize = 100;  % minimize the noise at this frequency

ifo.Optics.ETM.BeamRadius = NUMTOOLS.wBeam;
ifo.Optics.ITM.BeamRadius = NUMTOOLS.wBeam;

%% 
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
