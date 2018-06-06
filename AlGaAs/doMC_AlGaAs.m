function doMC_AlGaAs(filename, N, nDim, savename)

% OMG - I need a decent help comment + examples !!!! Aaahhhhhh!

% Function to perturb a given coating design to see how sensitive 
% derived parameters like Transmission etc depend on perturbations
% Since perturbations are assumed to be iid, no fancy MCMC sampler is used
% If we desire, we can adapt this code to use 
%   https://www.mathworks.com/matlabcentral/fileexchange/49820-ensemble-mcmc-sampler

% this is a hacky way of doing matlab coding
clc
close all

try % set path if its not there already
    op2phys(4, 2);
catch
    addpath(genpath('../'));
    addpath(genpath('../generic/'));
    %addpath(genpath('../../gwincDev'));  % add GWINCdev path to get TO coating noise
end

% Load the coating file...
if filename == 0
    % get time str for latest coating design
    d = dir('Data/*layers*.mat');
    [~,idx] = max([d.datenum]); % get index of newest file
    filename = d(idx).name;
    strfind(filename,'_18');
    tnowstr = filename(12:22)
    if strfind(d(idx).name,'ETM')
        NUMTOOLS.opt_name = 'ETM';
    else
        NUMTOOLS.opt_name = 'ITM';
    end
    funame = [NUMTOOLS.opt_name '_layers_' tnowstr];
    filename = ['Data/' funame];
end

load(filename);

% number of output variables
nVars = 4;

% Generate the perturbations...
means    = zeros(nDim, 1);
sigmas   = diag(diag((0.005^2)*ones(nDim), 0));
perturbs = mvnrnd(means, sigmas, N);

% Define arrays for output variables of interest
T_IR      = zeros(N,1);
surfField = zeros(N,1);
TOnoise   = zeros(N,1);
BRnoise   = zeros(N,1);

% Nominal values from coating design
aoi = 0;
n   = TNout.n;
L   = TNout.L;
L_phys = op2phys(L, n(2:end-1));

% Nominal params for calculation
Ei    = 27.46; % [V/m], for surface field calculation ???
f_to  = 100;   % [Hz]
wBeam = 0.065; % [meters]

%%%%%   Apply the perturbations   %%%%%%%%
for i = 1:N
    n_IRs   = n;
    Ls      = L_phys;
    perturb = 1 + perturbs(i,:);
    % All physical thicknesses with equal error
    Ls = Ls * perturb(1);   
    % Error in refractive index of low index layers
    n_IRs(2:2:end-1) = n_IRs(2:2:end-1) .* perturb(2);
    % Error in refractive index of high index layers
    n_IRs(3:2:end-1) = n_IRs(3:2:end-1) .* perturb(3);

    % Compute reflectivity, surface field, Brownian noise and TO noise.
    [Gamma, ~] = multidiel1(n_IRs, Ls.*n_IRs(2:end-1), TNout.lambda);
    T_IR(i)    = 1e6*(1 - abs(Gamma).^2);
    % Thermo-Optic 
    ifo.Materials.Coating.Indices = n_IRs;
    [StoZ, ~, ~, ~]  = getCoatThermoOptic(f_to, ifo, wBeam,...
                                           Ls'.*n_IRs(2:end-1)');
    TOnoise(i) = sqrt(StoZ);
    % Brownian
    SbrZ = getCoatBrownian(f_to, ifo, wBeam, Ls.*n_IRs(2:end-1));
    BRnoise(i) = sqrt(SbrZ);
    % Surface Field
    surfField(i) = Ei * abs(1+Gamma);
end

% Save everything for corner plotting with python 
% (why not use regular save command? Its HDF5 naturally)
MCout = horzcat(T_IR, 1e21*TOnoise, 1e21*BRnoise, surfField).';

disp(['Saving into ' funame])
save(['MCout/' funame], 'MCout')

if savename == 0
    savename = ['MCout/' funame '.h5'];
end
system(['rm ' savename]);
h5create(savename, '/MCout', [nVars N]);
h5write( savename, '/MCout', MCout);

end

