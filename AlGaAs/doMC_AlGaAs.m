function doMC_AlGaAs(varargin)

% Function to perturb a given coating design to see how sensitive
% derived parameters like Transmission etc depend on perturbations
% Since perturbations are assumed to be iid, no fancy MCMC sampler is used
% If we desire, we can adapt this code to use
% https://www.mathworks.com/matlabcentral/fileexchange/49820-ensemble-mcmc-sampler
% Input arguments:
%   - filename = Path to .mat file that is output from PSO optimization, for which MC calculation is to be done
%   - N        = number of MC samples to generate
%   - nDim     = number of variables being perturbed in this MC study
%   - savename = Path to .hdf5 file to which the MC output is to be saved. This will be used for nice corner plotting with Python
% Outputs: NONE
%
% Example usage:
%	doMC_AlGaAs('Data/ETM_layers_180607_0028.mat', 1e5, 5, 'Data/ETM_layers_180607_0028.hdf5')


if nargin == 0
    filename = 0;
    N        = 1e5;
    nDim     = 5;
    savename = 0;
elseif nargin == 4
    filename = varargin{1};
    N        = varargin{2};
    nDim     = varargin{3};
    savename = varargin{4};
else
    warning('Illegal Number of INput Args !!')
end


% this is a hacky way of doing matlab coding
%clc
%close all

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
nVars = 5;

% Generate the perturbations...
means    = zeros(nDim, 1);
sigmas   = diag(diag((0.005^2)*ones(nDim), 0));
perturbs = mvnrnd(means, sigmas, N);

% Define arrays for output variables of interest
T_IR      = zeros(N,1);
surfField = zeros(N,1);
TOnoise   = zeros(N,1);
BRnoise   = zeros(N,1);
absorp    = zeros(N,1);

% Nominal values from coating design
aoi    = 0;  % the beam is at normal incidence
n      = TNout.n;
L      = TNout.L;
L_phys = op2phys(L, n(2:end-1));

zifo = ifo;  % weird that this is necessary for parfor to use this variable

% Nominal params for calculation
% Corresponds to 1 W/m^2 peak intensity incident Gaussian beam. 
Ei     = 27.46;      % [V/m], for surface field calculation. 
f_to   = 100;        % [Hz]
wBeam  = 0.065;      % [meters]
lam    = 1064e-9;    % [m], laser wavelength
nPts   = 10;         % [m], number of points inside each layer at which to evaluate E field squared
alpha_GaAs   = 1.5;  % [m^-1], Absorption of GaAs   layers
alpha_AlGaAs = 4.5;  % [m^-1], Absorption of AlGaAs layers

%reverseStr = '';
%kk = 0;

%%%%%   Apply the perturbations   %%%%%%%%
parfor i = 1:N

    % Display the progress
% $$$     if  mod(kk, 100) < 1                
% $$$       percentDone = 100 * i / N;
% $$$       msg = sprintf('Percent done: %3.1f', percentDone);
% $$$       fprintf([reverseStr, msg]);
% $$$       reverseStr = repmat(sprintf('\b'), 1, length(msg));
% $$$     end
    
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
    [Gamma, ~] = multidiel1(n_IRs, Ls .* n_IRs(2:end-1), lam);
    T_IR(i)    = 1e6*(1 - abs(Gamma).^2);

    % Thermo-Optic
    [StoZ, ~, ~, ~]  = getCoatThermoOptic(f_to,...
                          myfoe(zifo, n_IRs), wBeam, Ls.' .* n_IRs(2:end-1).');
    TOnoise(i) = sqrt(StoZ);
    
    % Brownian
    SbrZ = getCoatBrownian(f_to,...
               myfoe(zifo, n_IRs), wBeam, Ls .* n_IRs(2:end-1));
    BRnoise(i) = sqrt(SbrZ);

    % Surface Field
    surfField(i) = Ei * abs(1+Gamma);

    % Absorption
    [zz, E_prof] = calcEField_AlGaAs(lam*Ls, n_IRs,...
                                     length(Ls), lam, 0,'p',nPts);
    absorp(i) = calcAbsorption_AlGaAs(E_prof,...
                          lam*Ls, nPts, alpha_GaAs, alpha_AlGaAs);
end

% Save everything for corner plotting with python
% (why not use regular save command? Its HDF5 naturally)
MCout = horzcat(T_IR, 1e21*TOnoise, 1e21*BRnoise, surfField, absorp).';

disp('  ')
disp(['Saving into ' funame])
%save(['MCout/' funame], 'MCout')

if savename == 0
    savename = ['MCout/' funame '.h5'];
end
[a,b] = system(['rm ' savename]);
h5create(savename, '/MCout', [nVars N]);
h5write( savename, '/MCout', MCout);

end


function ifo = myfoe(ifo, n_IRs)
        ifo.Materials.Coating.Indices = n_IRs;
end
