% Example code for coating optimization
% currently uses optimal ETM coating,
global NUMTOOLS

addpath(genpath('../../lib'));

% Initial guess vector of layer thicknesses in units of lambda
no_of_stacks = 17;
x0 = [1/2 7/32];
for kk = 1:no_of_stacks
    x0 = [x0 8/32 8/32];
end
x0 = x0 + randn(size(x0))/2000;    % Add a small random perturb
x0 = x0(:);                      % Make it a column

NUMTOOLS.func      = 'optETM';
NUMTOOLS.maxiter   = 150;
NUMTOOLS.maxerr    = .001;
NUMTOOLS.stepsize  = 1e-4 * ones(size(x0));

% Do the optimization
x0 = conjGradSolve(x0);

% Run it one final time with the flag option turned on
ye = optETM(x0,1);

gak = ye.L';

save etm_layers.txt gak -ascii -double


