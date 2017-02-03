% Example code for coating optimization
% currently uses optimal ETM coating,
global NUMTOOLS

addpath(genpath('../../lib'));

% Initial guess vector of layer thicknesses in units of lambda
no_of_stacks = 8;
x0 = [1/2 1/4];
for kk = 1:no_of_stacks
    x0 = [x0 3/8 1/8];
end
x0 = x0 + rand(size(x0))/1500;    % Add a small random perturb
x0 = x0(:);                       % Make it a column

NUMTOOLS.func      = 'optITM';
NUMTOOLS.maxiter   = 5150;
NUMTOOLS.maxerr    = 1e-6;
NUMTOOLS.stepsize  = 1e-3 * ones(size(x0));

% Do the optimization
xout = conjGradSolve(x0);

% Run it one final time with the flag option turned on
ye = optITM(xout,1);

gak = ye.L';

save itm_layers.txt gak -ascii -double


