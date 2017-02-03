%% Optimization code for 40m ETM HR
%
% Rana Dec 16, 2009
%
%
global NUMTOOLS T_1064 T_532

addpath(genpath('../../lib'));

layer_matrix = [];  % Array of layer thicknesses

% Make sure to start matlabpool before doing parfor
%matlabpool(4)

% Initial guess vector of layer thicknesses in units of lambda
no_of_stacks = 17;
x0 = [1/2 7/32];
for kk = 1:no_of_stacks
    x0 = [x0 9/32 7/32];
end
x0 = x0(:);

NUMTOOLS.func      = 'optETM';
NUMTOOLS.maxiter   = 150;
NUMTOOLS.maxerr    = 0.01;
NUMTOOLS.stepsize  = 3e-3 * ones(size(x0));

N = 10000;  % number of trials
xs = zeros(length(x0),N); % matrix of layer thicknesses
errs = zeros(1,N);        % array for the scores

ys = cell(N,1);           % cell array of structs for optETM


%% Do many iterations with different initial guesses
tic
for kk = 1:N
    
    xs(:,kk) = x0 + randn(size(x0))/100;    % Add a small random perturb

    % Do the optimization
    xs(:,kk) = conjGradSolve(xs(:,kk));

    % Run it one final time with the flag option turned on
    ys{kk} = optETM(xs(:,kk),1);


    % This is a different quantity than the 'err'
    % which is returned by conjGradSolve
    err = 3 * abs(ys{kk}.T1(2) - T_1064)/T_1064 +...
          1 * abs(ys{kk}.T1(1) - T_532) / T_532
    errs(kk) = err;
    disp(['Iteration # ' num2str(kk)'])       
    drawnow    
end
toc

%% collect solutions
for jj = 1:N
    if errs(jj) < 2
        layer_matrix = [layer_matrix xs(:,jj)];
    end
end

gak = layer_matrix;

save etm_layer_matrix.txt gak -ascii -double


