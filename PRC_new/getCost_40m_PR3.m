function y = getCost_40m_PR3(x, params, flag)
% getCost_40m_PR3
% This function evaluates a SCALAR cost function which is minimized 
% by the particle swarm optimizing script runSwarm.m
% x ------ initial guess for dielectric stack layer thicknesses
% params - structure with various parameters for evaluating the cost function
% flag --- telling the function whether to save data or not
% terms in the cost function are taking into consideration the design
% requirements on the 40m PR3, see https://dcc.ligo.org/E1800089.

%Importing the parameter called "NUMTOOLS" in do<Optic>.m
glob_param = params;

%Making this a row vector for multidiel1.m
L = x(:)';

ifo = glob_param.ifo;

% Laser wavelength
lambda_0 = glob_param.lambda;
% Angle of incidence
aoi = glob_param.aoi;

%Dispersion taken from Ramin
na = 1.0000;           %Vacuum
n1_IR = glob_param.n1_IR;
n2_IR = glob_param.n2_IR;
nb_IR = glob_param.nb_IR;

%Weight vector for scalarization of cost function
weights = glob_param.wt;

no_of_stacks = floor(length(x)/2);    % use floor if x is not even
% set up the (alternating) array of index of refractions

if strcmp(glob_param.coatingType, 'HR') %HR coating
    n_c = [];
    for kk = 1:(no_of_stacks)
        n_c  = [n_c n1_IR n2_IR];
    end
    n_IR = [na n_c nb_IR];  % add the indices of the vacuum and the substrate at either end
end

% list of wavelengths to use for the optimization
lambda_IR = [lambda_0] / lambda_0;

% Calculate reflectivities for each polarization at the angle of incidence
% ['tm' = 'p-pol', 'te' = 's-pol'], for the given layer structure
[Gammap, ~]   = multidiel1(n_IR, L, lambda_IR, aoi, 'tm');
Rp_IR = abs(Gammap).^2;
Tp_IR = 1 - Rp_IR;

% minimize reflected E-field (usually by 1/2 wave cap on top)
% only at main wavelength
r_refl = abs(1 + Gammap);

% define the error function
T_1 = glob_param.T_1;         % desired  transmission @ lambda_0
R_1 = 1 - T_1;

%surfTarg = glob_param.surfaceField;

% cost function which gets minimized 
yy = [];
% Indices of the yy term
% 1 = Transmission at 1064, p-pol
% 2 = sensitivity to change in layer thickness @1064nm p- +/-1%...
% 3 = sensitivity to change in n1 @1064nm p- +/-1%...
% 4 = sensitivity to change in n2 @1064nm p- +/-1%...
% 5 = sensitivity to change in aoi +/-5%...

%Don't penalize for being better than the spec...
%cost_T_IR = heaviside(Tp_IR - T_1) * ((Tp_IR - T_1)/T_1)^2;
cost_T_IR = (Tp_IR/T_1)^2;
if cost_T_IR > 20
	yy = [yy (20./log10(20))*log10(cost_T_IR)];    % match the T @ lambda
else
	yy = [yy cost_T_IR];
end

y = sum(yy(1) .* weights(1));   %Terms included by default...


if glob_param.include_sens
    [err_IR] = doSens_40m_PR3(n_IR, L, lambda_IR, aoi, Tp_IR, Rp_IR);
    %%Add sensitivity function to the cost function...
    e1 = 1*abs(err_IR.Tp.coatLayer_plus) + 1*abs(err_IR.Tp.coatLayer_minus);
    if e1>20
	    yy = [yy (20./log10(20))*log10(e1)];
    else
	    yy = [yy e1];
    end 

    e1 = 1*abs(err_IR.Tp.n1_plus) + 1*abs(err_IR.Tp.n1_minus);
    if e1>20
	    yy = [yy (20./log10(20))*log10(e1)];
    else
	    yy = [yy e1];
    end

    e1 = 1*abs(err_IR.Tp.n2_plus) + 1*abs(err_IR.Tp.n2_minus);
    if e1>20
	    yy = [yy (20./log10(20))*log10(e1)];
    else
	    yy = [yy e1];
    end
    
    e1 = 1*abs(err_IR.Tp.aoi_plus) + 1*abs(err_IR.Tp.aoi_minus);
    if e1>20
	    yy = [yy (20./log10(20))*log10(e1)];
    else
	    yy = [yy e1];
    end
    y = y + sum(yy(2:5) .* weights(2:5));
end

if flag == 1 %Make an output structure with all variables of importance..
    clear y;
    y.aoi = aoi;
    y.L = L;
    y.Tp_IR = Tp_IR;
    y.Rp_IR = Rp_IR;
    y.yy = yy;
    y.n_IR = n_IR;
    y.weights = weights;
end

end

