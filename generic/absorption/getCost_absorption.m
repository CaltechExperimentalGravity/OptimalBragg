function y = getCost_absorption(x, params, flag)
% getCost_absorption
% This function evaluates a SCALAR cost function which is minimized 
% by the particle swarm optimizing script runSwarm_absorption.m
% x ------ initial guess for dielectric stack layer thicknesses
% params - structure with various parameters for evaluating the cost function
% flag --- telling the function whether to save data or not
% terms in the cost function are taking into consideration the design
% requirements on the aLGIO ETM.

%Importing the parameter called "NUMTOOLS" in runSwarm_absorption.m
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
    
elseif strcmp(glob_param.coatingType, 'AR') %AR coating
    n_c = [];
    for kk = 1:(no_of_stacks)
        n_c  = [n_c n2_IR n1_IR];
    end
    n_IR = [nb_IR n_c na];  % add the indices of the vacuum and the substrate at either end
    
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

surfTarg = glob_param.surfaceField;

% cost function which gets minimized 
yy = [];
% Indices of the yy term
% 1 = Transmission at 1064, p-pol
% 2 = Absorption (integrated) at 1064nm
% 3 = HR surface field

%Don't penalize for being better than the spec...
yy = [yy heaviside(Tp_IR - T_1)*abs((Tp_IR - T_1)/T_1)^1];    % match the T @ lambda

if glob_param.minAbsorp
        [zz, EE] = calcEField(lambda_0*op2phys(L,n_IR(2:end-1)), n_IR, 10, lambda_0, 0, 'p', glob_param.nPts);
        absorp = calcAbsorption(EE, lambda_0*op2phys(L(1:10),n_IR(2:11)), glob_param.nPts, glob_param.alpha1, glob_param.alpha2);
        yy = [yy absorp];
end


if strcmp(glob_param.coatingType,'HR')
    %Term for minimizing surface E-field, see T0900626
    yy     = [yy heaviside(27.46*r_refl - surfTarg)*abs((27.46*r_refl -surfTarg) / surfTarg)^1]; %this 27.46 is assuming Ein = 1W/m^2
else
    yy = [yy 0];
end

%Don't bother computing other terms unless we are within factor of 2 of design goals
%if yy(1)<2 && yy(2)<2 && yy(3)<2
y = sum(yy(1:3) .* weights(1:3));   %Terms included by default...

if glob_param.include_brownian
    % ----Brownian Thermal noise - parameters from GWINC IFOmodel file
    % Formula from E0900068-v3, see also T0900161
    phi_high = ifo.Materials.Coating.Phihighn;
    n_high   = n2_IR; 
    Y_high   = ifo.Materials.Coating.Yhighn;

    phi_low  = ifo.Materials.Coating.Philown; 
    n_low    = n1_IR;
    Y_low    = ifo.Materials.Coating.Ylown; 

    Y_sub  = ifo.Materials.Substrate.MirrorY;

    % Total thickness
    z_low  = sum(L(1:2:length(L)));  %odd layers
    z_high = sum(L(2:2:length(L)));  %even layers

    a = phi_high / phi_low;
    b = n_low / n_high;
    c = Y_high/Y_sub + Y_sub/Y_high;
    d = Y_low/Y_sub  + Y_sub/Y_low;
    little_gamma = a*b*c/d;

    % Brownian thermal noise estimate (only used for optimization)
    %This is a kind of proxy function for minimizing coating thermal noise
    S = z_low + (little_gamma * z_high);
    % ---------------------------------------------------
    y = y + weights(4) * S;
    yy(4) = S;                           % minimize the Brownian noise
    
else
    S = 0;
    yy(4) = 0;
end

if glob_param.include_TO
    % Thermo-Optic Noise (there should be a flag to use this or not)
    f_to  = glob_param.f_optimize;     %Frequency to evaluate this
    wBeam = glob_param.wBeam;
    dOpt  = L';
    [StoZ, SteZ, StrZ, T]  = getCoatThermoOptic(f_to, ifo, wBeam, dOpt);
    TO = sqrt(StoZ); %m/rtHz
    yy(5) = TO;
    y = y + weights(5) * TO;
else
    TO = 0;
    %y = y + weights(5) * TO;
    yy(5) = TO;
end

if glob_param.include_sens
    [err_IR] = doSens(n_IR, L, lambda_IR, aoi, Tp_IR, Rp_IR);
    [err_green] = doSens(n_green, L_green, lambda_green, aoi_green, Tp_green, Rp_green);

    %%Add sensitivity function to the cost function... 10x weight for IR
    %%than green
    e1 = 1*abs(err_IR.Tp.coatLayer_plus) + 1*abs(err_IR.Tp.coatLayer_minus);
    e2 = 1*abs(err_green.Rp.coatLayer_plus) + 1*abs(err_green.Rp.coatLayer_minus);
    yy = [yy e1 e2]; %1064nmp, 532nm p-pol

    e1 = 1*abs(err_IR.Tp.n1_plus) + 1*abs(err_IR.Tp.n1_minus);
    e2 = 1*abs(err_green.Rp.n1_plus) + 1*abs(err_green.Rp.n1_minus);
    yy = [yy e1 e2];

    e1 = 1*abs(err_IR.Tp.n2_plus) + 1*abs(err_IR.Tp.n2_minus);
    e2 = 1*abs(err_green.Rp.n2_plus) + 1*abs(err_green.Rp.n2_minus);
    yy = [yy e1 e2];

    %e1 = 1*abs(err_IR.Tp.aoi_plus) + 1*abs(err_IR.Tp.aoi_minus);
    %e2 = 1*abs(err_green.Rp.aoi_plus) + 1*abs(err_green.Rp.aoi_minus);
    %yy = [yy e1 e2];
    
	yy = [yy err_IR.totSurfField];

    y = y + sum(yy(6:12) .* weights(6:12));
end

if flag == 1 %Make an output structure with all variables of importance..
    clear y;
    y.aoi = aoi;
    y.L = L;
    y.Tp_IR = Tp_IR;
    y.absorp = absorp;
    y.surfField = yy(3);
    y.yy = yy;
    y.n_IR = n_IR;
    y.weights = weights;
end

end

