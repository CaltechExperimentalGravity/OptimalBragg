function y = getCost_aLIGO_ETM(x, params, flag)
% getCost_aLIGO_ETM
% This function evaluates a SCALAR cost function which is minimized 
% by the particle swarm optimizing script runSwarm.m
% x ------ initial guess for dielectric stack layer thicknesses
% params - structure with various parameters for evaluating the cost function
% flag --- telling the function whether to save data or not
% terms in the cost function are taking into consideration the design
% requirements on the aLGIO ETM.

%Importing the parameter called "NUMTOOLS" in do<Optic>.m
glob_param = params;

%Making this a row vector for multidiel1.m
L = x(:)';

ifo = glob_param.ifo;

% Laser wavelength
lambda_0 = glob_param.lambda;
% Angle of incidence
aoi = glob_param.aoi;
aoi_green = glob_param.aoi;

%Dispersion taken from Ramin
na = 1.0000;           %Vacuum
n1_IR = glob_param.n1_IR;
n2_IR = glob_param.n2_IR;
nb_IR = glob_param.nb_IR;

n1_green = glob_param.n1_green;
n2_green = glob_param.n2_green;
nb_green = glob_param.nb_green;

n1_ratio = glob_param.n1_ratio;
n2_ratio = glob_param.n2_ratio;

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
    
    n_c = [];
    L_green = [];
    for kk = 1:(no_of_stacks)
        n_c  = [n_c n1_green n2_green];
        L_green(2*kk - 1) = L(2*kk - 1)*n1_ratio;
        L_green(2*kk) = L(2*kk)*n2_ratio;
    end
    n_green = [na n_c nb_green];  % add the indices of the vacuum and the substrate at either end
elseif strcmp(glob_param.coatingType, 'AR') %AR coating
    n_c = [];
    for kk = 1:(no_of_stacks)
        n_c  = [n_c n2_IR n1_IR];
    end
    n_IR = [nb_IR n_c na];  % add the indices of the vacuum and the substrate at either end
    
    n_c = [];
    L_green = [];
    for kk = 1:(no_of_stacks)
        n_c  = [n_c n2_green n1_green];
        L_green(2*kk - 1) = L(2*kk - 1)*n2_ratio;
        L_green(2*kk) = L(2*kk)*n1_ratio;
    end
    n_green = [nb_green n_c na];  % add the indices of the vacuum and the substrate at either end
end

% list of wavelengths to use for the optimization
lambda_IR = [lambda_0] / lambda_0;
lambda_green = [lambda_0/2] / lambda_0;

% Calculate reflectivities for each polarization at the angle of incidence
% ['tm' = 'p-pol', 'te' = 's-pol'], for the given layer structure
[Gammap, ~]   = multidiel100(n_IR, L, lambda_IR, aoi, 'tm');
Rp_IR = abs(Gammap).^2;
Tp_IR = 1 - Rp_IR;

% minimize reflected E-field (usually by 1/2 wave cap on top)
% only at main wavelength
r_refl = abs(1 + Gammap);

%Now for green
[Gammap, ~]   = multidiel100(n_green, L_green, lambda_green, aoi_green, 'tm');
Rp_green = abs(Gammap).^2;
Tp_green = 1 - Rp_green;

% define the error function
T_1 = glob_param.T_1;         % desired  transmission @ lambda_0
R_1 = 1 - T_1;

T_2p = glob_param.T_2p;         % desired transmission @ lambda_0 / 2
R_2p = 1 - T_2p;      % currently not used

surfTarg = glob_param.surfaceField;

% cost function which gets minimized 
yy = [];
% Indices of the yy term
% 1 = Transmission at 1064, p-pol
% 2 = Transmission at 532, p-pol
% 3 = HR surface field
% 4 = Brownian Noise
% 5 = TO noise (Thermo-refractive + Thermoelastic)
% 6 = sensitivity to change in layer thickness @1064nm p- +/-1%...
% 7 = sensitivity to change in layer thickness @532nm p-pol +/-1%...
% 8 = sensitivity to change in n1 @1064nm p- +/-1%...
% 9 = sensitivity to change in n1 @532nm p-pol +/-1%...
% 10 = sensitivity to change in n2 @1064nm p- +/-1%...
% 11 = sensitivity to change in n2 @532nm p-pol +/-1%...
% 12 = sensitivity to change in aoi @1064nm p- +/-1%...
% 13 = sensitivity to change in aoi @532nm p-pol +/-1%...

%Don't penalize for being better than the spec...
yy = [yy heaviside(Tp_IR - T_1)*abs((Tp_IR - T_1)/T_1)^1];    % match the T @ lambda

yy = [yy heaviside(R_2p - Rp_green)*abs((Rp_green - R_2p)/R_2p)^1];  % match the T @ lambda/2, p-pol
%yy = [yy abs((Rp_green - R_2p)/R_2p)^1];  % match the T @ lambda/2, p-pol


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
    y = y + weights(5) * TO;
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
    y.Rp_IR = Rp_IR;
    y.Tp_green = Tp_green;
    y.Rp_green = Rp_green;
    y.Sthermal = S;
    y.TO = TO;
    y.yy = yy;
    y.n_IR = n_IR;
    y.n_green = n_green;
    y.weights = weights;
end

end
