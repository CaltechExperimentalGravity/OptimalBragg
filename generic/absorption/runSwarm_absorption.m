% Code for optimizing dielectric coating design using MATLAB particle swarm
% This wrapper runs the swarm to try and minimize absorption in a hypothetical 
% bilayer consisting of Silica, and X, where X has 10x the absorption of silica
% at 1064nm. Target is to have T=5ppm.
% https://git.ligo.org/40m/Coatings
function OUT = runSwarm_aLIGO_ETM(settings)

    addpath('../');
    addpath('../thermalNoiseFuncs/');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                    DEFAULT SETTINGS                                             
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                              
    %defaults.Workers          = 4; % For laptop        
    defaults.Workers          = 20;        
    defaults.SwarmSize        = 155;      
    defaults.MaxIter          = 1111;      
    defaults.SelfAdjustment   = 1.49;
    defaults.SocialAdjustment = 1.49;
    defaults.TolFun           =  1e-4;
                                                                                    
                                                                                    
    if nargin ~= 0                                                                  
            settings = setstructfields(defaults,settings);                          
    else                                                                            
            settings = defaults;                                                    
    end                                                                             
                                                                                    
                                                                                                                                                                                                                                                        
    %Start a parpool with N  workers                                                
    delete(gcp('nocreate'));  
    
    if settings.Workers > 0
    parpool(settings.Workers);  
    parallel_toggle = 1;
    pctRunOnAll warning off
    warning('off');
    else
        parallel_toggle = 0;
        warning('off');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                    USER CONFIG
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%Load dispersion data...
    load dispersion_revised.mat;
    na = 1.0000;           %Vacuum
    NUMTOOLS.n1_IR = interp1(SiO2(:,1),SiO2(:,2),1064,'pchip');
    NUMTOOLS.n2_IR = interp1(Ta2O5(:,1),Ta2O5(:,2),1064,'pchip');
    NUMTOOLS.nb_IR = 1.449641;   %From Ramin
    NUMTOOLS.aoi = 0;            %Angle of incidence
    NUMTOOLS.alpha1 = 1e-3;      %m^-1, as per https://journals.aps.org/prd/pdf/10.1103/PhysRevD.91.042002
    NUMTOOLS.alpha2 = 1e-2;      %m^-1 (hypothetically x10 of the SiO2 layers)
    NUMTOOLS.nPts = 5;           % Number of points at which to evaluate field in each coating layer, for absorption calc
    % Initial guess vector of layer thicknesses in units of lambda
    no_of_stacks = 20;
    x0 = [];
    for kk = 1:no_of_stacks
        x0 = [x0 2/8 2/8];
    end

    clear costOut
    clear ifo

    NUMTOOLS.lambda         = 1064e-9;  % Wavelength at which R and T are to be evaluated
    NUMTOOLS.opt_name       = 'absorp_ETM';
    NUMTOOLS.coatingType    = 'HR';
    NUMTOOLS.T_1            = 5e-6;     % Target transmission at 1064 nm --- 5ppm for HR
    NUMTOOLS.surfaceField   = 0.01;     % Target Surface E Field, [V/m]
    NUMTOOLS.wBeam        = 6.5e-2;       % Beam size on optic being optimized
    NUMTOOLS.f_optimize   = 100;        % Frequency at which to evaluate noise being optimized
    NUMTOOLS.noise_weight = 1e21;       % to equate Brownian and TO noise
    NUMTOOLS.include_sens   = 0;        % Include derivatives in sensitivity function (1 or 0)
    NUMTOOLS.include_brownian  = 0;     % Include coating brownian term in sensitivity function
    NUMTOOLS.include_TO   = 0;          % Include thermo-optic term in sensitivity function
    NUMTOOLS.minAbsorp   = 1;          % Include thermo-optic term in sensitivity function

    % Load a standard gwinc style parameter file with various material properties defined
    ifo = SilicaTantala300;
    ifo.Laser.Wavelength = NUMTOOLS.lambda;
    ifo.Materials.Coating.Phihighn = 2.3e-4;
    %ifo.Materials.Coating.Phihighn = loss_highn;
    % load this lookup table for speedup
    load besselzeros.mat
    ifo.Constants.BesselZeros = besselzeros;
    ifo.Optics.ETM.BeamRadius = NUMTOOLS.wBeam; %using the gwinc infrastructure for calculating TO noise for ease...
    ifo.Optics.ITM.BeamRadius = NUMTOOLS.wBeam;
    %set some of the material params in this structure to match that used to optimize coating, i.e. values from Ramin
    ifo.Materials.Coating.Indexhighn = NUMTOOLS.n2_IR;
    ifo.Materials.Coating.Indexlown  = NUMTOOLS.n1_IR;
    ifo.Materials.Substrate.RefractiveIndex = 1.449641; % Corning datasheet
    NUMTOOLS.ifo = ifo;

    %set up the weighting to convert MO cost to scalar cost
    % Indices of the yy term
    % 1 = Transmission at 1064, p-pol
    % 2 = Transmission at 1064, p-pol
    % 3 = HR surface field
    
    %weights = [11 0 0]; % For optimizing T only
    weights = [11 11 0]; % For optimizing T and absorption
    %weights = [11 0 11]; % For optimizing T and surface field

    NUMTOOLS.wt = weights;

    % setting the bounds for the variables to be searched over
    LB = 0.030 * ones(size(x0));     % lower bounds on the layer thicknesses, lambda/50
    UB = 0.510 * ones(size(x0));     % upper bound, lambda/2
    nvars = length(UB);

    if strcmp(NUMTOOLS.coatingType,'HR')
    %    LB(1) = 0.45;  %Constrain topost layer to be nearly halfwave thick...
    %    UB(1) = 0.55;
        NUMTOOLS.firstLayer = 'SiO2';     %Depending whether coating is on the HR side or AR side
    else
        NUMTOOLS.firstLayer = 'Ta2O5';
    end

    % set particle swarm optimization options
    hybridopts = optimoptions('fmincon',...
                    'Display','iter',...
                    'MaxIter', 1000,...
                    'TolFun', 1e-4,...
                    'MaxFunEvals', 5000);

    options = optimoptions('particleswarm',...
                'SwarmSize', nvars*settings.SwarmSize,...
                'UseParallel', parallel_toggle,...
                'MaxIter', settings.MaxIter,...
                'SelfAdjustment',   settings.SelfAdjustment,...
                'SocialAdjustment', settings.SocialAdjustment,...
                'TolFun', settings.TolFun,...
                'Display', 'iter',...
                'MaxTime',3600,...
                'HybridFcn',{@fmincon, hybridopts});

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   END USER CONFIG
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   Run particle swarm optimization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x0 = x0(:);
    %% Do the optimization
    % RUNS the Particle Swarm ========
    tic
    [xout, fval, exitflag] =...
        particleswarm(@(x) getCost_absorption([x], NUMTOOLS, 0),...
                                        nvars, LB, UB, options);
    OUT.executionTime = toc                                                         
    OUT.xout = xout;                                                                
    OUT.fval = fval;                                                                                                                                                                                
    OUT.exitflag = exitflag;

    % Check for thin layers
    if find(xout < 0.05)
        disp('Bad Layer Thickness: invalid results')
    end
    %% SAVE layer structure, IFOmodel, and noise
    %costOut = getCost_aLIGO_ETM(xout, NUMTOOLS, 1);
    costOut = getCost_absorption(xout, NUMTOOLS, 1);
    tnowstr = datestr(now, 'yymmdd_HHMM');
    savename = ['Data/' NUMTOOLS.opt_name '_' NUMTOOLS.coatingType '_' num2str(no_of_stacks) '_layers_' tnowstr];
    save(savename, 'costOut', 'ifo');
    OUT.tnowstr = tnowstr;
end
