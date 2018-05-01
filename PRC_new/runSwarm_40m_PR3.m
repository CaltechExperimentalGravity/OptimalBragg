% Code for optimizing dielectric coating design using MATLAB particle swarm
% Specs for the 40m PR3/SR3 replacement project.
% https://git.ligo.org/40m/Coatings
function OUT = runSwarm_40m_PR3(settings)

    addpath('../');
    addpath('../generic/');
    addpath('../generic/thermalNoiseFuncs/');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                    DEFAULT SETTINGS                                             
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                              
    defaults.Workers          = 20;        
    defaults.SwarmSize        = 55;      
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

    NUMTOOLS.n1_green = interp1(SiO2(:,1),SiO2(:,2),532,'pchip');
    NUMTOOLS.n2_green = interp1(Ta2O5(:,1),Ta2O5(:,2),532,'pchip');
    NUMTOOLS.nb_green = 1.4607; %Estimate

    NUMTOOLS.n1_ratio = NUMTOOLS.n1_green / NUMTOOLS.n1_IR;
    NUMTOOLS.n2_ratio = NUMTOOLS.n2_green / NUMTOOLS.n2_IR;

    % Initial guess vector of layer thicknesses in units of lambda
    no_of_stacks = 20;
    x0 = [];
    for kk = 1:no_of_stacks
        x0 = [x0 2/8 2/8];
    end

    clear costOut
    clear ifo

    NUMTOOLS.lambda         = 1064e-9;  % Wavelength at which R and T are to be evaluated
    NUMTOOLS.opt_name       = 'PR3';
    NUMTOOLS.coatingType    = 'HR';
    NUMTOOLS.T_1            = 50e-6;     % Target transmission at 1064 nm --- 50ppm for HR
    NUMTOOLS.aoi            = 41.1;       % Angle of incidence (degrees)  --- 41.1deg for HR, 24.8deg for AR
    %NUMTOOLS.surfaceField   = 0.01;     % Target surface electric field, [V/m]
    NUMTOOLS.include_sens   = 1;        % Include derivatives in sensitivity function (1 or 0)
    NUMTOOLS.include_brownian  = 0;     % Include coating brownian term in sensitivity function
    NUMTOOLS.include_TO   = 0;          % Include thermo-optic term in sensitivity function

    % Load a standard gwinc style parameter file with various material properties defined
    ifo = SilicaTantala300;
    ifo.Laser.Wavelength = NUMTOOLS.lambda;
    % load this lookup table for speedup
    load besselzeros.mat
    ifo.Constants.BesselZeros = besselzeros;
    %set some of the material params in this structure to match that used to optimize coating, i.e. values from Ramin
    ifo.Materials.Coating.Indexhighn = NUMTOOLS.n2_IR;
    ifo.Materials.Coating.Indexlown  = NUMTOOLS.n1_IR;
    ifo.Materials.Substrate.RefractiveIndex = 1.449641; % Corning datasheet
    NUMTOOLS.ifo = ifo;

    %set up the weighting to convert MO cost to scalar cost
    % Indices of the yy term
    % 1 = Transmission at 1064, p-pol
    % 2 = sensitivity to change in layer thickness @1064nm p- +/-1%...
    % 3 = sensitivity to change in n1 @1064nm p- +/-1%...
    % 4 = sensitivity to change in n2 @1064nm p- +/-1%...
    % 5 = sensitivity to change in aoi +/-5%...


    %weights = [10000 1000 1 10 1e22 5 1 5 1 5 1 0.01];

    weights = [100/60 10/40 10/30 10/40 10/30];

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
        particleswarm(@(x) getCost_40m_PR3([x], NUMTOOLS, 0),...
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
    costOut = getCost_40m_PR3(xout, NUMTOOLS, 1);
    tnowstr = datestr(now, 'yymmdd_HHMM');
    savename = ['Data/' NUMTOOLS.opt_name '_' NUMTOOLS.coatingType '_' num2str(no_of_stacks) '_layers_' tnowstr];
    save(savename, 'costOut', 'ifo');
    OUT.tnowstr = tnowstr;
end
