% runOptim.m -  Code for optimizing dielectric coating design using MATLAB Simulated Annealing
%
% Adapted to optimize aLIGO ETM coating requirements - E0900068-v5
% https://git.ligo.org/40m/Coatings
%
% Usage: [OUT] = runOptim;
%        [OUT] = runOptim(settings);
%
% To change default values pass structure "settings", for example:
% settings.Workers = 40; 
% settings.XO = ...; % Initial Guess
%
% OUT.XO   = Optimized Parameter Values
% OUT.fval = Obtained global minima solution for the cost function 
%
% NM & GV 31 March 2018
function [OUT] = runOptim(settings)
% clc
% clear
% close all
addpath('../');
addpath('thermalNoiseFuncs/');
addpath('./SA_utils/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    DEFAULT SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defaults.Workers = 40;
defaults.parallel_toggle = 1;

if nargin ~= 0
    settings = setstructfields(defaults,settings);
else
    settings = defaults;
end



OUT = [];


%Start a parpool with N  workers
delete(gcp('nocreate'));

if (settings.parallel_toggle == 1)
    parpool('local',settings.Workers)   ;
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
NUMTOOLS.opt_name       = 'aLIGO_ETM';
NUMTOOLS.coatingType    = 'HR';
NUMTOOLS.T_1            = 5e-6;     % Target transmission at 1064 nm --- 5ppm for HR
NUMTOOLS.T_2p           = 0.02;     % Target transmission at 532 nm, p-pol ---- same as above
NUMTOOLS.surfaceField   = 0.01;     % Target Surface E Field, [V/m]
NUMTOOLS.aoi            = 0;       % Angle of incidence (degrees)  --- 41.1deg for HR, 24.8deg for AR
NUMTOOLS.aoi_green      = 0;       % Angle of incidence (degrees)  ---24.746 deg for green
NUMTOOLS.wBeam        = 6e-2;       % Beam size on optic being optimized
NUMTOOLS.f_optimize   = 100;        % Frequency at which to evaluate noise being optimized
NUMTOOLS.noise_weight = 1e21;       % to equate Brownian and TO noise
NUMTOOLS.include_sens   = 1;        % Include derivatives in sensitivity function (1 or 0)
NUMTOOLS.include_brownian  = 1;     % Include coating brownian term in sensitivity function
NUMTOOLS.include_TO   = 1;          % Include thermo-optic term in sensitivity function

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
% 2 = Transmission at 532, p-pol
% 3 = HR surface field
% 4 = Brownian Noise
% 5 = TO Noise
% 6 = sensitivity to change in layer thickness @1064nm p- +/-1%...
% 7 = sensitivity to change in layer thickness @532nm p-pol +/-1%...
% 8 = sensitivity to change in n1 @1064nm p- +/-1%...
% 9 = sensitivity to change in n1 @532nm p-pol +/-1%...
% 10 = sensitivity to change in n2 @1064nm p- +/-1%...
% 11 = sensitivity to change in n2 @532nm p-pol +/-1%...
% 12 = sensitivity of surface field term to perturbations in n1, n2 and L


%weights = [10000 1000 1 10 1e22 5 1 5 1 5 1 0.01];

weights = [11111. 11111. 1. 10 1e22 50. 10. 50. 10. 50. 10. 0.01];

NUMTOOLS.wt = weights;

% setting the bounds for the variables to be searched over
LB = 0.030 * ones(size(x0));     % lower bounds on the layer thicknesses, lambda/50
UB = 0.510 * ones(size(x0));     % upper bound, lambda/2
nvars = length(UB);


XO = 0.5*(LB+UB);

if isfield(settings,'XO')
    XO = settings.XO;
end


Temp_set = logspace(4,-2,3);
RI_set   = floor(linspace(5,30,2));

save('RI_set','RI_set');
save('Temp_set','Temp_set');


if strcmp(NUMTOOLS.coatingType,'HR')
    %    LB(1) = 0.45;  %Constrain topost layer to be nearly halfwave thick...
    %    UB(1) = 0.55;
    NUMTOOLS.firstLayer = 'SiO2';     %Depending whether coating is on the HR side or AR side
else
    NUMTOOLS.firstLayer = 'Ta2O5';
end

ObjectiveFunction = @(x) getCost_aLIGO_ETM_mex([x], NUMTOOLS, 0);


save('X_initial','XO');

Residual = ObjectiveFunction(XO);

tic;
diary('SA_Diary.txt'); % text file created to save the results
disp(['%%%======== Parallelized Cascaded  Simulated Annealing =========%%%'])
disp(['%%%=============                                  ==============%%%'])
disp(['%%%=============','       ', datestr(clock),'       ','=============%%%'])
disp(['%%%=============                                  ==============%%%'])
disp(['%%%=============================================================%%%'])


disp(['Initial_Residual = ',num2str(Residual)])
disp( ['Initial  Values = ']);
disp([XO' ]);

parsave(sprintf('Optim.mat'), XO,Residual);


% set particle swarm optimization options
hybridopts = optimoptions('fmincon',...
    'Display','off',...
    'MaxIter', 2500,...
    'TolFun', 1e-8,...
    'MaxFunEvals', 2500);


for Init_temp_index = 1:length(Temp_set)
    for RI_index = 1:length(RI_set)
        
        [Xinp] = parload(sprintf('Optim.mat'));
        disp(['Current Best optimum (input to workers) = ',num2str(ObjectiveFunction(Xinp) )])
        
        parfor parallel_worker = 1:settings.Workers
            
            % Prepare input for workers (slightly modify the previous optimal parameters)
            XO = Xinp + Xinp.*randn(size(Xinp))*0.001;
            
            RI = SA_RI_load(RI_index);
            Temp = SA_Temp_load(Init_temp_index);
            
            SA_options = saoptimset('InitialTemperature',Temp,'ReannealInterval',RI,'Display','final','AnnealingFcn',...
                @annealingfast,'TemperatureFcn', @temperatureexp,'AcceptanceFcn',...
                @acceptancesa,'TolFun', 1.0000e-006 ,'StallIterLimit',...
                500,'MaxFunEvals',1000,'HybridFcn',{@fmincon,hybridopts}  );%'HybridFcn',{@patternsearch,hybridopts}
            
            
            [xout,fval,exitFlag,output] = simulannealbnd(ObjectiveFunction,XO,LB,UB,SA_options);
            
            Res = ObjectiveFunction(xout);
            
            parsave(sprintf('output%d.mat', parallel_worker), xout,Res);
            
        end  % parfor loop end
        
        
        % Update optimal values
        
        for worker_index = 1:settings.Workers
            [Xout,Res] = parload(sprintf('output%d.mat', worker_index));
            disp([' Residual for Worker ',num2str(worker_index),' =  ',num2str(Res)])
            % disp( ['Worker ',num2str(worker_index),' optimized parameters = ']);
            % disp([Xout' ]);
            
            if (Res < Residual)
                disp(['better solution than previous best(',num2str(Residual) ,')found by worker',num2str(worker_index),' with value ', num2str(Res)])
                Xopt = Xout;
                Residual = Res;
                parsave(sprintf('Optim.mat'), Xopt,  Residual);
            end  % Residual update / if ends
            
        end  % Optim Coeff set / for loop ends
        
        disp(['Initial Temperature =  ',num2str(SA_Temp_load(Init_temp_index))]);
        disp(['Reanneal Interval =   ',num2str(SA_RI_load(RI_index))]);
        disp(['Residual = ', num2str(ObjectiveFunction(Xopt) )]);
        disp( ['Optimized Parameters = ']);
        disp([Xopt' ]);
        toc
    end  % Hierarchical SA - RI-for loop ends
    
end % Hierarchical SA - Init_temp-for loop ends

OUT.executionTime = toc;

if (settings.parallel_toggle == 1)
    delete(gcp('nocreate'))
end
system('rm -r output*.mat RI_set.mat Temp_set.mat X_initial.mat'); % deletes individual  output files
diary off

OUT.fval = Residual;
OUT.xout = Xopt;

% Check for thin layers
%if find(xout < 0.05)
%    disp('Bad Layer Thickness: invalid results')
%end

%% SAVE layer structure, IFOmodel, and noise
costOut = Residual;
tnowstr = datestr(now, 'yymmdd_HHMM');
savename = ['Data/' NUMTOOLS.opt_name '_' NUMTOOLS.coatingType '_' num2str(no_of_stacks) '_layers_' tnowstr];
save(savename, 'costOut', 'ifo');
OUT.tnowstr = tnowstr;
end
