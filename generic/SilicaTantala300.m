% aSiModel returns a structure describing an IFO for use in
% benchmark programs and noise simulator. Part of the gwinc
% package, which provides science-grounded figures of merit for
% comparing interferometric gravitational wave detector designs. 

% in this case its just used for thermal noise calcs


function ifo = SilicaTantala300(varargin)


  
  %% Physical and other constantMaterialss; All in SI units------------------
  
  ifo.Constants.E0   = 8.8541878176e-12;            % F/m; Permittivity of Free Space
  ifo.Constants.hbar = 1.054572e-34;                % J-s; (Plancks constant)/(2*pi)
  ifo.Constants.c    = 299792458;                   % m/s; Speed of light in Vacuum
  ifo.Constants.G    = 6.67259e-11;                 % m^3/Kg/s^2; Grav. Constant
  ifo.Constants.kB   = 1.380658e-23;                % J/K; Boltzman Constant
  ifo.Constants.h    = ifo.Constants.hbar * 2*pi;   % J-s; Planks constant
  ifo.Constants.R    = 8.31447215;                  % J/(K*mol); Gas Constant
  ifo.Constants.Temp = 290;                         % K; Temperature of the Vacuum
  ifo.Constants.yr   = 365.2422 * 86400;            % sec; Seconds in a year

  ifo.Constants.g = 9.81;                         % m/s^2; grav. acceleration

%% Dielectric coating material parameters----------------------------------

%% high index material: tantala
ifo.Materials.Coating.Yhighn = 140e9;
ifo.Materials.Coating.Sigmahighn = 0.23;
ifo.Materials.Coating.CVhighn = 2.1e6;               % Crooks et al, Fejer et al
ifo.Materials.Coating.Alphahighn = 3.6e-6;           % 3.6e-6 Fejer et al, 5e-6 from Braginsky
ifo.Materials.Coating.Betahighn = 1.4e-5;              % dn/dT, value Gretarrson (G070161)
ifo.Materials.Coating.ThermalDiffusivityhighn = 33;  % Fejer et al
ifo.Materials.Coating.Phihighn = 2.3e-4;
ifo.Materials.Coating.Indexhighn = 2.06539;

%% low index material: silica
ifo.Materials.Coating.Ylown = 72e9;
ifo.Materials.Coating.Sigmalown = 0.17;
ifo.Materials.Coating.CVlown = 1.6412e6;             % Crooks et al, Fejer et al
ifo.Materials.Coating.Alphalown = 5.1e-7;            % Fejer et al
ifo.Materials.Coating.Betalown = 8e-6;             % dn/dT,  (ref. 14)
ifo.Materials.Coating.ThermalDiffusivitylown = 1.38; % Fejer et al
ifo.Materials.Coating.Philown = 4.0e-5;
ifo.Materials.Coating.Indexlown = 1.45;

%%Substrate Material parameters--------------------------------------------

ifo.Materials.Substrate.c2  = 7.6e-12;                 % Coeff of freq depend. term for bulk mechanical loss, 7.15e-12 for Sup2
ifo.Materials.Substrate.MechanicalLossExponent=0.77;   % Exponent for freq dependence of silica loss, 0.822 for Sup2
ifo.Materials.Substrate.Alphas = 5.2e-12;              % Surface loss limit (ref. 12)
ifo.Materials.Substrate.MirrorY    = 7.27e10;          % N/m^2; Youngs modulus (ref. 4)
ifo.Materials.Substrate.MirrorSigma = 0.167;           % Kg/m^3; Poisson ratio (ref. 4)
ifo.Materials.Substrate.MassDensity = 2.2e3;           % Kg/m^3; (ref. 4)
ifo.Materials.Substrate.MassAlpha = 3.9e-7;            % 1/K; thermal expansion coeff. (ref. 4)
ifo.Materials.Substrate.MassCM = 739;                  % J/Kg/K; specific heat (ref. 4)
ifo.Materials.Substrate.MassKappa = 1.38;              % J/m/s/K; thermal conductivity (ref. 4)
ifo.Materials.Substrate.RefractiveIndex = 1.45;        % mevans 25 Apr 2008
  
  %% Silicon Substrate Material parameters--------------------------------------------
  % Silicon @ 120K (http://www.ioffe.ru/SVA/NSM/Semicond/Si/index.html)
                                                         %  phi_sub = c2 * f^(MechLossExp)
%   ifo.Materials.Substrate.c2                = 3e-13;   % Coeff of freq dep. term for bulk loss (Lam & Douglass, 1981)
%   ifo.Materials.Substrate.MechanicalLossExponent = 1;   % Exponent for freq dependence of silicon loss
%   ifo.Materials.Substrate.Alphas            = 5.2e-12;   % Surface loss limit ???
%   ifo.Materials.Substrate.MirrorY           = 155.8e9;   % N/m^2; Youngs modulus (ioffe)
%   ifo.Materials.Substrate.MirrorSigma       = 0.27;      % kg/m^3; Poisson ratio (ioffe)
%   ifo.Materials.Substrate.MassDensity       = 2329;      % kg/m^3; (ioffe)
%   ifo.Materials.Substrate.MassAlpha         = 1e-9;      % 1/K; CTE = 0 @ 120 K
%   ifo.Materials.Substrate.MassCM            = 0.3*1000;  % J/kg/K; specific heat (ioffe @ 120K)
%   ifo.Materials.Substrate.MassKappa         = 700;       % W/(m*K); thermal conductivity (ioffe @ 120)
%   ifo.Materials.Substrate.RefractiveIndex   = 3.5;       % 3.38 * (1 + 4e-5 * T)   (ioffe)



  ifo.Materials.MassRadius    = 0.0254 * 1/2;             % m  
  ifo.Materials.MassThickness = 0.0254 * 1/4;
  
  ifo.Materials.Substrate.Temp = 290;            % mirror temperature [K]
  %% Laser-------------------------------------------------------------------
  ifo.Laser.Wavelength                   = 1064e-9;                                  % m;
  ifo.Laser.Power                        = 125;                                       % W;
  
  %% Optics------------------------------------------------------------------

  ifo.Optics.Curvature.ITM = 1970;               % ROC of ITM
  ifo.Optics.Curvature.ETM = 2192;               % ROC of ETM
  ifo.Optics.SubstrateAbsorption = 0.5e-4;       % 1/m; bulk absorption coef (ref. 2)

  
  % factor of 2.5 added to simulate LNG modes - remove after new LNG code is added
  ifo.Optics.ITM.BeamRadius = 0.055 * 2.5;                     % m; 1/e^2 power radius
  ifo.Optics.ETM.BeamRadius = 0.062 * 2.5;                     % m; 1/e^2 power radius
  
  % coating layer optical thicknesses - mevans June 2008
  ifo.Optics.ITM.CoatingThicknessLown = 0.308;
  ifo.Optics.ITM.CoatingThicknessCap = 0.5;
  
  ifo.Optics.ETM.CoatingThicknessLown = 0.27;
  ifo.Optics.ETM.CoatingThicknessCap = 0.5;
  
  
  
end
