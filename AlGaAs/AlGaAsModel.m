% this is a bunch of mirror params to be used for mirror thermal nois calcs

function ifo = AlGaAsModel(varargin)

  %% Infrastructure----------------------------------------------------------
  
  ifo.Infrastructure.Length                     = 0.05;      % m;
  ifo.Infrastructure.ResidualGas.pressure       = 4.0e-7;    % Pa;
  ifo.Infrastructure.ResidualGas.mass           = 3.35e-27;  % kg; Mass of H_2 (ref. 10)
  ifo.Infrastructure.ResidualGas.polarizability = 7.8e-31 ;  % m^3; Gas polarizability
  
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
  ifo.Constants.Mpc  = ifo.Constants.yr * ifo.Constants.c * 3.26e6;    % m;
  
  ifo.Constants.MSol = 1.989e30*ifo.Constants.G/ifo.Constants.c^2;% m;
  ifo.Constants.g = 9.81;                         % m/s^2; grav. acceleration
  ifo.Constants.fs = 16384;                       % Sampling frequency (Hz)
  
  ifo.Constants.fInspiralMin = 3;  % cut-off for inspiral range (Hz, see int73)
  

  %% AlGaAs coating material parameters----------------------------------
  %
  % from http://www.ioffe.ru/SVA/NSM/Semicond/AlGaAs/  for T = 300 K
  % kappa = 55-212*x+248*x^2   (W / (m K))
  % sigma = 0.31 + 0.1*x
  % Y = (85.3 -1.8*x)*1e9    (Pa)
  % Density = 5320 - 1560*x    (kg / m^3)
  % Heat Capacity = (330 + 120*x)
  % CV = (Heat Capacity) * Density
  % alpha = (5.73-0.53*x)*1e-6     (dL/L)/K
  
  % high index material: GaAs
  ifo.Materials.Coating.Yhighn = 85.3e9;
  ifo.Materials.Coating.Sigmahighn = 0.31;            % Poisson's ratio
  ifo.Materials.Coating.CV = 330;                      % Heat Cap per kg
  ifo.Materials.Coating.Density = 5320;
  ifo.Materials.Coating.CVhighn =...
    ifo.Materials.Coating.CV * ifo.Materials.Coating.Density; % Heat Capacity per Volume
  ifo.Materials.Coating.Alphahighn = 5.73e-6;          % 5.73 ppm/C
  ifo.Materials.Coating.Betahighn = 2.67e-4;           % Appl. Phys. Lett. 66, 335 (1995) 
  % name is wrong; should be Thermal Conductivity
  ifo.Materials.Coating.ThermalDiffusivityhighn = 55;  % 
  ifo.Materials.Coating.Phihighn = 2e-5;             % guess
  ifo.Materials.Coating.Indexhighn = 3.48041;            % G.D.Cole
  
  % low index material: Al_{0.92}Ga_{0.08}As
  ifo.Materials.Coating.Ylown = 83.6e9;
  ifo.Materials.Coating.Sigmalown = 0.32;
  ifo.Materials.Coating.CV = 440;                      % Heat Cap per kg
  ifo.Materials.Coating.Density = 3885;
  ifo.Materials.Coating.CVlown =...
    ifo.Materials.Coating.CV * ifo.Materials.Coating.Density; % Heat Capacity per Volume
  ifo.Materials.Coating.Alphalown = 5.24e-6;             % 
  ifo.Materials.Coating.Betalown = 1.53e-4;             % dn/dT
  ifo.Materials.Coating.ThermalDiffusivitylown = 70;  % 
  ifo.Materials.Coating.Philown = 2e-5;               % guess
  ifo.Materials.Coating.Indexlown = 2.97717;             % 
  
  
  %%   Substrate Material parameters for fused silica --------------------
  
  ifo.Materials.Substrate.c2  = 7.6e-12;                 % Coeff of freq depend. term for bulk mechanical loss,
                                                         %  7.15e-12 for Sup2
  ifo.Materials.Substrate.MechanicalLossExponent=0.77;   % Exponent for freq dependence of silica loss, 0.822 for Sup2
  ifo.Materials.Substrate.Alphas = 5.2e-12;              % Surface loss limit (ref. 12)
  ifo.Materials.Substrate.MirrorY    = 72.7e9;           % N/m^2; Youngs modulus (ref. 4)
  ifo.Materials.Substrate.MirrorSigma = 0.167;           % Kg/m^3; Poisson ratio (ref. 4)
  ifo.Materials.Substrate.MassDensity = 2200;            % Kg/m^3; (ref. 4)
  ifo.Materials.Substrate.MassAlpha = 3.9e-7;            % 1/K; thermal expansion coeff. (ref. 4)
  ifo.Materials.Substrate.MassCM = 739;                  % J/Kg/K; specific heat (ref. 4)
  ifo.Materials.Substrate.MassKappa = 1.38;              % J/m/s/K; thermal conductivity (ref. 4)
  ifo.Materials.Substrate.RefractiveIndex = 1.45;        % mevans 25 Apr 2008
  
  ifo.Materials.MassRadius    = 0.0254 / 2;              % m  (1" dia)
  ifo.Materials.MassThickness = 0.0254 / 4;              % m  (1/4" thick)
  
  %% Laser-------------------------------------------------------------------
  ifo.Laser.Wavelength                   = 1.064e-6;                                  % m;
  ifo.Laser.Power                        = 1;                                       % W;
  
  %% Optics------------------------------------------------------------------
  ifo.Optics.Type = 'SignalRecycled';
  
  ifo.Optics.SRM.CavityLength         = 55;      % m; ITM to SRM distance
  ifo.Optics.PhotoDetectorEfficiency  = 0.95;     % photo-detector quantum efficiency
  ifo.Optics.Loss                     = 37.5e-6; % average per mirror power loss
  ifo.Optics.BSLoss  = 0.5e-3;                   % power loss near beamsplitter
  ifo.Optics.coupling = 1.0;                   % mismatch btwn arms & SRC modes; used to
  % calculate an effective r_srm
  ifo.Optics.Curvature.ITM = 1970;               % ROC of ITM
  ifo.Optics.Curvature.ETM = 2192;               % ROC of ETM
  ifo.Optics.SubstrateAbsorption = 0.5e-4;       % 1/m; bulk absorption coef (ref. 2)
  ifo.Optics.pcrit = 10;                         % W; tolerable heating power (factor 1 ATC)
  
  
  ifo.Optics.SRM.Tunephase = 0.0;             % SRM tuning
  ifo.Optics.Quadrature.dc = pi/2;            % demod/detection/homodyne phase
  
  
  
end
