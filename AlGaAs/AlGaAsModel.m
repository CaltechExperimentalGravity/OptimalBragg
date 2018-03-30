% IFOMODEL returns a structure describing an IFO for use in
% benchmark programs and noise simulator. Part of the gwinc
% package, which provides science-grounded figures of merit for
% comparing interferometric gravitational wave detector designs. 
% 


% parameters for quad pendulum suspension updated 3rd May 2006, NAR
% References:
% LIGO-T000012-00-D
% 	* Differentiate between silica and sapphire substrate absorption
% 	* Change ribbon suspension aspect ratio
% 	* Change pendulum frequency
% * References:
% * 1. Electro-Optic Handbook, Waynant & Ediger (McGraw-Hill: 1993)
% * 2. LIGO/GEO data/experience
% * 3. Suspension reference design, LIGO-T000012-00
% * 4. Quartz Glass for Optics Data and Properties, Heraeus data sheet,
% *    numbers for suprasil
% * 5. Y.S. Touloukian (ed), Thermophysical Properties of Matter 
% *    (IFI/Plenum,1970)
% * 6. Marvin J. Weber (ed) CRC Handbook of laser science and technology, 
% *    Vol 4, Pt 2
% * 7. R.S. Krishnan et al.,Thermal Expansion of Crystals, Pergamon Press
% * 8. P. Klocek, Handbook of infrared and optical materials, Marcel Decker, 
% *    1991
% * 9. Rai Weiss, electronic log from 5/10/2006
% * 10. Wikipedia online encyclopedia, 2006
% * 11. D.K. Davies, The Generation and Dissipation of Static Charge on
% * dielectrics in a Vacuum, page 29
% * 12. Gretarsson & Harry, Gretarsson thesis
% * 13. Fejer
% * 14. Braginsky

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
  
  ifo.Materials.Substrate.c2  = 7.6e-12;                 % Coeff of freq depend. term for bulk mechanical loss, 7.15e-12 for Sup2
  ifo.Materials.Substrate.MechanicalLossExponent=0.77;   % Exponent for freq dependence of silica loss, 0.822 for Sup2
  ifo.Materials.Substrate.Alphas = 5.2e-12;              % Surface loss limit (ref. 12)
  ifo.Materials.Substrate.MirrorY    = 72.7e9;           % N/m^2; Youngs modulus (ref. 4)
  ifo.Materials.Substrate.MirrorSigma = 0.167;           % Kg/m^3; Poisson ratio (ref. 4)
  ifo.Materials.Substrate.MassDensity = 2200;            % Kg/m^3; (ref. 4)
  ifo.Materials.Substrate.MassAlpha = 3.9e-7;            % 1/K; thermal expansion coeff. (ref. 4)
  ifo.Materials.Substrate.MassCM = 739;                  % J/Kg/K; specific heat (ref. 4)
  ifo.Materials.Substrate.MassKappa = 1.38;              % J/m/s/K; thermal conductivity (ref. 4)
  ifo.Materials.Substrate.RefractiveIndex = 1.45;        % mevans 25 Apr 2008
  
  ifo.Materials.MassRadius    = (34e-2)/2;             % m  (1" dia)
  ifo.Materials.MassThickness = 20e-2;             % m  (1/4" thick)
  
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
  
  % factor of 2.5 added to simulate LNG modes - remove after new LNG code is added
  ifo.Optics.ITM.BeamRadius = 0.055 * 2.5;                     % m; 1/e^2 power radius
  ifo.Optics.ETM.BeamRadius = 0.062 * 2.5;                     % m; 1/e^2 power radius
  
  ifo.Optics.ITM.CoatingAbsorption = 0.5e-6;            % absorption of ITM
  ifo.Optics.ITM.Transmittance  = 0.014;                % Transmittance of ITM
  ifo.Optics.ETM.Transmittance  = 5e-6;                 % Transmittance of ETM
  ifo.Optics.SRM.Transmittance  = 0.20;                 % Transmittance of SRM
  ifo.Optics.PRM.Transmittance  = 0.03;
  
  % coating layer optical thicknesses - mevans June 2008
  ifo.Optics.ITM.CoatingThicknessLown = 0.308;
  ifo.Optics.ITM.CoatingThicknessCap = 0.5;
  
  ifo.Optics.ETM.CoatingThicknessLown = 0.27;
  ifo.Optics.ETM.CoatingThicknessCap = 0.5;
  
  %ifo.Optics.SRM.Tunephase = 0.23;           % SRM tuning, 795 Hz narrowband
  ifo.Optics.SRM.Tunephase = 0.0;             % SRM tuning
  ifo.Optics.Quadrature.dc = pi/2;            % demod/detection/homodyne phase
  
  
  
end
