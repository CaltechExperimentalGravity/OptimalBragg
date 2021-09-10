% aSiModel returns a structure describing an IFO for use in
% benchmark programs and noise simulator. Part of the gwinc
% package, which provides science-grounded figures of merit for
% comparing interferometric gravitational wave detector designs. 

% in this case its just used for thermal noise calcs


function ifo = aSiModel(varargin)



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

  %% high index material: a-Si
  %  http://www.mit.edu/~6.777/matprops/asi.htm
  %  https://wiki.ligo.org/OPT/AmorphousSilicon
  ifo.Materials.Coating.Yhighn     = 60e9;      % http://dx.doi.org/10.1063/1.344462
  ifo.Materials.Coating.Sigmahighn = 0.22;

  % http://dx.doi.org/10.1103/PhysRevLett.96.055902
  c_over_T3 = 0.23;
  ifo.Materials.Coating.CVhighn    = c_over_T3 * 123^3 * 0.001;
  ifo.Materials.Coating.Alphahighn = 1e-9;  % guess; offset form
                                            % zero CTE

  ifo.Materials.Coating.Betahighn  = 1.4e-4;           % dn/dT (http://dx.doi.org/10.1063/1.1383056)
  ifo.Materials.Coating.ThermalDiffusivityhighn = 1.03; % W/m/K |http://dx.doi.org/10.1103/PhysRevLett.96.055902
  ifo.Materials.Coating.Phihighn   = 1e-5;             % just a guess (depends on prep)
  ifo.Materials.Coating.Indexhighn = 3.5;


  %% low index material: Silicon-Nitride   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  %  https://wiki.ligo.org/OPT/CoatingSiN
  %
  %  1) doi:10.1016/j.ssc.2003.08.048
  %  2) doi.org/10.1364/AO.51.007229
  %  3) https://www.filmetrics.com/refractive-index-database/Si3N4/Silicon-Nitride-SiN-SiON
  %  4) https://doi.org/10.1063/1.1638886
  %  5) 10.1109/JPHOT.2016.2561622
  %     (dndT= 3.21e-7 - 1.99e-8 * T + 8.856e-10 * T.^2 - 1.375e-12 *...
  %            T.^3 - 1.1e-15 * T.^4;)
  %  6) 10.1103/PhysRevD.98.102001 (PRD 2018)
  %  7) https://doi.org/10.1103/PhysRevD.97.022004 (S. Chao 2018);
  %     Y varies from 103 - 137 GPa depending on dep method. Could
  %     use other methods to go from 85 - 210 GPa


  ifo.Materials.Coating.Ylown      = 140e9;               % ref[6]
  ifo.Materials.Coating.Sigmalown  = 0.25;                % ref[6]
  ifo.Materials.Coating.CVlown     = 1e-7 * 123^3 * 1000; % ref[1]
  ifo.Materials.Coating.Alphalown  = 3.3e-6;              % ref[2]
  ifo.Materials.Coating.Betalown   = 8.5e-6;              % dn/dT ref[5]

  % really thermal conductivity
  ifo.Materials.Coating.ThermalDiffusivitylown = 2;       % ref[1]
  ifo.Materials.Coating.Philown    = 1e-4;                % need to measure with low stress

  ifo.Materials.Coating.Indexlown  = 1.8;                % ref [3]
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  %% Silicon Substrate Material parameters--------------------------------------------
  % Silicon @ 120K (http://www.ioffe.ru/SVA/NSM/Semicond/Si/index.html)
                                                         %  phi_sub = c2 * f^(MechLossExp)
  ifo.Materials.Substrate.c2                = 3e-13;   % Coeff of freq dep. term for bulk loss (Lam & Douglass, 1981)
  ifo.Materials.Substrate.MechanicalLossExponent = 1;   % Exponent for freq dependence of silicon loss
  ifo.Materials.Substrate.Alphas            = 5.2e-12;   % Surface loss limit ???
  ifo.Materials.Substrate.MirrorY           = 155.8e9;   % N/m^2; Youngs modulus (ioffe)
  ifo.Materials.Substrate.MirrorSigma       = 0.27;      % kg/m^3; Poisson ratio (ioffe)
  ifo.Materials.Substrate.MassDensity       = 2329;      % kg/m^3; (ioffe)
  ifo.Materials.Substrate.MassAlpha         = 1e-9;      % 1/K; CTE = 0 @ 120 K
  ifo.Materials.Substrate.MassCM            = 0.3*1000;  % J/kg/K; specific heat (ioffe @ 120K)
  ifo.Materials.Substrate.MassKappa         = 700;       % W/(m*K); thermal conductivity (ioffe @ 120)
  ifo.Materials.Substrate.RefractiveIndex   = 3.5;       % 3.38 * (1 + 4e-5 * T)   (ioffe)



  ifo.Materials.MassRadius    = 0.430/2;             % m
  ifo.Materials.MassThickness = 0.55;

  ifo.Materials.Substrate.Temp = 123;            % mirror temperature [K]
  %% Laser-------------------------------------------------------------------
  ifo.Laser.Wavelength                   = 2128e-9;                                  % m;
  ifo.Laser.Power                        = 125;                                       % W;

  %% Optics------------------------------------------------------------------

  ifo.Optics.Curvature.ITM = 1970;               % ROC of ITM
  ifo.Optics.Curvature.ETM = 2192;               % ROC of ETM
  ifo.Optics.SubstrateAbsorption = 0.5e-4;       % 1/m; bulk absorption coef (ref. 2)


  % 
  ifo.Optics.ITM.BeamRadius = 0.059;                     % m; 1/e^2 power radius
  ifo.Optics.ETM.BeamRadius = 0.084;                     % m; 1/e^2 power radius

  % coating layer optical thicknesses - mevans June 2008
  ifo.Optics.ITM.CoatingThicknessLown = 0.308;
  ifo.Optics.ITM.CoatingThicknessCap  = 0.5;

  ifo.Optics.ETM.CoatingThicknessLown = 0.27;
  ifo.Optics.ETM.CoatingThicknessCap  = 0.5;

  % Defining some additional fields for the structure for pygwinc compatibility
  ifo.Optics.ETM.Coating       = ifo.Materials.Coating;
  ifo.Optics.ETM.Substrate     = ifo.Materials.Substrate;
  ifo.Optics.ETM.MassRadius    = ifo.Materials.MassRadius;
  ifo.Optics.ETM.MassThickness = ifo.Materials.MassThickness;

end

