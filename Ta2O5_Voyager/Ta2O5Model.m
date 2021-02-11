% Ta2O5Model returns a structure describing an IFO for use in
% benchmark programs and noise simulator. Part of the gwinc
% package, which provides science-grounded figures of merit for
% comparing interferometric gravitational wave detector designs. 

% in this case its just used for thermal noise calcs


function ifo = Ta2O5Model(varargin)



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
  % https://git.ligo.org/voyager/mariner40/-/wikis/optics/dichroic-coatings

  %% high index material: Ta2O5
  % https://doi.org/10.1364/AO.48.004536
  ifo.Materials.Coating.Yhighn     = 136e9;     
  ifo.Materials.Coating.Sigmahighn = 0.22;
  ifo.Materials.Coating.CVhighn    = 1.355e6;           % Specific heat (Cp * rho ~ 165.68 * 8180)
  ifo.Materials.Coating.Alphahighn = 0.09e-6;           % Lowest value
  ifo.Materials.Coating.Betahighn  = 0.4e-6;            % dn/dT (http://dx.doi.org/10.1063/1.1383056)
  ifo.Materials.Coating.ThermalDiffusivityhighn = 1.03; % As good as anyones guess
  ifo.Materials.Coating.Phihighn   = 5e-4;              % https://arxiv.org/pdf/1903.06094.pdf
  ifo.Materials.Coating.Indexhighn = 2.083;

  %% low index material: silica
  %  https://wiki.ligo.org/OPT/SilicaCoatingProp
  ifo.Materials.Coating.Ylown      = 72e9;  
  ifo.Materials.Coating.Sigmalown  = 0.17;
  ifo.Materials.Coating.CVlown     = 0.744e6;           
  ifo.Materials.Coating.Alphalown  = 0.0145e-6;            
  ifo.Materials.Coating.Betalown   = 4.2e-6;            
  ifo.Materials.Coating.ThermalDiffusivitylown = 1.05;  % http://dx.doi.org/10.1109/ITHERM.2002.1012450
  ifo.Materials.Coating.Philown    = 2e-4;              % reference ??

  % calculated for 123 K and 2128 nm following
  % Ghosh, et al (1994):  http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=317500
  ifo.Materials.Coating.Indexlown  = 1.43545; % calculated to 4 (RXA) and 6 (Paco) figures

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


  % factor of 2.5 added to simulate LNG modes - remove after new LNG code is added
  ifo.Optics.ITM.BeamRadius = 0.059;                     % m; 1/e^2 power radius
  ifo.Optics.ETM.BeamRadius = 0.084;                     % m; 1/e^2 power radius

  % coating layer optical thicknesses - mevans June 2008
  ifo.Optics.ITM.CoatingThicknessLown = 0.308;
  ifo.Optics.ITM.CoatingThicknessCap = 0.5;

  ifo.Optics.ETM.CoatingThicknessLown = 0.27;
  ifo.Optics.ETM.CoatingThicknessCap = 0.5;

  % Defining some additional fields for the structure for pygwinc compatibility
  ifo.Optics.ETM.Coating = ifo.Materials.Coating;
  ifo.Optics.ETM.Substrate = ifo.Materials.Substrate;
  ifo.Optics.ETM.MassRadius = ifo.Materials.MassRadius;
  ifo.Optics.ETM.MassThickness = ifo.Materials.MassThickness;

end

