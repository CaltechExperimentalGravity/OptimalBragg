function [StoZ, SteZ, StrZ, T] = getCoatThermoOptic(f, ifo, wBeam, dOpt)

% [StoZ, SteZ, StrZ, T] = getCoatThermoOptic(f, ifo, opticName);
% [StoZ, SteZ, StrZ, T] = getCoatThermoOptic(f, ifo, wBeam, dOpt);
%   return power spectra of coating thermo-optic noise for a single optic
%
% f = frequency vector in Hz
% ifo = parameter struct from IFOmodel.m
% opticName = name of the Optic struct to use for wBeam and dOpt
%   wBeam = ifoArg.Optics.(opticName).BeamRadius
%   dOpt = ifoArg.Optics.(opticName).CoatLayerOpticalThickness
%
% wBeam = beam radius (at 1 / e^2 power)
% dOpt   = optical thickness / lambda of each layer
%        = geometrical thickness * refractive index / lambda
%
% StoZ = power spectra of apparent mirror position in m^2 / Hz
% SteZ = thermo-elastic componenet of StoZ
% StrZ = thermo-refractive componenet of StoZ
% T = coating power transmission, made available as a cross-check

  
  % check arguments
  if nargin < 4
    % wBeam has been used as the name of an optic
    [wBeam, dOpt] = getCoatParFromName(ifo, wBeam);
  end
  
  % compute coefficients
  [dTO, dTR, dTE, T, ~] = getCoatTOPos(ifo, wBeam, dOpt);
  
  % compute correction factors
  gTO = getCoatThickCorr(f, ifo, dOpt, dTE, dTR);
  gTE = getCoatThickCorr(f, ifo, dOpt, dTE, 0);
  gTR = getCoatThickCorr(f, ifo, dOpt, 0, dTR);
  
  % compute thermal source spectrum
  SsurfT = getCoatThermal(f, ifo, wBeam);
  
  StoZ = SsurfT .* gTO * dTO^2;
  SteZ = SsurfT .* gTE * dTE^2;
  StrZ = SsurfT .* gTR * dTR^2;
end
