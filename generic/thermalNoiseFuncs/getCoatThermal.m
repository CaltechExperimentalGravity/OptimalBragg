% [SsurfT, rdel, gFC, xi] = getCoatThermal(f, ifo, wBeam)
%   returns the thermal noise spectra for a surface layer
%
% f = frequency vector in Hz
% ifo = parameter struct from IFOmodel.m
%
% wBeam = beam radius (at 1 / e^2 power)
%       = usually ifo.Optics.ITM.BeamRadius or ifo.Optics.ETM.BeamRadius
%
% SsurfT = power spectra of thermal fluctuations in K^2 / Hz
% rdel = thermal diffusion length at each frequency in m

function [SsurfT, rdel] = getCoatThermal(f, ifo, wBeam)
  
  % parameter extraction
  kBT2 = ifo.Constants.kB * ifo.Materials.Substrate.Temp^2;
  
  pS = ifo.Materials.Substrate;
  C_S = pS.MassCM * pS.MassDensity;
  K_S = pS.MassKappa;

  % omega
  w = 2 * pi * f;
  
  % thermal diffusion length
  rdel = sqrt(2 * K_S ./ (C_S * w));
  
  % noise equation
  SsurfT = 4 * kBT2 ./ (pi * w * C_S .* rdel * wBeam^2);
end
