% [dc, Cc, Kc, aSub] = getCoatAvg(ifo, dOpt)
%   get coating average properties
%   (see getCoatTOPos for example usage)
%
% ifo  = parameter struct from IFOmodel.m
% dOpt = optical thickness / lambda of each layer
%      = geometrical thickness * refractive index / lambda
%
% dc = total thickness (meters)
% Cc = heat capacity
% Kc = thermal diffusivity
% aSub = effective substrate thermal expansion (weighted by heat capacity)

function [dc, Cc, Kc, aSub] = getCoatAvg(ifo, dOpt)
  
  % coating parameters
  pS = ifo.Materials.Substrate;
  pC = ifo.Materials.Coating;
  
  alphaS = pS.MassAlpha;
  C_S = pS.MassCM * pS.MassDensity;
  sigS = pS.MirrorSigma;
  
  C_L = pC.CVlown;
  K_L = pC.ThermalDiffusivitylown;
  
  C_H = pC.CVhighn;
  K_H = pC.ThermalDiffusivityhighn;
  
  % compute refractive index, effective alpha and beta
  [~, ~, ~, dLayer] = getCoatLayers(ifo, dOpt);
  
  % heat capacity
  dc = sum(dLayer);
  dL = sum(dLayer(1:2:end));
  dH = sum(dLayer(2:2:end));
  Cc = (C_L * dL + C_H * dH) / dc;
  
  % thermal diffusivity
  KinvL = 1 / K_L;
  KinvH = 1 / K_H;
  Kc = dc / (KinvL * dL + KinvH * dH);
  
  % effective substrate thermal expansion
  aSub = 2 * alphaS * (1 + sigS) * Cc / C_S;
  
end
