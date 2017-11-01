% [dTO, dTR, dTE, T, R] = getCoatTOPos(ifo, opticName)
% [dTO, dTR, dTE, T, R] = getCoatTOPos(ifo, wBeam, dOpt)
%   returns mirror position derivative wrt thermal fluctuations
%
% ifo  = parameter struct from IFOmodel.m
% opticName = name of the Optic struct to use for wBeam and dOpt
%   wBeam = ifoArg.Optics.(opticName).BeamRadius
%   dOpt = ifoArg.Optics.(opticName).CoatLayerOpticalThickness
%
% wBeam = beam radius, for finite mirror correction (at 1 / e^2 power)
% dOpt = optical thickness / lambda of each layer
%      = geometrical thickness * refractive index / lambda
%
% dTO = total thermo-optic dz/dT
% dTR = thermo-refractive dz/dT
% dTE = thermo-elastic dz/dT
%
% compute thermal fluctuations with getCoatThermal
% (see also T080101)

function [dTO, dTR, dTE, T, R] = getCoatTOPos(ifo, wBeam, dOpt)
  
  if nargin < 3
    % wBeam has been used as the name of an optic
    [wBeam, dOpt] = getCoatParFromName(ifo, wBeam);
  end
  
  % parameters
  lambda = ifo.Laser.Wavelength;
  nS = ifo.Materials.Substrate.RefractiveIndex;  
  
  % compute refractive index, effective alpha and beta
  [nLayer, aLayer, bLayer, dLayer, sLayer] = getCoatLayers(ifo, dOpt);

  % compute coating average parameters
  [dc, Cc, Kc, aSub] = getCoatAvg(ifo, dOpt);
  
  % compute reflectivity and parameters
  [dphi_dT, dphi_TE, dphi_TR, rCoat] = ...
    getCoatTOPhase(1, nS, nLayer, dOpt, aLayer, bLayer, sLayer);
  R = abs(rCoat).^2;
  T = 1 - R;
  
  % for debugging
  %disp(sprintf('R = %.3f, T = %.0f ppm', R, 1e6 * T))
  
  % convert from phase to meters, subtracting substrate
  dTR = dphi_TR * lambda / (4 * pi);
  dTE = dphi_TE * lambda / (4 * pi) - aSub * dc;

  % mirror finite size correction
  Cfsm = getCoatFiniteCorr(ifo, wBeam, dOpt);
  dTE = dTE * Cfsm;
  
  % add TE and TR effects (sign is already included)
  dTO = dTE + dTR;
  
end
