% [nLayer, aLayer, bLayer, dLayer] = getCoatLayers(ifo, dOpt)
%   get layer vectors for refractive index, effective alpha and beta
%   and geometrical thickness
%   (see getCoatTOPos for example usage)
%
% ifo    = parameter struct from IFOmodel.m
% dOpt   = optical thickness / lambda of each layer
%        = geometrical thickness * refractive index / lambda
%
% nLayer = refractive index of each layer, ordered input to output (N x 1)
% aLayer = change in geometrical thickness with temperature
%        = the effective thermal expansion coeffient of the coating layer
% bLayer = change in refractive index with temperature
%        = dn/dT
% dLayer = geometrical thickness of each layer
% sLayer = Yamamoto thermo-refractive correction
%        = alpha * (1 + sigma) / (1 - sigma)

function [nLayer, aLayer, bLayer, dLayer, sLayer] = getCoatLayers(ifo, dOpt)
   
  % coating parameters
  lambda = ifo.Laser.Wavelength;
    
  pS   = ifo.Materials.Substrate;
  pC   = ifo.Materials.Coating;
  
  Y_S  = pS.MirrorY;
  sigS = pS.MirrorSigma;
  
  alphaL = pC.Alphalown;
  betaL  = pC.Betalown;
  Y_L    = pC.Ylown;
  sigL   = pC.Sigmalown;
  nL     = pC.Indexlown;
  
  alphaH = pC.Alphahighn;
  betaH  = pC.Betahighn;
  Y_H    = pC.Yhighn;
  sigH   = pC.Sigmahighn;
  nH     = pC.Indexhighn;

  Nlayer = numel(dOpt);

  try
      CoatType = ifo.Materials.Coating.Type;
  catch
      CoatType = 'standard';
  end

  % compute effective alpha
  aLayer = zeros(Nlayer, 1);
  switch CoatType
    case 'standard'
      aLayer(1:2:end) = alphaL * getExpansionRatio(Y_L, sigL, Y_S, sigS);
      aLayer(2:2:end) = alphaH * getExpansionRatio(Y_H, sigH, Y_S, sigS);
    case 'AlGaAs'
      aLayer(2:2:end) = alphaL * getExpansionRatio(Y_L, sigL, Y_S, sigS);
      aLayer(1:2:end) = alphaH * getExpansionRatio(Y_H, sigH, Y_S, sigS);
  end

  % and beta
  bLayer = zeros(Nlayer, 1);
  switch CoatType
    case 'standard'
      bLayer(1:2:end) = betaL;
      bLayer(2:2:end) = betaH;
    case 'AlGaAs'
      bLayer(2:2:end) = betaL;
      bLayer(1:2:end) = betaH;
  end

  % and refractive index
  nLayer = zeros(Nlayer, 1);
  switch CoatType
    case 'standard'
      nLayer(1:2:end) = nL;
      nLayer(2:2:end) = nH;
    case 'AlGaAs'
      nLayer(2:2:end) = nL;
      nLayer(1:2:end) = nH;
  end

  % and geometrical thickness
  dLayer = lambda * dOpt ./ nLayer;

  % and sigma correction
  sLayer = zeros(Nlayer, 1);
  switch CoatType
    case 'standard'
      sLayer(1:2:end) = alphaL * (1 + sigL) / (1 - sigL);
      sLayer(2:2:end) = alphaH * (1 + sigH) / (1 - sigH);
    case 'AlGaAs'
      sLayer(2:2:end) = alphaL * (1 + sigL) / (1 - sigL);
      sLayer(1:2:end) = alphaH * (1 + sigH) / (1 - sigH);
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y_C and sigC are for the coating material (can also be substrate)
% Y_S and sigS are for the substrate material
%
function ce = getExpansionRatio(Y_C, sigC, Y_S, sigS)
  
  ce = ((1 + sigS) / (1 - sigC)) ...
    * ( ((1 + sigC) / (1 + sigS)) + (1 - 2 * sigS) * Y_C / Y_S );
  
end
