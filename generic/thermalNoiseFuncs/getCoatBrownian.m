% SbrZ  = getCoatBrownian(f, ifo, opticName)
% SbrZ  = getCoatBrownian(f, ifo, wBeam, dOpt)
%  returns the coating brownian noise for a given collection of
%  coating layers.  The layers are assumed to be alernating low-n
%  high-n layers, with low-n first.
%
% f = frequency vector in Hz
% ifo = parameter struct from IFOmodel.m
%
% opticName = name of the Optic struct to use for wBeam and dOpt
%   wBeam = ifoArg.Optics.(opticName).BeamRadius
%   dOpt = ifoArg.Optics.(opticName).CoatLayerOpticalThickness
%
% wBeam = beam radius (at 1 / e^2 power)
% dOpt = coating layer thickness vector (Nlayer x 1)
%      = the optical thickness, normalized by lambda, of each coating layer.
%
% SbrZ = Brownian noise spectra for one mirror in m^2 / Hz
%
% adapted from bench62
% based on Harry et al., Class Quant Grav 24 (2007) 405-415

function SbrZ = getCoatBrownian(f, ifo, wBeam, dOpt)
  
  % check arguments
  if nargin < 4
    % wBeam has been used as the name of an optic
    [wBeam, dOpt] = getCoatParFromName(ifo, wBeam);
  end
  
  % Constants
  if isfield(ifo.Materials.Substrate, 'Temp')
    % use substrate temperature
    subTemp = ifo.Materials.Substrate.Temp;
  else
    % use system temperature
    subTemp = ifo.Constants.Temp;
  end
  
  kBT    = ifo.Constants.kB * subTemp;
  lambda = ifo.Laser.Wavelength;

  Ysub     = ifo.Materials.Substrate.MirrorY;
  sigmasub = ifo.Materials.Substrate.MirrorSigma;

  Yhighn     = ifo.Materials.Coating.Yhighn;
  sigmahighn = ifo.Materials.Coating.Sigmahighn;
  phihighn   = ifo.Materials.Coating.Phihighn;
  nH         = ifo.Materials.Coating.Indexhighn;

  Ylown     = ifo.Materials.Coating.Ylown;
  sigmalown = ifo.Materials.Coating.Sigmalown;
  philown   = ifo.Materials.Coating.Philown;
  nL        = ifo.Materials.Coating.Indexlown;

  % compute thickness of each material in the coating
  dlown  = sum(dOpt(1:2:end)) * lambda / nL;
  dhighn = sum(dOpt(2:2:end)) * lambda / nH;
  dCoat  = dlown + dhighn;

% for debugging, this is a rough but direct estimate
%   Llown = dlown * philown;
%   Lhighn = dhighn * phihighn;
%   Lsum = Llown + Lhighn;
%   Lall = [Lsum, Llown, Lhighn]

  %%%%%%%%%%%%%%%%% this part is directly from bench62 %%%%%%%%%%%%%%%%%
  Yperp = dCoat/(dhighn/Yhighn+dlown/Ylown);
  phiperp = Yperp/dCoat*(dlown*philown/Ylown + dhighn*phihighn/Yhighn);
  Ypara = 1/dCoat*(Yhighn*dhighn + Ylown*dlown);
  phipara = 1/(dCoat*Ypara)*(Ylown*philown*dlown + Yhighn*phihighn*dhighn);

  % This is a kludge, the real formula is very complicted but this
  % average works really well
  sigma1 = 1/2*(sigmahighn+sigmalown);

  % This is exact
  sigma2 = (sigmahighn*Yhighn*dhighn+sigmalown*Ylown*dlown)/(Yhighn*dhighn+Ylown*dlown);

  % Brownian contribution to coating thermal noise, low Poisson ratio limit
  % cITM = dITM/(pi*wITM^2)*(Ypara/Ysub^2*phipara+phiperp/Yperp);
  % cETM = dCoat/(pi*wETM^2)*(Ypara/Ysub^2*phipara+phiperp/Yperp);

  % Brownian contribution to coating thermal noise, full formula
  cbnoise = dCoat*(1-sigmasub^2)/(pi*wBeam^2)*((1/(Yperp*(1-sigmasub^2))-...
    2*sigma2^2*Ypara/(Yperp^2*(1-sigmasub^2)*(1-sigma1)))*phiperp +...
    Ypara*sigma2*(1-2*sigmasub)/(Yperp*Ysub*(1-sigma1)*(1-sigmasub))*(phipara-phiperp)+...
    Ypara*(1+sigmasub)*(1-2*sigmasub)^2/(Ysub^2*(1-sigma1^2)*(1-sigmasub))*phipara);

  % noise power spectrum
  SbrZ = 4 * kBT * cbnoise ./ (2 * pi * f);
  
end