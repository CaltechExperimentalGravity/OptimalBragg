function Cfsm = getCoatFiniteCorr(ifo, wBeam, dOpt)

% Cfsm = getCoatFiniteCorr(ifo, opticName)
% Cfsm = getCoatFiniteCorr(ifo, wBeam, dOpt)
%   finite mirror size correction
%     Uses correction factor from PLA 2003 vol 312 pg 244-255
%     "Thermodynamical fluctuations in optical mirror coatings"
%     by V. B. Braginsky and S. P. Vyatchanin
%      http://arxiv.org/abs/cond-mat/0302617
%
% ifo = parameter struct from IFOmodel.m
% opticName = name of the Optic struct to use for wBeam and dOpt
%   wBeam = ifoArg.Optics.(opticName).BeamRadius
%   dOpt = ifoArg.Optics.(opticName).CoatLayerOpticalThickness
%
% wBeam = beam radius (at 1 / e^2 power)
% dOpt   = optical thickness / lambda of each layer
%        = geometrical thickness * refractive index / lambda
%
% (see getCoatTOPos for example usage)
%
% version 1 by Sam Wald, 2008


  % check arguments
  if nargin < 3
    % wBeam has been used as the name of an optic
    [wBeam, dOpt] = getCoatParFromName(ifo, wBeam);
  end
  
  % parameter extraction
  R      = ifo.Materials.MassRadius;      %substrate radius
  H      = ifo.Materials.MassThickness;   %substrate thickness
  lambda = ifo.Laser.Wavelength;
  zeta   = ifo.Constants.BesselZeros;  % zeros of 1st order bessel function (J1)
    
  pS = ifo.Materials.Substrate;
  pC = ifo.Materials.Coating;
  
  alphaS = pS.MassAlpha;
  C_S    = pS.MassCM * pS.MassDensity;
  Y_S    = pS.MirrorY;
  sigS   = pS.MirrorSigma;
  
  alphaL = pC.Alphalown;
  C_L    = pC.CVlown;
  Y_L    = pC.Ylown;
  sigL   = pC.Sigmalown;
  nL     = pC.Indexlown;
  
  alphaH = pC.Alphahighn;
  C_H    = pC.CVhighn;
  Y_H    = pC.Yhighn;
  sigH   = pC.Sigmahighn;
  nH     = pC.Indexhighn;

  % coating sums
  dL = lambda * sum(dOpt(1:2:end)) / nL;
  dH = lambda * sum(dOpt(2:2:end)) / nH;
  dc = dH + dL;
  
  % AVERAGE SPECIFIC HEAT (simple volume average for coating)
  Cf = (C_L * dL + C_H * dH) / dc;
  Cr = Cf / C_S;

  % COATING AVERAGE VALUE X = ALPHAF*(1+POISSONf)/(1-POISSONf) avg
  xxL = alphaL * (1 + sigL) / (1 - sigL);
  xxH = alphaH * (1 + sigH) / (1 - sigH);
  Xf = (xxL * dL + xxH * dH) / dc;
  Xr = Xf / alphaS;

  % COATING AVERAGE VALUE Y = ALPHAF* YOUNGSF/(1-POISSONF) avg
  yyL = alphaL * Y_L / (1 - sigL);
  yyH = alphaH * Y_H / (1 - sigH);
  Yf = (yyL * dL + yyH * dH) / dc;
  Yr = Yf / (alphaS * Y_S);

  % COATING AVERAGE VALUE Z = 1/(1-POISSONF) avg
  zzL = 1 / (1 - sigL);
  zzH = 1 / (1 - sigH);
  Zf = (zzL * dL + zzH * dH) / dc;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FINITE SIZE CORRECTION CALCULATION

  % beam size parameter used by Braginsky
  r0 = wBeam / sqrt(2);

  % values of J0 at zeros of J1
  j0m = besselj(0, zeta);

  % between eq 77 and 78
  km = zeta / R;
  Qm = exp(-2 * km * H);
  pm = exp(-km.^2 * r0^2 / 4) ./ j0m;  % left out factor of pi * R^2 in denominator

  % eq 88
  Lm = Xr - Zf * (1 + sigS) + (Yr * (1 - 2 * sigS) + Zf - 2 * Cr) * ...
    (1 + sigS) * (1 - Qm).^2 ./ ((1 - Qm).^2 - 4 * km.^2 * H^2 .* Qm);

  % eq 90 and 91
  S1 = (12 * R^2 / H^2) * sum(pm ./ zeta.^2);
  S2 = sum(pm.^2 .* Lm.^2);
  P = (Xr - 2 * sigS * Yr - Cr + S1 * (Cr - Yr * (1 - sigS)))^2 + S2;
  
  % eq 60 and 70
  LAMB = -Cr + (Xr / (1 + sigS) + Yr * (1 - 2 * sigS)) / 2;

  % eq 92
  Cfsm = sqrt((r0^2 * P) / (2 * R^2 * (1 + sigS)^2 * LAMB^2));
  
end
