%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gTC = getCoatThickCorr(f, ifo, dOpt, dTE, dTR)
%   finite coating thickness correction
%     Uses correction factor from T080101, "Thick Coating Correction" (Evans)
%
% (see getCoatThermoOptic for example usage)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For comparison in the bTR = 0 limit, the
% equation from Fejer (PRD70, 2004)
% modified so that gFC -> 1 as xi -> 0
%  gTC = (2 ./ (R * xi.^2)) .* (sh - s + R .* (ch - c)) ./ ...
%    (ch + c + 2 * R * sh + R^2 * (ch - c));
% which in the limit of xi << 1 becomes
%  gTC = 1 - xi * (R - 1 / (3 * R));

function gTC = getCoatThickCorr(f, ifo, dOpt, dTE, dTR)

  % parameter extraction
  pS = ifo.Materials.Substrate;
  Cs = pS.MassCM * pS.MassDensity;
  Ks = pS.MassKappa;
  
  % compute coating average parameters
  [d_coat, Cc, Kc] = getCoatAvg(ifo, dOpt);
  
  % R and xi (from T080101, Thick Coating Correction)
  w  = 2 * pi * f;
  R  = sqrt(Cc * Kc / (Cs * Ks));
  xi = d_coat * sqrt(2 * w * Cc / Kc);
  
  % trig functions of xi
  s  = sin(xi);
  c  = cos(xi);
  sh = sinh(xi);
  ch = cosh(xi);
  
  % pR and pE (dTR = -\bar{\beta} lambda, dTE = \Delta \bar{\alpha} d)
  pR = dTR / (dTR + dTE);
  pE = dTE / (dTR + dTE);
  
  % various parts of gTC
  g0 = 2 * (sh - s) + 2 * R * (ch - c);
  g1 = 8 * sin(xi / 2) .* (R * cosh(xi / 2) + sinh(xi / 2));
  g2 = (1 + R^2) * sh + (1 - R^2) * s + 2 * R * ch;
  gD = (1 + R^2) * ch + (1 - R^2) * c + 2 * R * sh;

  % and finally, the correction factor
  gTC = (pE^2 * g0 + pE * pR * xi .* g1 + pR^2 * xi.^2 .* g2) ./ (R * xi.^2 .* gD);

end
