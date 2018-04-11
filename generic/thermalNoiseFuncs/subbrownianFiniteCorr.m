function [cftm, aftm] = subbrownianFiniteCorr(ifo, opticName)
% subbrownianFiniteCorr - Estimate amplitude coefficient of
% mirror thermal noise contribution for finite-size test masses.
%
% [cftm, aftm] = subbrownianFiniteCorr(ifo, opticName)
% cftm = finite mirror correction factor
% aftm = amplitude coefficient for thermal noise:
%        thermal noise contribution to displacement noise is
%         S_x(f) = (8 * kB * T / (2*pi*f)) * Phi(f) * aftm
%
% Equation references to Bondu, et al. Physics Letters A 246 (1998)
% 227-236 (hereafter BHV) and Liu and Thorne gr-qc/0002055 (hereafter LT)

  % get some numbers
  a = ifo.Materials.MassRadius;
  h = ifo.Materials.MassThickness;
  w = ifo.Optics.(opticName).BeamRadius;
  Y = ifo.Materials.Substrate.MirrorY;
  sigma = ifo.Materials.Substrate.MirrorSigma;
  zeta = ifo.Constants.BesselZeros;
  
  % do the work
  j0m = besselj(0,zeta);
  r0 = w / sqrt(2);			% LT uses e-folding of power
  km = zeta/a;

  Qm = exp(-2*km*h);			% LT eq. 35a

  Um = (1-Qm).*(1+Qm)+4*h*km.*Qm;
  Um = Um./((1-Qm).^2-4*(km*h).^2.*Qm);	% LT 53 (BHV eq. btwn 29 & 30)

  x = exp(-(zeta*r0/a).^2/4);
  s = sum(x./(zeta.^2.*j0m));		% LT 57

  x2 = x.*x;
  U0 = sum(Um.*x2./(zeta.*j0m.^2));
  U0 = U0*(1-sigma)*(1+sigma)/(pi*a*Y);	% LT 56 (BHV eq. 3)

  p0 = 1/(pi*a^2);			% LT 28
  DeltaU = (pi*h^2*p0)^2;
  DeltaU = DeltaU + 12*pi*h^2*p0*sigma*s;
  DeltaU = DeltaU + 72*(1-sigma)*s^2;
  DeltaU = DeltaU*a^2/(6*pi*h^3*Y);	% LT 54

  aftm = DeltaU + U0;			% LT 58 (eq. following BHV 31)

  % amplitude coef for infinite TM
  %   factored out: (8 * kB * T * Phi) / (2 * pi * f)
  aitm = (1 - sigma^2) / (2 * sqrt(2 * pi) * Y * r0); % LT 59
  
  % finite mirror correction
  cftm = aftm / aitm;