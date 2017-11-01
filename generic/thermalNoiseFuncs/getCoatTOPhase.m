% [dphi_dT, dphi_TE, dphi_TR, rCoat] = ...
%   getCoatTOPhase(nIn, nOut, nLayer, dOpt, aLayer, bLayer)
%   returns coating reflection phase derivatives w.r.t. temperature
%
% nIn = refractive index of input medium (e.g., vacuum = 1)
% nOut = refractive index of output medium (e.g., SiO2 = 1.45231 @ 1064nm)
% nLayer = refractive index of each layer, ordered input to output (N x 1)
% dOpt   = optical thickness / lambda of each layer
%        = geometrical thickness * refractive index / lambda
% aLayer = change in geometrical thickness with temperature
%        = the effective thermal expansion coeffient of the coating layer
% bLayer = change in refractive index with temperature
%        = dn/dT 
%        = dd/dT - n * a
%
% dphi_dT = total thermo-optic phase derivative with respect to temperature
%         = dphi_TE + dphi_TR
% dphi_TE = thermo-elastic phase derivative (dphi / dT)
% dphi_TR = thermo-refractive phase derivative (dphi / dT)
% rCoat = amplitude reflectivity of coating (complex)
%
% Note about aLayer: on a SiO2 substrate,
% a_Ta2O5 ~ 3.5 * alpha_Ta2O5
% a_SiO2 ~ 2.3 * alpha_SiO2
%
% see getCoatTOPos for more information
% (see also T080101)

function [dphi_dT, dphi_TE, dphi_TR, rCoat] = ...
    getCoatTOPhase(nIn, nOut, nLayer, dOpt, aLayer, bLayer, sLayer)

  % vector of all refractive indexs
  nAll = [nIn; nLayer(:); nOut];
  
  % reflectivity of each interface
  r = (nAll(1:end-1) - nAll(2:end)) ./ (nAll(1:end-1) + nAll(2:end));
  
  % combine reflectivities
  rbar = zeros(size(r));
  ephi = zeros(size(r));
  
  ephi(end) = exp(-4i * pi * dOpt(end));
  rbar(end) = ephi(end) * r(end);
  for n = numel(dOpt):-1:1
    % one-way phase in this layer
    if n > 1
      ephi(n) = exp(-4i * pi * dOpt(n - 1));
    else
      ephi(n) = 1;
    end
    
    % accumulate reflecitivity
    rbar(n) = ephi(n) * (r(n) + rbar(n + 1)) / (1 + r(n) * rbar(n + 1));
  end
  
  % reflectivity derivatives
  dr_dphi = zeros(size(dOpt));
  for n = numel(dOpt):-1:1
    dr_dphi(n) = -1i * rbar(n + 1);
    for m = n:-1:1
      dr_dphi(n) = dr_dphi(n) * ephi(m) * ...
        (1 - r(m).^2) / (1 + r(m) * rbar(m + 1)).^2;
    end
  end
  
  % geometrical distances
  dGeo = dOpt ./ nLayer;
  
  % phase derivatives
  dphi_dd = 4 * pi * imag(dr_dphi / rbar(1));
  
  % thermo-refractive coupling
  dphi_TR = sum(dphi_dd .* (bLayer + sLayer .* nLayer) .* dGeo);
  
  % thermo-elastic
  dphi_TE = 4 * pi * sum(aLayer .* dGeo);
  
  % total
  dphi_dT = dphi_TR + dphi_TE;
  
  % coating reflectivity
  rCoat = rbar(1);
end