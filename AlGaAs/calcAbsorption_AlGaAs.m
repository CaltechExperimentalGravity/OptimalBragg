% Function for calculating the Absorption given an electric field profile
% Made to work together with calcEField_AlGaAs.m
% Input arguments:
%   - Esq       = NORMALIZED electric field squared as a function of distance in a coating
%   - L         = PHYSICAL thickness of coating layers
%   - nPts      = # of points inside each layer at which field is to be evaluated
%   - alphaOdd  = Absorption coefficient of all odd layers [m ^-1] (top layer is assumed layer #1)
%   - alphaEven = Absorption coefficient of all even layers [m ^-1] (top layer is assumed layer #1)
% Output variables:
%   - absorp    = Absorption of coating [ppm]

function absorp = calcAbsorption_AlGaAs(Esq, L, nPts, alphaOdd, alphaEven)
	dL = L ./nPts;
	dz = [];
	%Define the grid for rectangular integration
	for ii = 1:length(L)
	       temp = repmat(dL(ii),nPts,1);
       	       dz = vertcat(dz,temp);
	end
	alph = ones(length(L)*nPts,1);
	alph(1:2:end) = alphaOdd;
	alph(2:2:end) = alphaEven;
	absorp = 1e6 * sum(alph .* dz .* Esq);
end
