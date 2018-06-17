% Function for calculating the E-field squared inside first nLayers of a dielectric stack coating
% E-field can then be used to compute absorption
% Input arguments:
%   - L         = vector of PHYSICAL layer thicknesses [m]
%   - n         = vector of refractive indices
%   - nLayers   = number of layers of coating in which to calculate E field
%   - lam       = wavelength of light [m]
%   - theta     = angle of incidence [degrees]
%   - pol       = polarization, 's' or 'p'
%   - nPts      = # of points inside each layer at which field is to be evaluated
% Output variables:
%   - z         = Vector at which E field is evaluated [m]
%   - E_prof    = Normalized E-field squared as a function of penetration depth, evaluated at z.

function [z, E_prof] = calcEfield_AlGaAs(L, n, nLayers, lam, theta, pol, nPts)
	function alpha = arrayTheta(n, theta0)
		alpha = ones(length(n)-1,1);
		alpha(1) = deg2rad(theta0);
		for ii=1:length(n)-1
			alpha(ii+1) = asin(n(ii)*sin(alpha(ii))/n(ii+1));
		end
		alpha = alpha(2:end)';
	end
	function mm = M_i(b_i, qq_i)
		mm = [cos(b_i), 1j*sin(b_i)/qq_i; 1j*sin(b_i)*qq_i, cos(b_i)];
	end

	function qq = q_i(n_i, theta_i)
		if strcmp(pol,'s')
			qq = n_i * cos(theta_i);
		elseif strcmp(pol,'p')
			qq = n_i / cos(theta_i);
		end
	end

	function bb = beta_i(tt_i, nn_i, hh_i, lam)
		bb = 2*pi*cos(tt_i)*nn_i*hh_i/lam;
	end

	function Epk = E0pk(Mtot)
		q0 = q_i(n(1), theta);
		qSub = q_i(n(end), qAngle);
		Epk = 0.25 * (abs(Mtot(1,1) + Mtot(2,2)*qSub/q0)^2 + abs(Mtot(2,1)/q0/1j + Mtot(1,2)*qSub/1j)^2);	
	end

	function dh = delta_h(bb_i, qq_i)
		dh = M_i(bb_i, -1.*qq_i);
	end

	%Calculate the array of angles, in radians
	angles = arrayTheta(n, theta);
	qAngle = angles(end);
	angles = angles(1:end-1);
	%Evaluate the characteristic matrix for the stack
	Mt = eye(2);
	for ii=1:length(L)
		Mt = Mt * M_i(beta_i(angles(ii),n(ii+1),L(ii),lam), q_i(n(ii+1),angles(ii)));
	end
	Mtotz = Mt;
	%Multiply by inverse of infinitesimal stacks
	E_profile = zeros(length(nLayers)*nPts,1);
	z = zeros(nLayers*nPts,1);
	Z = 0;
	correction = 1;
	qSub = q_i(n(end),qAngle);
	for ii=1:nLayers
		n_i = n(ii+1);
		dL = L(ii)/nPts;
		theta_i = angles(ii);
		if strcmp(pol,'p')
			correction=(cos(theta)/cos(theta_i))^2;
		end
		for jj=1:nPts
			Z = Z + dL;
			z((ii-1)*nPts + jj) = Z;
			Mtotz = delta_h(beta_i(theta_i, n_i, dL,lam), q_i(n_i, theta_i)) * Mtotz;
			E_profile((ii-1)*nPts+jj)=correction*abs(Mtotz(1,1)^2 + abs(qSub*Mtotz(1,2)/1j)^2);
		end	
	end	       
	
	E_prof = E_profile/E0pk(Mt);
end