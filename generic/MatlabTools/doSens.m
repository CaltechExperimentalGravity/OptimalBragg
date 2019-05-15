%%Function that calculates a bunch of derivatives given a layer structure
%%by the calling function. 
%%Derivatives computed:
%%(dT/T) / (dl/l) for dl/l = +/-1% uniform -- @1064nm
%%(dR/R) / (dl/l) for dl/l = +/-1% uniform -- @532nm
%%(dT/T) / (dn1/n1) for dn1/n1 = +/-1% 
%%(dT/T) / (dn2/n2) for dn2/n2 = +/-1% 
%%(dT/T) / (daoi/aoi) for daoi/aoi = +/-1%

function [err] = doSens(n, L, lambda, aoi, Tp, Rp, varargin)

%The L here is actually product of l*n, so we wish to disentangle this...
L_phys = op2phys(L,n(2:end-1));

%First, vary L by +/- 1%...
[temp, ~] = multidiel100(n, (1.01*L_phys).*n(2:end-1), lambda, aoi, 'tm');
tempT = 1-abs(temp).^2;
tempR = abs(temp).^2;
sField_LP = ((27.46*abs(1+temp) - 0.01) / 0.01)/0.01;
err.Tp.coatLayer_plus = ((tempT - Tp)/(Tp))/0.01; %(dT/T) / (dL/L)
err.Rp.coatLayer_plus = ((tempR - Rp)/(Rp))/0.01; %(dT/T) / (dL/L)
if length(varargin)
    Rs = varargin{1};
    [temp, ~] = multidiel100(n, (1.01*L_phys).*n(2:end-1), lambda, aoi, 'te');
    tempR = abs(temp).^2;
    err.Rs.coatLayer_plus = ((tempR - Rs)/(Rs))/0.01; %(dT/T) / (dL/L)
end

[temp, ~] = multidiel100(n, (0.99*L_phys).*n(2:end-1), lambda, aoi, 'tm');
tempT = 1-abs(temp).^2;
tempR = abs(temp).^2;
sField_LM = ((27.46*abs(1+temp) - 0.01) / 0.01)/0.01;
err.Tp.coatLayer_minus = ((tempT - Tp)/(Tp))/0.01; %(dT/T) / (dL/L)
err.Rp.coatLayer_minus = ((tempR - Rp)/(Rp))/0.01; %(dT/T) / (dL/L)
if length(varargin)
    Rs = varargin{1};
    [temp, ~] = multidiel100(n, (0.99*L_phys).*n(2:end-1), lambda, aoi, 'te');
    tempR = abs(temp).^2;
    err.Rs.coatLayer_minus = ((tempR - Rs)/(Rs))/0.01; %(dT/T) / (dL/L)
end

%Next, vary n1 by +/- 1%...
n_temp = n;
for i=2:length(n)-1
    if mod(i,2)==0
        n_temp(i) = 1.01*n_temp(i);
    end
end
[temp, ~] = multidiel100(n_temp, L_phys.*n_temp(2:end-1), lambda, aoi, 'tm');
tempT = 1-abs(temp).^2;
tempR = abs(temp).^2;
sField_n1P = ((27.46*abs(1+temp) - 0.01) / 0.01)/0.01;
err.Tp.n1_plus = ((tempT - Tp)/(Tp))/0.01; 
err.Rp.n1_plus = ((tempR - Rp)/(Rp))/0.01; 
if length(varargin)
    Rs = varargin{1};
    [temp, ~] = multidiel100(n_temp, L_phys.*n_temp(2:end-1), lambda, aoi, 'te');
    tempR = abs(temp).^2;
    err.Rs.n1_plus = ((tempR - Rs)/(Rs))/0.01; 
end

n_temp = n;
for i=2:length(n)-1
    if mod(i,2)==0
        n_temp(i) = 0.99*n_temp(i);
    end
end
[temp, ~] = multidiel100(n_temp, L_phys.*n_temp(2:end-1), lambda, aoi, 'tm');
tempT = 1-abs(temp).^2;
tempR = abs(temp).^2;
sField_n1M = ((27.46*abs(1+temp) - 0.01) / 0.01)/0.01;
err.Tp.n1_minus = ((tempT - Tp)/(Tp))/0.01; 
err.Rp.n1_minus = ((tempR - Rp)/(Rp))/0.01; 
if length(varargin)
    Rs = varargin{1};
    [temp, ~] = multidiel100(n_temp, L_phys.*n_temp(2:end-1), lambda, aoi, 'te');
    tempR = abs(temp).^2;
    err.Rs.n1_minus = ((tempR - Rs)/(Rs))/0.01; 
end


%Next, vary n2 by +/- 1%...
n_temp = n;
for i=2:length(n)-1
    if mod(i,2)==1
        n_temp(i) = 1.01*n_temp(i);
    end
end
[temp, ~] = multidiel100(n_temp, L_phys.*n_temp(2:end-1), lambda, aoi, 'tm');
tempT = 1-abs(temp).^2;
tempR = abs(temp).^2;
sField_n2P = ((27.46*abs(1+temp) - 0.01) / 0.01)/0.01;
err.Tp.n2_plus = ((tempT - Tp)/(Tp))/0.01; 
err.Rp.n2_plus = ((tempR - Rp)/(Rp))/0.01; 
if length(varargin)
    Rs = varargin{1};
    [temp, ~] = multidiel100(n_temp, L_phys.*n_temp(2:end-1), lambda, aoi, 'te');
    tempR = abs(temp).^2;
    err.Rs.n2_plus = ((tempR - Rs)/(Rs))/0.01; 
end

n_temp = n;
for i=2:length(n)-1
    if mod(i,2)==1
        n_temp(i) = 0.99*n_temp(i);
    end
end
[temp, ~] = multidiel100(n_temp, L_phys.*n_temp(2:end-1), lambda, aoi, 'tm');
tempT = 1-abs(temp).^2;
tempR = abs(temp).^2;
sField_n2M = ((27.46*abs(1+temp) - 0.01) / 0.01)/0.01;
err.Tp.n2_minus = ((tempT - Tp)/(Tp))/0.01; 
err.Rp.n2_minus = ((tempR - Rp)/(Rp))/0.01; 
if length(varargin)
    Rs = varargin{1};
    [temp, ~] = multidiel100(n_temp, L_phys.*n_temp(2:end-1), lambda, aoi, 'te');
    tempR = abs(temp).^2;
    err.Rs.n2_minus = ((tempR - Rs)/(Rs))/0.01; 
end

err.totSurfField=sqrt(sField_n2M^2 + sField_n2P^2 + ...
					  sField_n1M^2 + sField_n1P^2 + ...
					  sField_LM^2  + sField_LP^2);

if length(varargin)>1
	aoi=varargin{2};
	%Finally, vary AoI by +/- 1%...
	[temp, ~] = multidiel100(n, L, lambda, 1.01*aoi, 'tm');
	tempT = 1-abs(temp).^2;
	tempR = abs(temp).^2;
	err.Tp.aoi_plus = ((tempT - Tp)/(Tp))/0.01; %(dT/T) / (dL/L)
	err.Rp.aoi_plus = ((tempR - Rp)/(Rp))/0.01; %(dT/T) / (dL/L)
	if length(varargin)
		Rs = varargin{1};
		[temp, ~] = multidiel100(n, L, lambda, 1.01*aoi, 'te');
		tempR = abs(temp).^2;
		err.Rs.aoi_plus = ((tempR - Rs)/(Rs))/0.01; %(dT/T) / (dL/L)
	end

	[temp, ~] = multidiel100(n, L, lambda, 0.99*aoi, 'tm');
	tempT = 1-abs(temp).^2;
	tempR = abs(temp).^2;
	err.Tp.aoi_minus = ((tempT - Tp)/(Tp))/0.01; %(dT/T) / (dL/L)
	err.Rp.aoi_minus = ((tempR - Rp)/(Rp))/0.01; %(dT/T) / (dL/L)
	if length(varargin)
		Rs = varargin{1};
		[temp, ~] = multidiel100(n, L, lambda, 0.99*aoi, 'te');
		tempR = abs(temp).^2;
		err.Rs.aoi_minus = ((tempR - Rs)/(Rs))/0.01; %(dT/T) / (dL/L)
	end
end

