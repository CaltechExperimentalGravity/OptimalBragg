%%Function that calculates a bunch of derivatives given a layer structure
%%by the calling function. 
%%Derivatives computed:
%%(dT/T) / (dl/l) for dl/l = +/-1% uniform -- @1064nm
%%(dT/T) / (dn1/n1) for dn1/n1 = +/-1% 
%%(dT/T) / (dn2/n2) for dn2/n2 = +/-1% 
%%(dT/T) / (daoi/aoi) for daoi/aoi = +/-5%

function [err] = doSens(n, L, lambda, aoi, Tp, Rp)

%The L here is actually product of l*n, so we wish to disentangle this...
L_phys = op2phys(L,n(2:end-1));

%First, vary L by +/- 1%...
[temp, ~] = multidiel1(n, (1.01*L_phys).*n(2:end-1), lambda, aoi, 'tm');
tempT = 1-abs(temp).^2;
tempR = abs(temp).^2;
sField_LP = ((27.46*abs(1+temp) - 0.01) / 0.01)/0.01;
err.Tp.coatLayer_plus = ((tempT - Tp)/(Tp))/0.01; %(dT/T) / (dL/L)
err.Rp.coatLayer_plus = ((tempR - Rp)/(Rp))/0.01; %(dT/T) / (dL/L)

[temp, ~] = multidiel1(n, (0.99*L_phys).*n(2:end-1), lambda, aoi, 'tm');
tempT = 1-abs(temp).^2;
tempR = abs(temp).^2;
sField_LM = ((27.46*abs(1+temp) - 0.01) / 0.01)/0.01;
err.Tp.coatLayer_minus = ((tempT - Tp)/(Tp))/0.01; %(dT/T) / (dL/L)
err.Rp.coatLayer_minus = ((tempR - Rp)/(Rp))/0.01; %(dT/T) / (dL/L)

%Next, vary n1 by +/- 1%...
n_temp = n;
for i=2:length(n)-1
    if mod(i,2)==0
        n_temp(i) = 1.01*n_temp(i);
    end
end
[temp, ~] = multidiel1(n_temp, L_phys.*n_temp(2:end-1), lambda, aoi, 'tm');
tempT = 1-abs(temp).^2;
tempR = abs(temp).^2;
sField_n1P = ((27.46*abs(1+temp) - 0.01) / 0.01)/0.01;
err.Tp.n1_plus = ((tempT - Tp)/(Tp))/0.01; 
err.Rp.n1_plus = ((tempR - Rp)/(Rp))/0.01; 

n_temp = n;
for i=2:length(n)-1
    if mod(i,2)==0
        n_temp(i) = 0.99*n_temp(i);
    end
end
[temp, ~] = multidiel1(n_temp, L_phys.*n_temp(2:end-1), lambda, aoi, 'tm');
tempT = 1-abs(temp).^2;
tempR = abs(temp).^2;
sField_n1M = ((27.46*abs(1+temp) - 0.01) / 0.01)/0.01;
err.Tp.n1_minus = ((tempT - Tp)/(Tp))/0.01; 
err.Rp.n1_minus = ((tempR - Rp)/(Rp))/0.01; 


%Next, vary n2 by +/- 1%...
n_temp = n;
for i=2:length(n)-1
    if mod(i,2)==1
        n_temp(i) = 1.01*n_temp(i);
    end
end
[temp, ~] = multidiel1(n_temp, L_phys.*n_temp(2:end-1), lambda, aoi, 'tm');
tempT = 1-abs(temp).^2;
tempR = abs(temp).^2;
sField_n2P = ((27.46*abs(1+temp) - 0.01) / 0.01)/0.01;
err.Tp.n2_plus = ((tempT - Tp)/(Tp))/0.01; 
err.Rp.n2_plus = ((tempR - Rp)/(Rp))/0.01; 

n_temp = n;
for i=2:length(n)-1
    if mod(i,2)==1
        n_temp(i) = 0.99*n_temp(i);
    end
end
[temp, ~] = multidiel1(n_temp, L_phys.*n_temp(2:end-1), lambda, aoi, 'tm');
tempT = 1-abs(temp).^2;
tempR = abs(temp).^2;
sField_n2M = ((27.46*abs(1+temp) - 0.01) / 0.01)/0.01;
err.Tp.n2_minus = ((tempT - Tp)/(Tp))/0.01; 
err.Rp.n2_minus = ((tempR - Rp)/(Rp))/0.01; 

err.totSurfField=sqrt(sField_n2M^2 + sField_n2P^2 + ...
					  sField_n1M^2 + sField_n1P^2 + ...
					  sField_LM^2  + sField_LP^2);

%Finally, vary AoI by +/- 1%...
[temp, ~] = multidiel1(n, L, lambda, 1.05*aoi, 'tm');
tempT = 1-abs(temp).^2;
tempR = abs(temp).^2;
err.Tp.aoi_plus = ((tempT - Tp)/(Tp))/0.05; %(dT/T) / (dL/L)
err.Rp.aoi_plus = ((tempR - Rp)/(Rp))/0.05; %(dT/T) / (dL/L)

[temp, ~] = multidiel1(n, L, lambda, 0.95*aoi, 'tm');
tempT = 1-abs(temp).^2;
tempR = abs(temp).^2;
err.Tp.aoi_minus = ((tempT - Tp)/(Tp))/0.05; %(dT/T) / (dL/L)
err.Rp.aoi_minus = ((tempR - Rp)/(Rp))/0.05; %(dT/T) / (dL/L)
