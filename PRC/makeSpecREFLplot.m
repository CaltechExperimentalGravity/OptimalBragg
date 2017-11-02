%function to make spectral reflectivity plots for coating optimization tool

function [Rp,Tp,Rs,Ts] = makeSpecREFLplot(dataFile, firstLayer)
clear TNout;
addpath('../');
load(dataFile);
load dispersion_revised.mat;

lambda_0 = 1064e-9;
lambda = linspace(0.4,1.6,2200);

aoi = TNout.aoi;
L= TNout.L;
no_of_stacks = floor(length(L)/2);

na=1.000;
%Initialize arrays
Rp = zeros(length(lambda),1);
Tp = zeros(length(lambda),1);
Rs = zeros(length(lambda),1);
Ts = zeros(length(lambda),1);

n1_IR = interp1(SiO2(:,1),SiO2(:,2),1064,'pchip');
n2_IR = interp1(Ta2O5(:,1),Ta2O5(:,2),1064,'pchip');

for i = 1:length(lambda)
    n1 = interp1(SiO2(:,1),SiO2(:,2),lambda(i)*1064,'pchip');
    n2 = interp1(Ta2O5(:,1),Ta2O5(:,2),lambda(i)*1064,'pchip');
    nb = n1;
    if strcmp(firstLayer,'SiO2')
        n_c = [];
        L_temp = [];
        for kk = 1:(no_of_stacks)
            n_c  = [n_c n1 n2];
            L_temp(2*kk - 1) = L(2*kk - 1)*n1/n1_IR;
            L_temp(2*kk) = L(2*kk)*n2/n2_IR;
        end
        n = [na n_c nb];  % add the indices of the vacuum and the substrate at either end 
    [Gammap, ~] = multidiel1(n, L_temp, lambda(i), aoi, 'tm');
    [Gammas, ~] = multidiel1(n, L_temp, lambda(i), aoi, 'te');
    Rp(i) = abs(Gammap).^2;
    Tp(i) = 1 - Rp(i);
    Rs(i) = abs(Gammas).^2;
    Ts(i) = 1 - Rs(i);
    
    elseif strcmp(firstLayer, 'Ta2O5') %AR coating
        n_c = [];
        L_temp = [];
        for kk = 1:(no_of_stacks)
            n_c  = [n_c n2 n1];
            L_temp(2*kk - 1) = L(2*kk - 1)*n2/n2_IR;
            L_temp(2*kk) = L(2*kk)*n1/n1_IR;
        end
        n = [nb n_c na];  % add the indices of the vacuum and the substrate at either end 
        aoi_temp = theta2(n1,n2,41.1);
        [Gammap, ~] = multidiel1(n, L_temp, lambda(i), aoi_temp, 'tm');
        [Gammas, ~] = multidiel1(n, L_temp, lambda(i), aoi_temp, 'te');
        Rp(i) = abs(Gammap).^2;
        Tp(i) = 1 - Rp(i);
        Rs(i) = abs(Gammas).^2;
        Ts(i) = 1 - Rs(i);
    end
end
        

lambda_real = lambda * lambda_0 * 1e9;
figure(70711)
hold on
subplot(2,1,1)
h1=semilogy(lambda_real, Rp, 'b',...
lambda_real, Tp, 'r',...
'LineWidth', 4);
xlabel('Wavelength [nm]')
ylabel('R or T, p-pol')
legend('R_{p-pol}','T_{p-pol}', 'Location', 'SouthEast')
grid
axis([min(lambda_real) max(lambda_real) .9e-6 1.01])
line(lambda_0*1e9*ones(100,1), logspace(-7,0,100), 'Color','c','LineWidth',3)
text(1100, 3e-6, ['T @ ' num2str(lambda_0*1e9) ' = ' num2str(TNout.Tp_IR*1e2,3) '%'],...
'FontSize', 26)
% nice print
set( gca                       , ...
'FontName'   , 'Times'     );
% set([hXLabel, hYLabel], ...
%     'FontName'   , 'Times',...
%     'FontSize'   , 34         );
% set([hLegend, gca]             , ...
%     'FontSize'   , 22           );
%set( hTitle                    , ...
%    'FontSize'   , 12          , ...
%    'FontWeight' , 'bold'      );
set(gca, ...
'Box'         , 'on'     , ...
'TickDir'     , 'in'     , ...
'TickLength'  , [.02 .02] , ...
'XMinorTick'  , 'on'      , ...
'YMinorTick'  , 'on'      , ...
'YGrid'       , 'on'      , ...
'XColor'      , .1*[.3 .3 .3], ...
'YColor'      , .1*[.3 .3 .3], ...
'YTick'       , logspace(-6, 0, 7), ...
'FontSize'    , 26 ,...
'LineWidth'   , 2);
line(lambda_0*1e9/2*ones(100,1),logspace(-7,0,100),'Color','g','LineWidth',3)
text(532, 5e-6,['R @  ' num2str(lambda_0*1e9/2) ' nm = ' num2str(TNout.Rp_green*1e2) '%'],...
'FontSize', 26)
subplot(2,1,2)
h2=semilogy(lambda_real, Rs, 'b',...
lambda_real, Ts, 'r',...
'LineWidth', 4);
xlabel('Wavelength [nm]')
ylabel('R or T, s-pol')
legend('R_{s-pol}','T_{s-pol}', 'Location', 'SouthEast')
grid
axis([min(lambda_real) max(lambda_real) .9e-6 1.01])
%title('Calculated AlGaAs Mirror Reflectivity')
line(lambda_0*1e9/2*ones(100,1),logspace(-7,0,100),'Color','g','LineWidth',3)
line(lambda_0*1e9*ones(100,1), logspace(-7,0,100), 'Color','c','LineWidth',3)
%text(1600, 1e-3, ['R @ ' num2str(lambda_0*1e9) ' = ' num2str(R1(1),4)])
text(1100, 3e-6, ['T @ ' num2str(lambda_0*1e9) ' = ' num2str(TNout.Ts_IR*1e2,3) ' %'],...
'FontSize', 26)
%text(1300,0.1*1.3^-1,['R @  ' num2str(lambda_0*1e9/2) ' = ' num2str(R1(1))])
text(532, 5e-6,['R @  ' num2str(lambda_0*1e9/2) ' nm = ' num2str(TNout.Rs_green*1e2), '%'],...
'FontSize', 26)
% nice print
set( gca                       , ...
'FontName'   , 'Times'     );
% set([hXLabel, hYLabel], ...
%     'FontName'   , 'Times',...
%     'FontSize'   , 34         );
% set([hLegend, gca]             , ...
%     'FontSize'   , 22           );
%set( hTitle                    , ...
%    'FontSize'   , 12          , ...
%    'FontWeight' , 'bold'      );
set(gca, ...
'Box'         , 'on'     , ...
'TickDir'     , 'in'     , ...
'TickLength'  , [.02 .02] , ...
'XMinorTick'  , 'on'      , ...
'YMinorTick'  , 'on'      , ...
'YGrid'       , 'on'      , ...
'XColor'      , .1*[.3 .3 .3], ...
'YColor'      , .1*[.3 .3 .3], ...
'YTick'       , logspace(-6, 0, 7), ...
'FontSize'    , 26 ,...
'LineWidth'   , 2);
