function y = getMirrorCost_AR(x, params, flag)
% OPTETM This function calculates something about the coating
% it gets used by the do<Optic>.m program as the thing to
% minimize. The input argument 'x' is the initial guess for the layer structure
%
% Option to compute R/T for 532nm at multiple polarizations...
% This macro changes the order of the layers to account for designing an AR
% coating, putting it in the way multidiel1.m wants it...

%Importing the parameter called "NUMTOOLS" in do<Optic>.m
glob_param = params;

% This makes x into a column vector if it isn't already
x = x(:);

% WTF?
try, flag; catch flag = 0; end

ifo = glob_param.ifo;

% Laser wavelength
lambda_0 = ifo.Laser.Wavelength;

% Angle of incidence
aoi = glob_param.aoi;

na = 1.000;    % Index of vacuum
n1 = ifo.Materials.Coating.Indexlown;             % Index of SiO2  @ 1064 nm
n2 = ifo.Materials.Coating.Indexhighn;            % Index of Ta2O5 @ 1064 nm
nb = ifo.Materials.Substrate.RefractiveIndex;     % Substrate is made of SiO2

no_of_stacks = floor(length(x)/2);    % use floor if x is not even
% set up the (alternating) array of index of refractions
n_c = [];

for kk = 1:(no_of_stacks)
  n_c  = [n_c n2 n1];
end
       

n = [nb n_c na];  % add the indices of the vacuum and the substrate at either end

%Making this a row vector for multidiel1.m
L = transpose(x);

% list of wavelengths to use for the optimization
lambda = [lambda_0 lambda_0/2] / lambda_0;
%lambda = [lambda_0] / lambda_0;

% Calculate reflectivities for each polarization at the angle of incidence
% ['tm' = 'p-pol', 'te' = 's-pol'], for the given layer structure
[Gammap, ~]   = multidiel1(n, L, lambda, aoi, 'tm');
[Gammas, ~]   = multidiel1(n, L, lambda, aoi, 'te');

% Note the arrays Gammap/Gammas contains information at 1064nm as well as 532nm
% for use to compute error fcn...

Rp = abs(Gammap).^2;
Tp = 1 - Rp;

Rs = abs(Gammas).^2;
Ts = 1 - Rs;

% define the error function
T_1 = glob_param.T_1;         % desired  transmission @ lambda_0
R_1 = 1 - T_1;

T_2s = glob_param.T_2s;         % desired transmission @ lambda_0 / 2
R_2s = 1 - T_2s;      % currently not used

T_2p = glob_param.T_2p;         % desired transmission @ lambda_0 / 2
R_2p = 1 - T_2p;      % currently not used

% ----Brownian Thermal noise - parameters from GWINC IFOmodel file
% Formula from E0900068-v3, see also T0900161
phi_high = ifo.Materials.Coating.Phihighn;
n_high   = n2; 
Y_high   = ifo.Materials.Coating.Yhighn;

phi_low  = ifo.Materials.Coating.Philown; 
n_low    = n1;
Y_low    = ifo.Materials.Coating.Ylown; 

Y_sub  = ifo.Materials.Substrate.MirrorY;

% Total thickness
z_low  = sum(L(1:2:length(L)));  %odd layers
z_high = sum(L(2:2:length(L)));  %even layers

a = phi_high / phi_low;
b = n_low / n_high;
c = Y_high/Y_sub + Y_sub/Y_high;
d = Y_low/Y_sub  + Y_sub/Y_low;
little_gamma = a*b*c/d;

% Brownian thermal noise estimate (only used for optimization)
S = z_low + (little_gamma * z_high);
%This is a kind of proxy function for minimizing coating thermal noise
% ---------------------------------------------------


% Thermo-Optic Noise (there should be a flag to use this or not)
f_to  = glob_param.f_optimize;     %Frequency to evaluate this
wBeam = glob_param.wBeam;
dOpt  = L';
%[StoZ, SteZ, StrZ, T]  = getCoatThermoOptic(f_to, ifo, wBeam, dOpt);

% cost function which gets minimized - equal penalties given at 532nm and
% 1064nm
yy = [];
yy(1) = .0001*S;                           % minimize the Brownian noise
     
yy = [yy 111*abs((Tp(1) - T_1)/T_1)^1];    % match the T @ lambda

yy = [yy 111*abs((Tp(2) - T_2p)/T_2p)^1];  % match the T @ lambda/2, p-pol

yy = [yy 111*abs((Ts(2) - T_2s)/T_2s)^1];  % match the T @ lambda/2, s-pol

% minimize reflected E-field (usually by 1/2 wave cap on top)
% only at main wavelength
r_refl = abs(1 + Gammap(1));
yy     = [yy 10*r_refl^2];

if flag==2
    yy
elseif flag==3
    yy
    keyboard    % debug
end

% choose which terms in the cost function to use
% 1 = Brownian Noise
% 2 = Transmission at 1064, p-pol
% 3 = Transmission at 532, p-pol
% 4 = Transmission at 532, s-pol
% 5 = HR surface field
y = sum(yy([3 4]));  %For the AR, we only care about 3 and 4...

%This is the error function (evaluated) we are trying to minimize...
sss = y;

%Plot option...
if flag == 1
 yy
 clear y
 y.n = n;
 y.L = L;
 y.lambda = lambda;
 y.Tp = Tp;
 y.Ts = Ts;
 y.Rp = Rp;
 y.Rs = Rs;
 y.sss = sss;
 y.Sthermal = S;
 y.yy = yy;
 
 lambda = sort([linspace(0.4, 1.6, 2200), lambda]);

 [Gammap, ~] = multidiel1(n, L, lambda, aoi, 'tm');
 [Gammas, ~] = multidiel1(n, L, lambda, aoi, 'te');

 R = abs(Gammap).^2;
 T = 1 - R;
 
 R2 = abs(Gammas).^2;
 T2 = 1 - R2;

 lambda_real = lambda * lambda_0 * 1e9;
 
figure(70711)
hold on
subplot(2,1,1)
semilogy(lambda_real, R, 'b',...
         lambda_real, T, 'r',...
         'LineWidth', 4)
xlabel('Wavelength [nm]')
ylabel('R or T, p-pol')
legend('R','T', 'Location', 'SouthEast')
grid
axis([min(lambda_real) max(lambda_real) .9e-6 1.01])
%title('Calculated AlGaAs Mirror Reflectivity')

line(lambda_0*1e9/2*ones(100,1),logspace(-7,0,100),'Color','g','LineWidth',3)

% line(lambda_0*1e9*ones(100,1), logspace(-7,0,100), 'Color','c','LineWidth',3)
%text(1600, 1e-3, ['R @ ' num2str(lambda_0*1e9) ' = ' num2str(R1(1),4)])
% text(1100, 3e-6, ['T @ ' num2str(lambda_0*1e9) ' = ' num2str(Tp(1)*1e6,3) ' ppm'],...
%     'FontSize', 26)

%text(1300,0.1*1.3^-1,['R @  ' num2str(lambda_0*1e9/2) ' = ' num2str(R1(1))])
text(532, 5e-6,['R @  ' num2str(lambda_0*1e9/2) ' = ' num2str(Rp(2)*1e6,3) 'ppm'],...
    'FontSize', 26)
%S
% disp(['T @ ' num2str(lambda_0*1e9) ' nm, p-pol = ' num2str(Tp(1)*1e6,5) ' ppm'])
disp(['R @ ' num2str(lambda_0*1e9/2) ' nm, p-pol = ' num2str(Rp(2)*1e6,5) ' ppm'])


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
%%%%%%%%%%%%%%%  s-pol   %%%%%%%%%%%%%%%%%%%%
subplot(2,1,2)
semilogy(lambda_real, R2, 'b',...
         lambda_real, T2, 'r',...
         'LineWidth', 4)
xlabel('Wavelength [nm]')
ylabel('R or T, s-pol')
legend('R','T', 'Location', 'SouthEast')
grid
axis([min(lambda_real) max(lambda_real) .9e-6 1.01])
%title('Calculated AlGaAs Mirror Reflectivity')

line(lambda_0*1e9/2*ones(100,1),logspace(-7,0,100),'Color','g','LineWidth',3)

% line(lambda_0*1e9*ones(100,1), logspace(-7,0,100), 'Color','c','LineWidth',3)
%text(1600, 1e-3, ['R @ ' num2str(lambda_0*1e9) ' = ' num2str(R1(1),4)])
% text(1100, 3e-6, ['R @ ' num2str(lambda_0*1e9) ' = ' num2str(Rs(1)*1e6,3) ' ppm'],...
%     'FontSize', 26)

%text(1300,0.1*1.3^-1,['R @  ' num2str(lambda_0*1e9/2) ' = ' num2str(R1(1))])
text(532, 5e-6,['R @  ' num2str(lambda_0*1e9/2) ' = ' num2str(Rs(2)*1e6,3) 'ppm'],...
    'FontSize', 26)
%S
% disp(['T @ ' num2str(lambda_0*1e9) ' nm, s-pol = ' num2str(Ts(1)*1e6,5) ' ppm'])
disp(['R @ ' num2str(lambda_0*1e9/2) ' nm, s-pol = ' num2str(Rs(2)*1e6,5) ' ppm'])


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
end

