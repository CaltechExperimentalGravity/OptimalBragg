function y = getMirrorCost(x, params, flag)
% OPTETM This function calculates something about the ETM coating
% it gets used by the do_aSi.m program as the thing to
% minimize. The input argument 'x' is the initial guess for the layer structure
%

glob_param = params;

% This makes x into a column vector
x = x(:);

% WTF?
try, flag; catch flag = 0; end

ifo = glob_param.ifo;

% Laser wavelength
lambda_0 = ifo.Laser.Wavelength;

na = 1.000;    % Index of vacuum
n1 = ifo.Materials.Coating.Indexlown;             % Index of SiO2  @ 1064 nm
n2 = ifo.Materials.Coating.Indexhighn;            % Index of Ta2O5 @ 1064 nm
nb = ifo.Materials.Substrate.RefractiveIndex;     % Substrate is made of SiO2

no_of_stacks = floor(length(x)/2);    % use floor if x is not even
% set up the array of index of refractions
n_c = [];
n_c1 = [];
n_c2 = [];
for kk = 1:(no_of_stacks)
  n_c  = [n_c n1 n2];
  n_c1 = [n_c1 n1*1.01 n2];
  n_c2 = [n_c2 n1 n2*1.01];  
end
       

n = [na n_c nb];  % add the indices of the vacuum and the substrate

% Comment from Stefan here =>>
%L((length(L)-length(x))+1:length(L)) = transpose(x);
L = transpose(x);

% list of wavelengths to use for the optimization
lambda = [lambda_0 lambda_0/2] / lambda_0;
%lambda = [lambda_0] / lambda_0;

% Calculate reflectivities
[Gamma1, ~]   = multidiel1(n, L, lambda);
%[Gamma1e, Z1e] = multidiel1(n, L, lambda*1.01);

% (dT/T) / (dlambda / lambda)
%dTdlambda = (abs(Gamma1).^2 - abs(Gamma1e).^2)/0.01;
%dTdlambda = dTdlambda ./ (1 - abs(Gamma1).^2);

% (dT/T)/(dL/L)
% dT/T < 10% for a 1% change in L
%[Gamma5, Z5] = multidiel1(n, 1.01*L, lambda);
%dTdL = (abs(Gamma1).^2 - abs(Gamma5).^2)/0.01;
%dTdL = dTdL ./ (1 - abs(Gamma1).^2);

% (dT/T)/(dn1/n1)
% dT/T < 10% for a 1% change in n1
%[Gamma7, Z7] = multidiel1([na n_c1 nb], L, lambda);
%dTdn1 = (abs(Gamma7).^2 - abs(Gamma1).^2)/0.01;
%dTdn1 = dTdn1 ./ (1 - abs(Gamma1).^2);

% (dT/T)/(dn/n)
% dT/T < 10% for a 1% change in n2
%[Gamma9, Z9] = multidiel1([na n_c2 nb], L, lambda);
%dTdn2 = (abs(Gamma9).^2 - abs(Gamma1).^2)/0.01;
%dTdn2 = dTdn2 ./ (1 - abs(Gamma1).^2);

% reflection and transmission at main wavelength
R1 = abs(Gamma1).^2;
T1 = 1 - R1;

% define the error function
T_1 = glob_param.T_1;         % desired  transmission @ lambda_0
R_1 = 1 - T_1;

T_2 = glob_param.T_2;         % desired transmission @ lambda_0 / 2
R_2 = 1 - T_2;      % currently not used

% ----Brownian Thermal noise - parameters from GWINC IFOmodel file
% Formula from E0900068-v3
phi_high = ifo.Materials.Coating.Phihighn;
n_high   = n2; 
Y_high   = ifo.Materials.Coating.Yhighn;

phi_low  = ifo.Materials.Coating.Philown; 
n_low    = n1;
Y_low    = ifo.Materials.Coating.Ylown; 

Y_sub  = ifo.Materials.Substrate.MirrorY;

% Total thickness
z_low  = sum(L(1:2:length(L)));
z_high = sum(L(2:2:length(L)));

a = phi_high / phi_low;
b = n_low / n_high;
c = Y_high/Y_sub + Y_sub/Y_high;
d = Y_low/Y_sub  + Y_sub/Y_low;
little_gamma = a*b*c/d;

% Brownian thermal noise estimate (only used for optimization)
S = z_low + (little_gamma * z_high);
% ---------------------------------------------------


% Thermo-Optic Noise (there should be a flag to use this or not)
f_to  = glob_param.f_optimize;
wBeam = glob_param.wBeam;
dOpt  = L';
%[StoZ, SteZ, StrZ, T]  = getCoatThermoOptic(f_to, ifo, wBeam, dOpt);

% cost function which gets minimized
yy = [];
yy(1) = .0001*S;                       % minimize the Brownian noise
     
yy = [yy 111*abs((T1(1) - T_1)/T_1)^1];    % match the T @ lambda

yy = [yy 111*abs((T1(2) - T_2)/T_2)^1];    % match the T @ lambda/2
% Thermo-Optic noise cancellation
% the weight factor is chosen so that StoZ < Sbrown/10
%y = y + StoZ * 1e46;

% also add some costs to minimize the first derivatives
% minimize sensitivity to thickness cal of deposition
%yy = [yy 1*(abs(dTdL(1)))];

%yy = [yy 5*abs(dTdlambda(1))^2];

% minimize sensitivity to index of refraction mis-calibration
%yy = [yy 1*(abs(dTdn1(1)))];  % dn1 = 1% 
%yy = [yy 1*(abs(dTdn2(1)))];  % dn2 = 1%

% minimize reflected E-field (usually by 1/2 wave cap on top)
% only at main wavelength
r_refl = abs(1 + Gamma1(1));
yy     = [yy 10*r_refl^2];

if flag==2
    yy
elseif flag==3
    yy
    keyboard    % debug
end

% choose which terms in the cost function to use
% 1 = Brownian Noise
% 2 = Transmission at 1064
% 3 = Transmission at 532
% 4 = HR surface field
y = sum(yy([1 2 3 4]));

sss = y;

if flag == 1
 yy
 clear y
 y.n = n;
 y.L = L;
 y.lambda = lambda;
 y.T1 = T1;
 y.R1 = R1;
 y.sss = sss;
 y.Sthermal = S;
 y.yy = yy;
 
 lambda = sort([linspace(0.4, 1.6, 2200), lambda]);

 [Gamma1, ~] = multidiel1(n, L, lambda);

 R = abs(Gamma1).^2;
 T = 1 - R;

 lambda_real = lambda * lambda_0 * 1e9;
 
figure(70711)
semilogy(lambda_real, R, 'b',...
         lambda_real, T, 'r',...
         'LineWidth', 4)
xlabel('Wavelength [nm]')
ylabel('Reflectivity or Transmissivity')
legend('R','T', 'Location', 'SouthEast')
grid
axis([min(lambda_real) max(lambda_real) .9e-6 1.01])
%title('Calculated AlGaAs Mirror Reflectivity')

line(lambda_0*1e9/2*ones(100,1),logspace(-7,0,100),'Color','g')

line(lambda_0*1e9*ones(100,1), logspace(-7,0,100), 'Color','c')
%text(1600, 1e-3, ['R @ ' num2str(lambda_0*1e9) ' = ' num2str(R1(1),4)])
text(1100, 3e-6, ['T @ ' num2str(lambda_0*1e9) ' = ' num2str(T1(1)*1e6,3) ' ppm'],...
    'FontSize', 26)

%text(1300,0.1*1.3^-1,['R @  ' num2str(lambda_0*1e9/2) ' = ' num2str(R1(1))])
text(532, 5e-6,['T @  ' num2str(lambda_0*1e9/2) ' = ' num2str(T1(2))],...
    'FontSize', 26)
%S
disp(['T @ ' num2str(lambda_0*1e9) ' nm = ' num2str(T1(1)*1e6,5) ' ppm'])

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

