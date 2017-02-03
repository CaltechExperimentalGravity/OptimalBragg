function y = optETM2(x,flag)

% OPTETM This function calculates something about the ETM coating
% it gets used by the conjGradSolve.m program as the thing to
% minimize. The input argument 'x' is the initial guess for the layers.
%
% P.S. This is called a comment field and I encourage people to
% start thinking about using these habitually.


% This makes x into a column vector
x = x(:);

% WTF?
try, flag; catch flag = 0; end

% Laser wavelength
lambda_0 = 1064e-9;

na = 1.000;   % Index of vacuum
n1 = 1.450;   % Index of SiO2 @ 1064 nm    - both n1 & n2 from Bench 6.2
n2 = 2.065;   % Index of Ta2O5 @ 1064 nm
nb = n1;      % Substrate is made of SiO2

no_of_stacks = length(x)/2;
if rem(length(x),2)
    error('Odd number of layers')
end

n = [];
for kk = 1:(no_of_stacks)
  n = [n n1 n2];
end

n = [na n nb];

% Comment from Stefan here =>>
%L((length(L)-length(x))+1:length(L)) = transpose(x);
L = transpose(x);


lambda = [473,532,946,1064,1319,1550] / 1064;

[Gamma1, Z1] = multidiel1(n, L, lambda);

R1 = abs(Gamma1).^2;
T1 = 1-R1;

% define the error function
T_1064 = 10e-6;         % desired ETM transmission @ 1064
R_1064 = 1 - T_1064;

T_532 = 0.10;           % desired ETM transmission @ 532
R_532 = 1 - T_532;


y = 1 * (abs(R1(4)) - R_1064)^2 +...
    1 * (abs(R1(2)) - R_532)^2;


if flag
 clear y;
 y.n = n;
 y.L = L;
 y.lambda = lambda;
 y.T1 = T1;
 y.R1 = R1;
 lambda = sort([linspace(0.29,1.5,1000),lambda]);

 [Gamma1, Z1] = multidiel1(n, L, lambda);

 R = abs(Gamma1).^2;
 T = 1 - R;

 lambda_real = lambda * lambda_0 * 1e9;

semilogy(lambda_real, R, 'b',...
         lambda_real, T, 'r')
xlabel('Wavelength (nm)')
ylabel('Reflectivity or Transmissivity')
legend('R','T')
grid
axis([min(lambda_real) max(lambda_real) 1e-2 1])
title('Calculated ETM HR Reflectivity')
%line(473*ones(100,1),linspace(0,1,100),'Color','c')
line(532*ones(100,1),linspace(0,1,100),'Color','g')
%line(946*ones(100,1),linspace(0,1,100),'Color','m')
line(1064*ones(100,1),linspace(0,1,100),'Color','y')
text(1300,0.1*1.3^3,['R @ 1064 = ' num2str(R1(4))])
text(1300,0.1*1.3^2,['T @ 1064 = ' num2str(T1(4))])
%text(1500,0.2*1.3^1,['R @  946 = ' num2str(R1(3))])
%text(1500,0.2*1.3^0,['T @  946 = ' num2str(T1(3))])
text(1300,0.1*1.3^-1,['R @  532 = ' num2str(R1(2))])
text(1300,0.1*1.3^-2,['T @  532 = ' num2str(T1(2))])
%text(1500,0.2*1.3^-3,['R @  473 = ' num2str(R1(1))])
%text(1500,0.2*1.3^-4,['T @  473 = ' num2str(T1(1))])
%print -dpng AdvLIGO_etm_hr.png

end

