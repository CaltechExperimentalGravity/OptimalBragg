function y=optPRM(x,flag)
x=x(:);
try, flag; catch flag=0; end

lambda_0 = 1064e-9;

no_of_stacks = 6;  % How many pairs of high/low 1/4 wave stacks

na = 1.000; % Index of vacuum
n1 = 1.45231;  % Index of SiO2 @ 1064 nm
n2 = 2.02936;  % Index of Ta2O5 @ 1064 nm
nb = n1;    % Substrate is made of SiO2

n = [n1];
L = [1/2];
for kk = 1:no_of_stacks
  n = [n n2 n1];
  %L = [L 1/4 1/4];
  L = [L 1/8 3/8];
end

n = [na n n2 n1 n2 n1 n2 nb];
%L = [L 1/4 0.2961 0.2145 0.2231 0.3519];
L = [L 1/8 3/8 1/8 3/8 1/8];

L((length(L)-length(x))+1:length(L))=transpose(x);

lambda = [0.5,1,946/1064,1319/1064,1550/1064];

[Gamma1, Z1] = multidiel1(n, L, lambda);

R1 = abs(Gamma1).^2;
T1 = 1-R1;

% define the error function
y=10*(R1(2)-1+0.03)^2 + 1*(R1(3)-0.00)^2 + 0*sum(x(2:2:end)-1/8)^2;

%y=10*(R1(2)-0.9860)^2  -R1(1);
%y=10*(R1(2)-1+1e-5)^2 + 1*(R1(1)-0.8)^2;


if flag
clear y;
y.n=n;
y.L=L;
y.lambda=lambda;
y.T1=T1;
y.R1=R1;
lambda = sort([linspace(0.29,2.1,4000),lambda]);

[Gamma1, Z1] = multidiel1(n, L, lambda);

R = abs(Gamma1).^2;
T = 1-R;

lambda_real = lambda * lambda_0 * 1e9;

semilogy(lambda_real, R, 'b',...
         lambda_real, T, 'r')
xlabel('Wavelength (nm)')
ylabel('Reflectivity or Transmissivity')
legend('R','T')
grid
axis([min(lambda_real) max(lambda_real) 1e-2 1])
title('Calculated PRM HR Reflectivity')
line(532*ones(100,1),linspace(0,1,100),'Color','g')
line(946*ones(100,1),linspace(0,1,100),'Color','g')
line(1064*ones(100,1),linspace(0,1,100),'Color','y')
text(1500,0.1*1.3^3,['R @ 1064 = ' num2str(R1(2))])
text(1500,0.1*1.3^2,['T @ 1064 = ' num2str(T1(2))])
text(1500,0.1*1.3^1,['R @  946 = ' num2str(R1(3))])
text(1500,0.1*1.3^0,['T @  946 = ' num2str(T1(3))])
text(1500,0.1*1.3^-1,['R @  532 = ' num2str(R1(1))])
text(1500,0.1*1.3^-2,['T @  532 = ' num2str(T1(1))])
print -dpng AdvLIGO_prm_hr.png

end

