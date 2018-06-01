%Function to perturb a given coating design to see how sensitive 
%derived parameters like Transmission etc depend on perturbations
%Since perturbations are assumed to be iid, no fancy MCMC sampler is used
%If we desire, we can adapt this code to use 
%https://www.mathworks.com/matlabcentral/fileexchange/49820-ensemble-mcmc-sampler
function doMC_AlGaAs(filename, N, nDim, savename)
clc
close all

%Load the coating file...
load(filename);

%number of output variables
nVars = 4;

%Generate the perturbations...
means = zeros(nDim,1);
sigmas = diag(diag((0.005^2)*ones(nDim),0));
perturbs = mvnrnd(means,sigmas,N);

%Define arrays for output variables of interest
T_IR = zeros(N,1);
surfField = zeros(N,1);
TOnoise = zeros(N,1);
BRnoise = zeros(N,1);

%Nominal values from coating design
aoi = 0;
n = TNout.n;
L = TNout.L;
L_phys = op2phys(L,n(2:end-1));

%Nominal params for calculation
Ei = 27.46; %[V/m], for surface field calculation
f_to = 100; % [Hz]
wBeam = 0.065; % [meters]

%%%%%   Apply the perturbations   %%%%%%%%
for i=1:N
    n_IRs = n;
    Ls = L_phys;
    perturb = 1 + perturbs(i,:);
    %All physical thicknesses with equal error
    Ls = Ls * perturb(1);   
    %Error in refractive index of low index layers
    n_IRs(2:2:end-1) = n_IRs(2:2:end-1) .* perturb(2);
    %Error in refractive index of high index layers
    n_IRs(3:2:end-1) = n_IRs(3:2:end-1) .* perturb(3);

    %Compute reflectivity, surface field, Brownian noise and TO noise.
    [Gamma, ~] = multidiel1(n_IRs, Ls.*n_IRs(2:end-1), TNout.lambda);
    T_IR(i) =  1e6*(1 - abs(Gamma).^2);
    %Thermo-Optic 
    ifo.Materials.Coating.Indices = n_IRs;
    [StoZ, ~, ~, Tto]  = getCoatThermoOptic(f_to, ifo, wBeam, Ls'.*n_IRs(2:end-1)');
    TOnoise(i) = sqrt(StoZ);
    %Brownian
    SbrZ = getCoatBrownian(f_to, ifo, wBeam, Ls.*n_IRs(2:end-1));
    BRnoise(i) = sqrt(SbrZ);
    %Surface Field
    surfField(i) = Ei * abs(1+Gamma);
end

%Save everything for corner plotting with python
h5create(savename,'/MCout',[nVars N]);
h5write(savename,'/MCout',horzcat(T_IR, 1e21*TOnoise, 1e21*BRnoise, surfField)');

end