%%Run some kind of MC to figure out how params will be distributed...
%Number of MC samples...
clc
clear
close all

N = 5000;

load Data/PR3_layers_170412_1530.mat
savename = 'PR3_HR_MC_170412_1530.dat'
% load Data/PR3_AR_layers_170411_1140.mat

aoi = 41.1;
mu_aoi = 0.;%Conservative estimate on error on angle of incidence...
sig_aoi = 0.01;

AOIs = aoi.*(1+normrnd(mu_aoi,sig_aoi,N,1));


%Evaluate derivatives for varying each of these params...
T_IR_P_aoi = zeros(N,1);
R_green_P_aoi = zeros(N,1);
R_green_S_aoi = zeros(N,1);

T_IR_P_n1 = zeros(N,1);
R_green_P_n1 = zeros(N,1);
R_green_S_n1 = zeros(N,1);

T_IR_P_n2 = zeros(N,1);
R_green_P_n2 = zeros(N,1);
R_green_S_n2 = zeros(N,1);

T_IR_P_SiO2 = zeros(N,1);
R_green_P_SiO2 = zeros(N,1);
R_green_S_SiO2 = zeros(N,1);

T_IR_P_L1 = zeros(N,1);
R_green_P_L1 = zeros(N,1);
R_green_S_L1 = zeros(N,1);

T_IR_P_L2 = zeros(N,1);
R_green_P_L2 = zeros(N,1);
R_green_S_L2 = zeros(N,1);

T_IR_P_combined = zeros(N,1);
R_green_P_combined = zeros(N,1);
R_green_S_combined = zeros(N,1);

%First, perturb individually...
%%%%%   AOI   %%%%%%%%
for i=1:N
    n_IRs = TNout.n_IR;
    n_greens = TNout.n_green;
    Ls = TNout.L;
    [Gamma5p, ~] = multidiel1(n_IRs, Ls, TNout.lambda, AOIs(i), 'tm');
    T_IR_P_aoi(i) =  1 - abs(Gamma5p).^2;
    [Gamma5p, ~] = multidiel1(n_greens, 2*Ls, 0.5*TNout.lambda, AOIs(i), 'tm');
    R_green_P_aoi(i) =  abs(Gamma5p).^2;
    [Gamma5s, ~] = multidiel1(n_greens, 2*Ls, 0.5*TNout.lambda, AOIs(i), 'te');
    R_green_S_aoi(i) =  abs(Gamma5s).^2; 
end

%%%%%   n1   %%%%%%%%
for i=1:N
    n_IRs = TNout.n_IR;
    n_greens = TNout.n_green;
    Ls = TNout.L;
    error = normrnd(0.,0.01,1);
    for jj=2:length(n_IRs)-1
        if mod(jj+1,2) %i.e. is a SiO2 layer
            n_IRs(jj) = n_IRs(jj)*(1+error);
            n_greens(jj) = n_greens(jj)*(1+error);
        end
    end
    [Gamma5p, ~] = multidiel1(n_IRs, Ls, TNout.lambda, aoi, 'tm');
    T_IR_P_n1(i) =  1 - abs(Gamma5p).^2;
    [Gamma5p, ~] = multidiel1(n_greens, 2*Ls, 0.5*TNout.lambda, aoi, 'tm');
    R_green_P_n1(i) =  abs(Gamma5p).^2;
    [Gamma5s, ~] = multidiel1(n_greens, 2*Ls, 0.5*TNout.lambda, aoi, 'te');
    R_green_S_n1(i) =  abs(Gamma5s).^2; 
end 

%%%%%   n2   %%%%%%%%
for i=1:N
    n_IRs = TNout.n_IR;
    n_greens = TNout.n_green;
    Ls = TNout.L;
    error = normrnd(0.,0.01,1);
    for jj=2:length(n_IRs)-1
        if mod(jj,2) %i.e. is a Ta2O5 layer
            n_IRs(jj) = n_IRs(jj)*(1+error);
            n_greens(jj) = n_greens(jj)*(1+error);
        end
    end
    [Gamma5p, ~] = multidiel1(n_IRs, Ls, TNout.lambda, aoi, 'tm');
    T_IR_P_n2(i) =  1 - abs(Gamma5p).^2;
    [Gamma5p, ~] = multidiel1(n_greens, 2*Ls, 0.5*TNout.lambda, aoi, 'tm');
    R_green_P_n2(i) =  abs(Gamma5p).^2;
    [Gamma5s, ~] = multidiel1(n_greens, 2*Ls, 0.5*TNout.lambda, aoi, 'te');
    R_green_S_n2(i) =  abs(Gamma5s).^2; 
end 

%%%%%   SiO2   %%%%%%%%
for i=1:N
    n_IRs = TNout.n_IR;
    n_greens = TNout.n_green;
    Ls = TNout.L;
    error = normrnd(0.,0.01,1);
    n_IRs(end) = n_IRs(end)*(1+error);
    n_greens(end) = n_greens(end)*(1+error);
    [Gamma5p, ~] = multidiel1(n_IRs, Ls, TNout.lambda, aoi, 'tm');
    T_IR_P_SiO2(i) =  1 - abs(Gamma5p).^2;
    [Gamma5p, ~] = multidiel1(n_greens, 2*Ls, 0.5*TNout.lambda, aoi, 'tm');
    R_green_P_SiO2(i) =  abs(Gamma5p).^2;
    [Gamma5s, ~] = multidiel1(n_greens, 2*Ls, 0.5*TNout.lambda, aoi, 'te');
    R_green_S_SiO2(i) =  abs(Gamma5s).^2; 
end 

%%%%%   Thickness of SiO2 layers   %%%%%%%%
for i=1:N
    n_IRs = TNout.n_IR;
    n_greens = TNout.n_green;
    Ls = TNout.L;
    error = normrnd(0.,0.01,1);
    for jj=1:length(Ls)
        if mod(jj,2) %is even i.e. is a SiO2 layer
            Ls(jj) = Ls(jj)*(1+error);
        end
    end
    [Gamma5p, ~] = multidiel1(n_IRs, Ls, TNout.lambda, aoi, 'tm');
    T_IR_P_L1(i) =  1 - abs(Gamma5p).^2;
    [Gamma5p, ~] = multidiel1(n_greens, 2*Ls, 0.5*TNout.lambda, aoi, 'tm');
    R_green_P_L1(i) =  abs(Gamma5p).^2;
    [Gamma5s, ~] = multidiel1(n_greens, 2*Ls, 0.5*TNout.lambda, aoi, 'te');
    R_green_S_L1(i) =  abs(Gamma5s).^2; 
end

%%%%%   Thickness of Ta2O5 layers   %%%%%%%%
for i=1:N
    n_IRs = TNout.n_IR;
    n_greens = TNout.n_green;
    Ls = TNout.L;
    error = normrnd(0.,0.01,1);
    for jj=1:length(Ls)
        if mod(jj+1,2) %is odd, i.e. is a Ta2O5 layer
            Ls(jj) = Ls(jj)*(1+error);
        end
    end
    [Gamma5p, ~] = multidiel1(n_IRs, Ls, TNout.lambda, aoi, 'tm');
    T_IR_P_L2(i) =  1 - abs(Gamma5p).^2;
    [Gamma5p, ~] = multidiel1(n_greens, 2*Ls, 0.5*TNout.lambda, aoi, 'tm');
    R_green_P_L2(i) =  abs(Gamma5p).^2;
    [Gamma5s, ~] = multidiel1(n_greens, 2*Ls, 0.5*TNout.lambda, aoi, 'te');
    R_green_S_L2(i) =  abs(Gamma5s).^2; 
end

%%%%%   Combine everything   %%%%%%%%
for i=1:N
    n_IRs = TNout.n_IR;
    n_greens = TNout.n_green;
    Ls = TNout.L;
    e1 = normrnd(0.,0.03,1);
    e2 = normrnd(0.,0.01,1);
    e3 = normrnd(0.,0.01,1);
    e4 = normrnd(0.,0.01,1);
    e5 = normrnd(0.,0.01,1);
    e6 = normrnd(0.,0.01,1);
    
    aoi_temp = aoi*(1+e1);
    for jj=2:length(n_IRs)-1
        if mod(jj+1,2) %i.e. is a SiO2 layer
            n_IRs(jj) = n_IRs(jj)*(1+e2);
            n_greens(jj) = n_greens(jj)*(1+e2);
        end
    end
    for jj=2:length(n_IRs)-1
        if mod(jj,2) %i.e. is a Ta2O5 layer
            n_IRs(jj) = n_IRs(jj)*(1+e3);
            n_greens(jj) = n_greens(jj)*(1+e3);
        end
    end
    n_IRs(end) = n_IRs(end)*(1+e4);
    n_greens(end) = n_greens(end)*(1+e4);
    for jj=1:length(Ls)
        if mod(jj,2) %is even i.e. is a SiO2 layer
            Ls(jj) = Ls(jj)*(1+e5);
        end
    end    
    for jj=1:length(Ls)
        if mod(jj+1,2) %is odd, i.e. is a Ta2O5 layer
            Ls(jj) = Ls(jj)*(1+e5);
        end
    end
    [Gamma5p, ~] = multidiel1(n_IRs, Ls, TNout.lambda, aoi_temp, 'tm');
    T_IR_P_combined(i) =  1 - abs(Gamma5p).^2;
    [Gamma5p, ~] = multidiel1(n_greens, 2*Ls, 0.5*TNout.lambda, aoi_temp, 'tm');
    R_green_P_combined(i) =  abs(Gamma5p).^2;
    [Gamma5s, ~] = multidiel1(n_greens, 2*Ls, 0.5*TNout.lambda, aoi_temp, 'te');
    R_green_S_combined(i) =  abs(Gamma5s).^2; 
end

%Save everything for pretty plotting with python
data = [T_IR_P_aoi R_green_P_aoi R_green_S_aoi ...
    T_IR_P_n1 R_green_P_n1 R_green_S_n1 ...
    T_IR_P_n2 R_green_P_n2 R_green_S_n2 ...
    T_IR_P_SiO2 R_green_P_SiO2 R_green_S_SiO2 ...
    T_IR_P_L1 R_green_P_L1 R_green_S_L1 ...
    T_IR_P_L2 R_green_P_L2 R_green_S_L2 ...
    T_IR_P_combined R_green_P_combined R_green_S_combined ...
    ];

csvwrite(savename,data)

% %%%%%%%%%%%%%%%%%%%
% %PLOT everything %%
% %%%%%%%%%%%%%%%%%%%
% figure(539)
% hold on
% subplot(2,2,1)
% histogram(T_IR_P_aoi*1e6,50,'Normalization','cdf');
% xlabel('T @1064nm, p-pol [%]','FontSize',16,'FontWeight','bold')
% ylabel('CDF','FontSize',16,'FontWeight','bold')
% title('Varying AOI','FontSize',16,'FontWeight','bold');
% set(gca,'FontSize',16,'FontWeight','bold')
% subplot(2,2,2)
% histogram(R_green_P_aoi*1e2,50,'Normalization','cdf');
% xlabel('R @532nm, p-pol [%]','FontSize',16,'FontWeight','bold')
% ylabel('CDF','FontSize',16,'FontWeight','bold')
% set(gca,'FontSize',16,'FontWeight','bold')
% subplot(2,2,3)
% histogram(R_green_S_aoi*1e2,50,'Normalization','cdf');
% xlabel('R @532nm, s-pol [%]','FontSize',16,'FontWeight','bold')
% ylabel('CDF','FontSize',16,'FontWeight','bold')
% set(gca,'FontSize',16,'FontWeight','bold')
% 
% figure(8972)
% hold on
% subplot(2,2,1)
% histogram(T_IR_P_n1*1e6,50,'Normalization','cdf');
% xlabel('T @1064nm, p-pol [%]','FontSize',16,'FontWeight','bold')
% ylabel('CDF','FontSize',16,'FontWeight','bold')
% title('Varying n1','FontSize',16,'FontWeight','bold');
% set(gca,'FontSize',16,'FontWeight','bold')
% subplot(2,2,2)
% histogram(R_green_P_n1*1e2,50,'Normalization','cdf');
% xlabel('R @532nm, p-pol [%]','FontSize',16,'FontWeight','bold')
% ylabel('CDF','FontSize',16,'FontWeight','bold')
% set(gca,'FontSize',16,'FontWeight','bold')
% subplot(2,2,3)
% histogram(R_green_S_n1*1e2,50,'Normalization','cdf');
% xlabel('R @532nm, s-pol [%]','FontSize',16,'FontWeight','bold')
% ylabel('CDF','FontSize',16,'FontWeight','bold')
% set(gca,'FontSize',16,'FontWeight','bold')
% 
% figure(5309)
% hold on
% subplot(2,2,1)
% histogram(T_IR_P_n2*1e6,50,'Normalization','cdf');
% xlabel('T @1064nm, p-pol [%]','FontSize',16,'FontWeight','bold')
% ylabel('CDF','FontSize',16,'FontWeight','bold')
% title('Varying n2','FontSize',16,'FontWeight','bold');
% set(gca,'FontSize',16,'FontWeight','bold')
% subplot(2,2,2)
% histogram(R_green_P_n2*1e2,50,'Normalization','cdf');
% xlabel('R @532nm, p-pol [%]','FontSize',16,'FontWeight','bold')
% ylabel('CDF','FontSize',16,'FontWeight','bold')
% set(gca,'FontSize',16,'FontWeight','bold')
% subplot(2,2,3)
% histogram(R_green_S_n2*1e2,50,'Normalization','cdf');
% xlabel('R @532nm, s-pol [%]','FontSize',16,'FontWeight','bold')
% ylabel('CDF','FontSize',16,'FontWeight','bold')
% set(gca,'FontSize',16,'FontWeight','bold')
% 
% figure(4882)
% hold on
% subplot(2,2,1)
% histogram(T_IR_P_SiO2*1e6,50,'Normalization','cdf');
% xlabel('T @1064nm, p-pol [%]','FontSize',16,'FontWeight','bold')
% ylabel('CDF','FontSize',16,'FontWeight','bold')
% title('Varying nSiO2','FontSize',16,'FontWeight','bold');
% set(gca,'FontSize',16,'FontWeight','bold')
% subplot(2,2,2)
% histogram(R_green_P_SiO2*1e2,50,'Normalization','cdf');
% xlabel('R @532nm, p-pol [%]','FontSize',16,'FontWeight','bold')
% ylabel('CDF','FontSize',16,'FontWeight','bold')
% set(gca,'FontSize',16,'FontWeight','bold')
% subplot(2,2,3)
% histogram(R_green_S_SiO2*1e2,50,'Normalization','cdf');
% xlabel('R @532nm, s-pol [%]','FontSize',16,'FontWeight','bold')
% ylabel('CDF','FontSize',16,'FontWeight','bold')
% set(gca,'FontSize',16,'FontWeight','bold')
% 
% figure(5532)
% hold on
% subplot(2,2,1)
% histogram(T_IR_P_L1*1e6,50,'Normalization','cdf');
% xlabel('T @1064nm, p-pol [%]','FontSize',16,'FontWeight','bold')
% ylabel('CDF','FontSize',16,'FontWeight','bold')
% title('Varying L1','FontSize',16,'FontWeight','bold');
% set(gca,'FontSize',16,'FontWeight','bold')
% subplot(2,2,2)
% histogram(R_green_P_L1*1e2,50,'Normalization','cdf');
% xlabel('R @532nm, p-pol [%]','FontSize',16,'FontWeight','bold')
% ylabel('CDF','FontSize',16,'FontWeight','bold')
% set(gca,'FontSize',16,'FontWeight','bold')
% subplot(2,2,3)
% histogram(R_green_S_L1*1e2,50,'Normalization','cdf');
% xlabel('R @532nm, s-pol [%]','FontSize',16,'FontWeight','bold')
% ylabel('CDF','FontSize',16,'FontWeight','bold')
% set(gca,'FontSize',16,'FontWeight','bold')
% 
% 
% figure(9478)
% hold on
% subplot(2,2,1)
% histogram(T_IR_P_L2*1e6,50,'Normalization','cdf');
% xlabel('T @1064nm, p-pol [%]','FontSize',16,'FontWeight','bold')
% ylabel('CDF','FontSize',16,'FontWeight','bold')
% title('Varying L2','FontSize',16,'FontWeight','bold');
% set(gca,'FontSize',16,'FontWeight','bold')
% subplot(2,2,2)
% histogram(R_green_P_L2*1e2,50,'Normalization','cdf');
% xlabel('R @532nm, p-pol [%]','FontSize',16,'FontWeight','bold')
% ylabel('CDF','FontSize',16,'FontWeight','bold')
% set(gca,'FontSize',16,'FontWeight','bold')
% subplot(2,2,3)
% histogram(R_green_S_L2*1e2,50,'Normalization','cdf');
% xlabel('R @532nm, s-pol [%]','FontSize',16,'FontWeight','bold')
% ylabel('CDF','FontSize',16,'FontWeight','bold')
% set(gca,'FontSize',16,'FontWeight','bold')
% 
% figure(956)
% hold on
% subplot(2,2,1)
% histogram(T_IR_P_combined*1e6,50,'Normalization','cdf');
% xlabel('T @1064nm, p-pol [%]','FontSize',16,'FontWeight','bold')
% ylabel('CDF','FontSize',16,'FontWeight','bold')
% title('Varying everything','FontSize',16,'FontWeight','bold');
% set(gca,'FontSize',16,'FontWeight','bold')
% subplot(2,2,2)
% histogram(R_green_P_combined*1e2,50,'Normalization','cdf');
% xlabel('R @532nm, p-pol [%]','FontSize',16,'FontWeight','bold')
% ylabel('CDF','FontSize',16,'FontWeight','bold')
% set(gca,'FontSize',16,'FontWeight','bold')
% subplot(2,2,3)
% histogram(R_green_S_combined*1e2,50,'Normalization','cdf');
% xlabel('R @532nm, s-pol [%]','FontSize',16,'FontWeight','bold')
% ylabel('CDF','FontSize',16,'FontWeight','bold')
% set(gca,'FontSize',16,'FontWeight','bold')