addpath('../');
% load Data/PR3_AR_layers_170411_1140.mat
% load Data/PR3_AR_layers_170410_1321.mat
load Data/PR3_20_layers_170416_1409.mat
% load Data/PR3_layers_170411_1823.mat

dOpt = TNout.L(:);
n = TNout.n_IR(:);
dReal = dOpt ./ n(2:end-1);
figure(20413)
clf
subplot('Position', [0.085 0.5 0.85 0.4])
colrs = {'r','k'};
hold on
for k=1:length(dOpt)
    bar(k, dOpt(k), 0.95, colrs{rem(k,2)+1})
end
hold off
%xlabel('Layer #')
ylabel('Optical Thickness / \lambda')
title('SiO_2/Ta_2O_5 Coating')
%axis tight
legend('SiO_2')
ylim([0 .71])
xlim([0 length(dReal)+1])
set(gca,'XTick',[])
set(gca,'YTick',[0:0.1:1])
axis tight
%line(1:length(dOpt2),0.25*ones(size(dOpt2)),'Color','y')
grid

subplot('Position', [0.085 0.085 0.85 0.4])
hold on
for k=1:length(dReal)
    bar(k, dReal(k), 0.95, colrs{rem(k,2)+1})
end
hold off
xlabel('Layer #')
ylabel('Physical Thickness / \lambda')
axis tight
set(gca,'XTick',[1:length(dReal)])
xlim([0 length(dReal)+1]);
%legend('SiO2:a-Si 123 K');
grid
