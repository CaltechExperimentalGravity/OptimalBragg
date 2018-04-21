% PLOTTO
%
% plot thermo-optic noise for some mirror
%

% try
%     precompIFO;
% catch me
%     addpath(genpath('../'));
%     addpath(genpath('../../gwincDev'));
% end

f = logspace(0, 4, 300);

%%
% get time str for latest coating design
d = dir('Data/*.mat');
[~,idx] = max([d.datenum]);
filename = d(idx).name;
strfind(filename,'_18');
tnowstr = filename(12:22)
if strfind(d(idx).name,'ETM')
    NUMTOOLS.opt_name = 'ETM';
else
    NUMTOOLS.opt_name = 'ITM';
end

fname = ['Data/' NUMTOOLS.opt_name '_layers_' tnowstr];

load(fname);

dOpt  = TNout.L(:);
n     = TNout.n(:);
%divide out by index of refraction (omitting vacuum and substrate)
dReal = dOpt ./ n(2:end-1);
%
wBeam = eval(['ifo.Optics.' NUMTOOLS.opt_name '.BeamRadius']);
wfac  = 0.06/wBeam;

[StoZ, SteZ, StrZ, T]  = getCoatThermoOptic(f, ifo, wBeam, dOpt);

% Coating Brownian
SbrZ = getCoatBrownian(f, ifo, wBeam, dOpt);

% Mirror Brownian
c2  = ifo.Materials.Substrate.c2;
yy  = ifo.Materials.Substrate.MechanicalLossExponent;
kBT = ifo.Constants.kB * ifo.Materials.Substrate.Temp;

% Bulk substrate contribution
phibulk = c2 .* (f.^yy);

[cETM, aETM] = subbrownianFiniteCorr(ifo, 'ETM');
cbulk = 8 * kBT * aETM .* phibulk ./ (2 * pi * f);

% sub brown
sbr = sqrt(cbulk');

to  = sqrt(StoZ); % thermo optic
te  = sqrt(SteZ); % thermo elastic
tr  = sqrt(StrZ); % thermo refractive
tbr = sqrt(SbrZ); % coating brownian

%%
figure(415)
loglog(f, sbr,...
       f, tbr,...
       f,  te,...
       f,  tr,...
       f,  to,...
       'LineWidth', 6)

xlabel('Frequency [Hz]')
ylabel('Dispacement Noise [m/\surdHz]')
title('Single Mirror Thermo-Optic Noise (AlGaAs coating)')
grid
%axis tight
axis([1 1e4 1e-23*wfac 1e-19*wfac])
text(23, 1.3e-20*wfac, ['T = ' num2str(T*1e6,4) ' ppm'])
text(23, 5.1e-20*wfac, ['# of layers = ' num2str(length(dReal)) ' '])
text(23, 2.5e-20*wfac, ['Thickness = ' num2str(sum(dReal), 3) ' \mum'])
legend('Substrate Brownian','Brownian',...
    'Thermo-Elastic','Thermo-Refractive','Thermo-Optic',...
       'Location','NorthEast')

   
% nice print
set( gca                       , ...
    'FontName'   , 'Times'     , ...
    'FontSize'   , 28          );

set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , .1*[.3 .3 .3], ...
  'YColor'      , .1*[.3 .3 .3], ...
  'FontSize'    , 28 ,...
  'LineWidth', 1);
% print pretty plot
orient landscape
set(gcf,'Position',[200 80 1000 700])
set(gcf,'PaperPositionMode','auto')
fname = ['Figures/' NUMTOOLS.opt_name '_AlGaAs_TOnoise_' tnowstr];
print('-depsc','-r300',fname)
[a,~] = system(['Figures/makePDF.sh ' fname '.eps']);
if a ~= 0
    disp('PDF Generation Error')
end

%% plot coating design

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
title('GaAs:AlAs Coating')
%axis tight
legend('GaAs')
ylim([0 .71])
xlim([0 length(dReal)+1])
set(gca,'XTick',[])
set(gca,'YTick',[0:0.1:1])
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
set(gca,'XTick',[1:2:length(dReal)])
xlim([0 length(dReal)+1]);
%legend('SiO2:a-Si 123 K');
grid

% print pretty plot
orient landscape
set(gcf,'Position',[500 250 1000 700])
set(gcf,'PaperPositionMode','auto')
%print -depsc -r600 AlGaAs_Layers_60000.eps
fname = ['Figures/' NUMTOOLS.opt_name '_AlGaAs_Layers_' tnowstr];
print('-depsc','-r300', fname)
[a,~] = system(['Figures/makePDF.sh ' fname '.eps']);
if a ~= 0
    disp('PDF Generation Error')
else
    system('rm Figures/*.eps');
end

