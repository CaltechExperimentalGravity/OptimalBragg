%Run multiple calls to runSwarm_aLIGO_ETM.m
%for different values of the loss angle of tantala
clear
clc

ifo = SilicaTantala300; 
phi_Ta = [1., 1.5, 2.0, 2.5, 3.0];
phi_Ta = 1./phi_Ta;
phi_Ta = phi_Ta .* ifo.Materials.Coating.Phihighn;

for i=1:length(phi_Ta)
    loss_highn=phi_Ta(i);
    runSwarm_aLIGO_ETM;
end

