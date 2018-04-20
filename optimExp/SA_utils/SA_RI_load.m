function [RI] = SA_RI_load(RI_index)
RI_file = load('RI_set.mat');
RI = RI_file.RI_set(RI_index);
end

