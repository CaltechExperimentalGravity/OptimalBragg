function [Temp] = SA_Temp_load(Initial_temp_index)
Temp_file = load('Temp_set.mat');
Temp = Temp_file.Temp_set(Initial_temp_index);
end

