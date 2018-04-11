function [x,z] = parload(fname)
file = load(fname);

try
    x = file.x; z = file.z; 
catch
    x = file.x; z = [];
end

end
