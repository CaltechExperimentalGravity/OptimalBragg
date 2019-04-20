%Function to convert optical thickness to physical thickness (in units of reference wavelength)
function phys = op2phys(L,n)
phys = L./n;
end
