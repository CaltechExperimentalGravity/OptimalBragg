 function parsave(fname,x,z)

 if nargin == 3
     save(fname, 'x','z')
 elseif nargin == 2
     save(fname, 'x')
 end

 end

