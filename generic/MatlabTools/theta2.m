%function that calculates theta2 using Snell's law...
%For use by the coating optimization tool, takes into account the wedge angle
%as well...

function t2 = theta2(n1, n2, theta1)
t2 = asind(n1*sind(theta1)/n2);
t2 = t2 - 2;
