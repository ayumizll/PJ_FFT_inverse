function y = quad_argmax(index,u,v,w)
% Long project 2015
% Quentin Biache
% This function returns the position of the maximum of 3 points fitted by a
% parabola. 
% Arguments

% index = position of the maximum (integer)
% u = f[index-1]
% v = f[index]
% w = f[index+1]


y = u-w;
y = y./(u-(2*v)+w);
y = 0.5*y;
y = y + index;
end

   
