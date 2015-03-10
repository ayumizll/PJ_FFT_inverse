function y = quad_max(u,v,w)
% Long project 2015
% Quentin Biache
% This function returns the maximum of 3 points fitted by a
% parabola. 
% Arguments

% u = f[p-1]
% v = f[p]
% w = f[p+1]
% p is a given integer (its value doesn't matter)

y = (u.^2) + (16.0*(v.^2)) + (w.^2);
y = y + (-8*u.*v) + (-8*v.*w) + (-2*u.*w);
y = y./(u-(2*v)+w);
y = y*(-1./8.);
end

