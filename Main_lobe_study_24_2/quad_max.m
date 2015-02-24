function y = quad_max(u,v,w)
y = (u.^2) + (16.0*(v.^2)) + (w.^2);
y = y + (-8*u.*v) + (-8*v.*w) + (-2*u.*w);
y = y./(u-(2*v)+w);
y = y*(-1./8.);
end

