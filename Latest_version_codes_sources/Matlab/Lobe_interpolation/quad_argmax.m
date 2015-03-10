function y = quad_argmax(index,u,v,w)
y = u-w;
y = y./(u-(2*v)+w);
y = 0.5*y;
y = y + index;
end

   
