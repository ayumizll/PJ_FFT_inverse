function [a,b,c] = quad_coeffs(x1,x2,x3,y1,y2,y3)

denom1 = (x3-x1).*(x2-x1);
denom2 = (x3-x2).*(x2-x1);
denom3 = (x3-x1).*(x3-x2);

a = (y1./denom1) + (-y2./denom2) + (y3./denom3);
b = (-y1.*(x2+x3)./denom1) + (y2.*(x1+x3)./denom2) + (-y3.*(x1+x2)./denom3);
c = (y1.*x2.*x3./denom1) + (-y2.*x1.*x3./denom2) + (y3.*x1.*x2./denom3);

end

