function out = derive_1(x0,x1,x2,y0,y1,y2,x)
    out = (x0.^2.*y1 - x1.^2.*y0 - x0.^2.*y2 + x2.^2.*y0 + x1.^2.*y2 - x2.^2.*y1)./((x0 - x1).*(x0 - x2).*(x1 - x2)) - (2*x.*(x0*y1 - x1.*y0 - x0.*y2 + x2.*y0 + x1.*y2 - x2.*y1))./((x0 - x1).*(x0 - x2).*(x1 - x2));
end