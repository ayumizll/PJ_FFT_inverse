function out = derive_2(x0,x1,x2,y0,y1,y2,x)
    out = -(2*(x0.*y1 - x1.*y0 - x0.*y2 + x2.*y0 + x1.*y2 - x2.*y1))./((x0 - x1).*(x0 - x2).*(x1 - x2));
end