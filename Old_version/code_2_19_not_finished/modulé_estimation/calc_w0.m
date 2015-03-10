function w0 = calc_w0(x1,y0,y1,y2)
    w0 = x1+0.5*(y0-y2)./(y0-2*y1+y2);
end