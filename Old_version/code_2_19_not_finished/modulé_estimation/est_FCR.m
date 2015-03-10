function y=est_FCR(x0,x1,x2,y0,y1,y2,p0,p1,p2,w0_tild,p)
    y=p*(derive_2(x0,x1,x2,p0,p1,p2,w0_tild))/derive_2(x0,x1,x2,y0,y1,y2,w0_tild);
end