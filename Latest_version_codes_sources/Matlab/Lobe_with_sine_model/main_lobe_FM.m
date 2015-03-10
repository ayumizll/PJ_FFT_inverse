function y = main_lobe_FM(x,FCR)
% Projet long 2015
% Quentin Biache - Lili Zheng - Emilie Abia

% This function reconstructs a lobe using a sinusoidal approximation. 
% The constants generated in the file FM_coefficients.mat were generated 
% using the function sine_interp.m in the same folder.

% x = normalized frequency
% FCR = frequency change rate (in Hz/s)

% Load constants
load FM_coefficients.mat

% Third order polynomial fitting of the constants versus the FCR
u1 = a1.p1*abs(FCR).^3 + a1.p2*FCR.^2 + a1.p3*abs(FCR) + a1.p4;
u2 = a2.p1*abs(FCR).^3 + a2.p2*FCR.^2 + a2.p3*abs(FCR) + a2.p4;

v1 = b1.p1*abs(FCR).^3 + b1.p2*FCR.^2 + b1.p3*abs(FCR) + b1.p4;
v2 = b2.p1*abs(FCR).^3 + b2.p2*FCR.^2 + b2.p3*abs(FCR) + b2.p4;

w1 = c1.p1*FCR + c1.p2;
w2 = c2.p1*FCR + c2.p2;

% Build the lobe
y = u1.*sin(v1.*x+w1)+u2.*sin(v2.*x+w2);


end

