function [lobe_freq,lobe_mag,lobe_phase] = lobe_reconst(ACR,FCR,LUT_ACR,LUT_FCR,x_lobe,mag_lobe,arg_lobe,N_lobe)
%% Projet long 2015
% Lobe reconstruction

% Arguments
% ACR = Amplitude Change Rate
% FCR = Frequency Change Rate

% LUT_ACR   = ACR of the reference points
% LUT_FCR   = FCR of the reference points
% x_lobe    = Frequency of the reference points
% mag_lobe  = Magnitude of the reference points
% arg_lobe  = Argument of the reference points

lobe_freq = zeros(N_lobe,1);
lobe_mag = zeros(N_lobe,1);
lobe_phase = zeros(N_lobe,1);

% Perform bilinear interpolation to find the new points of the lobe
for v = 1:N_lobe
    lobe_freq(v,1) = griddata(LUT_ACR,LUT_FCR,x_lobe(v,:),ACR,FCR);
    lobe_mag(v,1) = griddata(LUT_ACR,LUT_FCR,mag_lobe(v,:),ACR,FCR);
    lobe_phase(v,1) = griddata(LUT_ACR,LUT_FCR,arg_lobe(v,:),ACR,FCR);
end    
    
end