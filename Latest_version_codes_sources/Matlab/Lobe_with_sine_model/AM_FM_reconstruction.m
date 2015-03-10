%% Projet long 2015
% Quentin Biache - Lili Zheng - Emilie Abia

% Lobe with different AM/FM modulation animation
% This script plot different parameters of the main lobe and their
% evolution as the FCR (or the ACR) changes.

% In this version, the signal is reconstructed now with a sinusoidal model 
% for the main lobe.

% In this case, only the magnitude of the lobe was reconstructed since the
% research were lead for the magnitude only.

%% Initialize workspace
clear all
close all
clc

%% Parameters definition

% Sampling parameters
fe = 44100;             % Sampling frequency

% FFT parameters
w_time = 0.023;         % Window duration (s)
N_padding = 7;          % Zero-padding factor
N_pts = round(fe*w_time);               % Total number of points
if (mod(N_pts,2) == 0)
    N_pts = N_pts + 1;
end

N_fft = 2^(nextpow2(N_pts)+N_padding);  % FFT size
t = (0:(N_pts-1))'/fe;                  % Time vector
f = (0:(N_fft-1))'*fe/N_fft;            % Frequency vector

w = window(@hanning,N_pts);             % Create window
w = w./sum(w);                          % Normalize window

% Put a '1' to apply FM
AM = 0;

% Put a '1' here to plot magnitudes in dB ('0' for linear scale)
dB = 1;

%% Compute FFT

% Number of simulations
N = 51;

if AM
    alpha0_dB = linspace(-0.1,-12,N)/0.023;
    FCR = 0;
else
    alpha0_dB = 0;
    FCR = linspace(-10000,10000,N);
end    
    
% Matrix to store lobe data
levels = zeros(9,N);
phase_levels = zeros(9,N);
levels_normalized = zeros(9,N);
levels_loc = zeros(9,N);

figure
pause

for u = 1:N

    % Parameters for the model
    f0 = 1000;              % Instantaneous frequency at t = 0
    A0 = 1;                 % Signal level
    phi0 = 0;               % Initial phase

    % Adjust parameters with the definition
    if AM        
        alpha0 = log(10^(alpha0_dB(u)/20));
        beta0 = FCR/pi;
    else
        alpha0 = log(10^(alpha0_dB/20));
        beta0 = FCR(u)/pi;
    end    
        
    omega0 = 2*pi*f0;
    lambda0 = log(A0);

    % Create signal
    s = exp(alpha0.*t).*exp(lambda0).*exp(1i*((beta0*t.^2)+(omega0.*t)+phi0));
 
    % Apply window
    s = s.*w;    
    
    % Plot the signal
    subplot(2,2,1)
        plot(t,real(s))
        title('Real part of the signal')
        xlabel('Time (s)')
        ylabel('Signal')
        grid on

    % Window phase correction        
    s_w = zeros(N_fft,1);
    s_w(1:(N_pts-1)/2) = s(((N_pts+1)/2)+1:end);
    s_w(N_fft-(N_pts-1)/2:end) = s(1:((N_pts+1)/2));
    
    % Compute the FFT
    N_half = (N_fft/2)+1;                   % Half FFT size
    fft_s = fft(s_w,N_fft);
    mod_fft_s = abs(fft_s);

    % Define the plot region (centered around f0)
    plot_index = find((f > (f0-200)) & (f < (f0+200)));
    
    subplot(2,2,3)
        if dB
            plot(f(plot_index),20*log(mod_fft_s(plot_index)))
        else
            plot(f(plot_index),mod_fft_s(plot_index))
        end
        title('|X_{DFT}(f)|^2')
        xlabel('Frequency (Hz)')
        ylabel('|X_{DFT}(f)|^2')
        grid on
    
    % Search for the main lobe
    [mod_max,index] = max(mod_fft_s);
    mod_argmax = f(index);

    % Find a more accurate value of the main lobe peak
    interp_index = quad_argmax(index,mod_fft_s(index-1),mod_fft_s(index),mod_fft_s(index+1));
    interp_max = quad_max(mod_fft_s(index-1),mod_fft_s(index),mod_fft_s(index+1));
    interp_argmax = interp1(1:length(f),f,interp_index);

    left_peak = mod_fft_s(1:end-2) > mod_fft_s(2:end-1);
    right_peak = mod_fft_s(3:end) > mod_fft_s(2:end-1);

    % Find location of zeros in the FFT
    mod_zeros_loc = find(left_peak & right_peak) + 1; 

    % Find the closest zero to the lobe peak
    [~,index_min] = min(abs(index-mod_zeros_loc)); 

    if (mod_zeros_loc(index_min) > index)
        lower_zero_loc = mod_zeros_loc(index_min - 1); 
        upper_zero_loc = mod_zeros_loc(index_min); 
    else
        lower_zero_loc = mod_zeros_loc(index_min); 
        upper_zero_loc = mod_zeros_loc(index_min + 1); 
    end    

    % Create vectors to plot the reconstructed main lobe
    x_lobe_reconst = (lower_zero_loc:upper_zero_loc)';
    f_lobe_reconst = (x_lobe_reconst-1)*fe/N_fft;
    y_lobe_reconst = main_lobe_FM((f_lobe_reconst-f0)/fe,FCR(u));
    
    % Analysis to synthesis index conversion factor
    C = (N_pts-1)/(N_fft-1);
    
    % Frequency indexes of the lower and upper bound of the main lobe in
    % the synthesis domain
    a = round(1+(C*(lower_zero_loc-1)));
    b = round(1+(C*(upper_zero_loc-1)));
    
    vec_ab = linspace(C*(lower_zero_loc-interp_index-1),C*(upper_zero_loc-interp_index-1),b-a+1)';
    
    p_index = (lower_zero_loc:upper_zero_loc)';
        
    % Phase in the synthesis domain (less FFT points)
    phase_reconst = interp1(1+(C*(p_index-1)),unwrap(phase(fft_s(p_index))),(1:N_pts)',[],0);
    
    % Reconstruct signal
    Y = zeros(N_pts,1);
    Y(a:b) = main_lobe_FM(vec_ab/N_pts,FCR(u)).*exp(1i*phase_reconst(a:b));
    s_reconst = real(ifft(Y,N_pts));

    % Undo window phase correction
    s_reconst_w = zeros(N_pts,1);
    s_reconst_w(1:(N_pts+1)/2) = s_reconst(((N_pts+1)/2):end);
    s_reconst_w(1+((N_pts+1)/2):end) = s_reconst(1:((N_pts-1)/2));
    
    s_reconst = s_reconst_w;
      
    % Plot reconstructed signal
    subplot(2,2,1)
        hold on
        plot(t,s_reconst,'r')
        hold off      
    
    % Plot the FFT and reconstructed FFT (synthesis)
    subplot(2,2,2)
        Y_in = fft(s,N_pts);
        f_in = (0:(N_pts-1))'*fe/N_pts;
        plot(f_in,20*log10(abs(Y_in)),'b',f_in,20*log10(abs(Y)),'r')
        grid on
        legend('Original','Reconstructed')
    
    % Plot the reconstructed lobe (analysis)
    subplot(2,2,3)
        hold on
        if dB
            plot(f_lobe_reconst,20*log(abs(y_lobe_reconst)),'r') 
        else
            plot(f_lobe_reconst,y_lobe_reconst,'r')
        end
        hold off
                
    % Plot the phase
    subplot(2,2,4)
        plot(f(plot_index),unwrap(phase(fft_s(plot_index))))
        title('arg(X_{DFT}(f))')
        xlabel('Frequency (Hz)')
        ylabel('arg(X_{DFT}(f))')
        grid on        
        
    pause(0.1)
end


subplot(2,2,2)
    title('Reconstructed FFT magnitude (synthesis domain)')
    