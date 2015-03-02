%% Projet long 2015
% IEEE algorithm for AM/FM parameters estimation


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
N_fft = 2^(nextpow2(N_pts)+N_padding);  % FFT size
t = (0:(N_pts-1))'/fe;                  % Time vector
f = (0:(N_fft-1))'*fe/N_fft;            % Frequency vector

w = window(@hanning,N_pts);             % Create window
w = w./sum(w);                          % Normalize window

AM = 1;

%% Compute FFT

% Number of simulations
N = 85;

if AM
    alpha0_dB = linspace(-0.1,-20,N)/0.023;
    FCR = 0;
else
    alpha0_dB = 0;
    FCR = linspace(-1000,1000,N);
end    
    
% Matrix to store lobe data
levels = zeros(9,N);
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

    subplot(2,2,1)
        plot(t,real(s))
        title('Real part of the signal')
        xlabel('Time (s)')
        ylabel('Signal')
        grid on

    % Apply window
    s = s.*w;

    % Compute the FFT
    N_half = (N_fft/2)+1;                   % Half FFT size
    fft_s = fft(s,N_fft);
    mod_fft_s = abs(fft_s);
    arg_fft_s = 0;%unwrap(phase(fft_s));

    plot_index = find((f > 800) & (f < 1200));
    
    subplot(2,2,3)
        plot(f(plot_index),20*log(mod_fft_s(plot_index)))
        title('|X_{DFT}(f)|^2')
        xlabel('Frequency (Hz)')
        ylabel('|X_{DFT}(f)|^2')
        grid on
        %hold on
    
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

    x = linspace(lower_zero_loc,upper_zero_loc,9)';

    levels_normalized(:,u) = interp1((1:N_fft)',mod_fft_s,x)/interp_max;
    levels(:,u) = interp1((1:N_fft)',mod_fft_s,x);
    
    subplot(2,2,3)
        hold on
        plot(mod_argmax,20*log(mod_max),'r+')
        plot(interp_argmax,20*log(interp_max),'g+')
        plot((x-1)*fe/N_fft,20*log(levels(:,u)),'r+')
        hold off
        
    levels_loc(:,u) = (x-1)*fe/N_fft;
        
    subplot(2,2,2)
        if AM        
            plot(alpha0_dB,20*log10(levels_normalized(1:5,:)'))
            title('Amplitude of the point versus the amplitude modulation rate')
            xlabel('Amplitude modulation rate (dB/s)')
            ylabel('Relative point level (dB)')
        else
            plot(FCR,20*log10(levels_normalized(1:5,:))')
            title('Relative amplitude versus the frequency modulation rate')
            xlabel('Frequency modulation rate (Hz/s)')
            ylabel('Relative point level (dB)')
        end 
        %legend('1','2','3','4','5')
        grid on        

    subplot(2,2,4)
        if AM        
            plot(alpha0_dB,levels_loc(1:5,:)'-f0)
            xlabel('Amplitude modulation rate (dB/s)')
            ylabel('Relative location (Hz)')
        else
            plot(FCR,levels_loc(1:5,:)'-f0)
            xlabel('Amplitude modulation rate (dB/s)')
            ylabel('Relative location (Hz)')
        end 
        %legend('1','2','3','4','5')
        grid on        
        
    disp([num2str(round(100*u/N)),' % achieved'])
    pause(0.01)
end

subplot(2,2,2)
    legend('1','2','3','4','5')

subplot(2,2,4)
    legend('1','2','3','4','5')
