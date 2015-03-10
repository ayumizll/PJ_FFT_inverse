%% Projet long 2015
% Quentin Biache - Lili Zheng - Emilie Abia

% Lobe with different AM/FM modulation animation
% This script plot different parameters of the main lobe and their
% evolution as the FCR (or the ACR) changes.

% The main lobe is sampled with 9 points (it can be adjusted with N_lobe)
% At the end, the script returns a .mat file containing the main lobe
% sampled coordinates. This script allows also to study either FM or AM,
% and to plot the values in dB or in linear scale.

% The main lobe points are stored for further studies (like sine_interp.m)
%% Initialize workspace
clear all
close all
clc

%% Parameters definition

% Sampling parameters
fe = 44100;             % Sampling frequency

% FFT parameters
w_time = 0.0115;         % Window duration (s)
N_padding = 8;          % Zero-padding factor
N_pts = round(fe*w_time);               % Total number of points
if (mod(N_pts,2) == 0)
    N_pts = N_pts + 1;
end

N_fft = 2^(nextpow2(N_pts)+N_padding);  % FFT size
t = (0:(N_pts-1))'/fe;                  % Time vector
f = (0:(N_fft-1))'*fe/N_fft;            % Frequency vector

w = window(@hanning,N_pts);             % Create window
w = w./sum(w);                          % Normalize window

AM = 0;
dB = 1;

%% Compute FFT

% Number of simulations
N = 31;

% Number of points on the main lobe
N_lobe = 9;

if AM
    alpha0_dB = linspace(+12,-12,N)/0.023;
    FCR = 0;
else
    alpha0_dB = 0;
    FCR = linspace(-30000,30000,N);
end    
    
% Matrix to store lobe data
levels = zeros(N_lobe,N);
levels_normalized = zeros(N_lobe,N);
levels_loc = zeros(N_lobe,N);

x_lobe = zeros(N_lobe,N);
y_lobe = zeros(N_lobe,N);

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
        
    s_w = zeros(N_fft,1);
    s_w(1:(N_pts-1)/2) = s(((N_pts+1)/2)+1:end);
    s_w(N_fft-(N_pts-1)/2:end) = s(1:((N_pts+1)/2));
    
    
    % Correct phase
    
    
    % Compute the FFT
    N_half = (N_fft/2)+1;                   % Half FFT size
    fft_s = fft(s_w,N_fft);
    mod_fft_s = abs(fft_s);
    arg_fft_s = 0;%unwrap(phase(fft_s));

    plot_index = find((f > 800) & (f < 1200));
    %plot_index = find((f > 964) & (f < 966));
    
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

    x = linspace(lower_zero_loc,upper_zero_loc,N_lobe)';

    levels_normalized(:,u) = interp1((1:N_fft)',mod_fft_s,x)/interp_max;
    levels(:,u) = interp1((1:N_fft)',mod_fft_s,x);
    
    x_lobe_interp = linspace(lower_zero_loc,upper_zero_loc,N_pts)';
    f_lobe_interp = (x_lobe_interp-1)*fe/N_fft;
    y_lobe_interp = quad_interp(x,levels(:,u),x_lobe_interp);
    
    % Reconstruct signal
    Y = zeros(N_pts,1);     % Complex spectrum
    Y = quad_interp(x*N_pts/N_fft,levels(:,u),(1:N_pts)')*phase(fft_s(index));
    s_reconst = real(ifft(Y,N_pts));

    s_reconst_w = zeros(N_pts,1);
    s_reconst_w(1:(N_pts+1)/2) = s_reconst(((N_pts+1)/2):end);
    s_reconst_w(1+((N_pts+1)/2):end) = s_reconst(1:((N_pts-1)/2));
    
    s_reconst = s_reconst_w;
    
%     subplot(2,2,1)
%         hold on
%         plot(t,s_reconst,'r')
%         hold off    
    
    subplot(2,2,3)
        hold on
        if dB
            plot(mod_argmax,20*log(mod_max),'r+')
            plot(interp_argmax,20*log(interp_max),'g+')
            plot((x-1)*fe/N_fft,20*log(levels(:,u)),'r+')
            plot(f_lobe_interp,20*log(y_lobe_interp),'r') 
        else
            plot(mod_argmax,mod_max,'r+')
            plot(interp_argmax,interp_max,'g+')
            plot((x-1)*fe/N_fft,levels(:,u),'r+')
            plot(f_lobe_interp,(y_lobe_interp),'r')
        end
        hold off
        
    levels_loc(:,u) = (x-1)*fe/N_fft;
        
    subplot(2,2,2)
        if AM        
            if dB
                plot(alpha0_dB,20*log10(levels_normalized'))
            else
                plot(alpha0_dB,levels_normalized')
            end
            title('Amplitude of the point versus the amplitude modulation rate')
            xlabel('Amplitude modulation rate (dB/s)')
            ylabel('Relative point level (dB)')
        else
            if dB
                plot(FCR,20*log10(levels_normalized'))
            else
                plot(FCR,levels_normalized')
            end
            title('Relative amplitude versus the frequency modulation rate')
            xlabel('Frequency modulation rate (Hz/s)')
            ylabel('Relative point level (dB)')
        end 
        %legend('1','2','3','4','5')
        grid on        

    subplot(2,2,4)
        if AM        
            plot(alpha0_dB,levels_loc'-f0)
            xlabel('Amplitude modulation rate (dB/s)')
            ylabel('Relative location (Hz)')
            title('Relative frequency of the points versus the amplitude modulation rate')
        else
            plot(FCR,levels_loc'-f0)
            xlabel('Frequency modulation rate (Hz/s)')
            ylabel('Relative location (Hz)')
            title('Relative frequency of the points versus the frequency modulation rate')
        end 
        grid on        
        
    pause(0.1)
end

subplot(2,2,2)
    legend('1','2','3','4','5')

subplot(2,2,4)
    legend('1','2','3','4','5')

figure
if AM        
    if dB
        plot(alpha0_dB,20*log10(levels_normalized(1:(N_lobe+1)/2,:)'))
    else
        plot(alpha0_dB,levels_normalized(1:(N_lobe+1)/2,:)')
    end
    title('Amplitude of the point versus the amplitude modulation rate')
    xlabel('Amplitude modulation rate (dB/s)')
    ylabel('Relative point level (dB)')
else
    if dB
        plot(FCR,20*log10(levels_normalized(1:(N_lobe+1)/2,:)'))
    else
        plot(FCR,levels_normalized(1:(N_lobe+1)/2,:)')
    end
    title('Relative amplitude versus the frequency modulation rate')
    xlabel('Frequency modulation rate (Hz/s)')
    ylabel('Relative point level (dB)')
end 
grid on     
    
figure
if AM        
    plot(alpha0_dB,levels_loc(1:(N_lobe+1)/2,:)'-f0)
    xlabel('Amplitude modulation rate (dB/s)')
    ylabel('Relative location (Hz)')
    title('Relative frequency of the points versus the amplitude modulation rate')
else
    plot(FCR,levels_loc(1:(N_lobe+1)/2,:)'-f0)
    xlabel('Frequency modulation rate (Hz/s)')
    ylabel('Relative location (Hz)')
    title('Relative frequency of the points versus the frequency modulation rate')
end 
grid on     


x_lobe = (levels_loc-f0)/fe;
y_lobe = levels;


if AM        
    parameter = alpha0_dB;
else
    parameter = FCR;
end 

save('data.mat','x_lobe','y_lobe','parameter')    
    
    
    