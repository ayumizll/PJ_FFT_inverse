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

PLOT_ON = 0;

N_pts = round(fe*w_time);               % Total number of points
N_fft = 2^(nextpow2(N_pts)+N_padding);  % FFT size
t = (0:(N_pts-1))'/fe;                  % Time vector
f = (0:(N_fft-1))'*fe/N_fft;            % Frequency vector

w = window(@hanning,N_pts);             % Create window
w = w./sum(w);                          % Normalize window

%% Compute FFT

% Number of simulations
N = 201;

FCR = linspace(-5000,5000,N);
%alpha0_dB = linspace(-3,-20,N)/0.023;

% Matrix to store lobe data
levels = zeros(9,N);

for u = 1:N

    % Parameters for the model
    f0 = 1000;              % Instantaneous frequency at t = 0
    A0 = 1;                 % Signal level
    alpha0_dB = -0/0.023;   % Amplification rate (dB/s)
    %FCR = 000;              % Frequency change rate (Hz/s)
    phi0 = 0;               % Initial phase

    % Adjust parameters with the definition
    alpha0 = log(10^(alpha0_dB/20));
    beta0 = FCR(u)/pi;
    omega0 = 2*pi*f0;
    lambda0 = log(A0);

    % Create signal
    s = exp(alpha0.*t).*exp(lambda0).*exp(1i*((beta0*t.^2)+(omega0.*t)+phi0));

    if (PLOT_ON)
            plot(t,real(s))
            title('Real part of the signal')
            xlabel('Time (s)')
            ylabel('Signal')
            grid on
            hold on
    end

    % Apply window
    s = s.*w;

    % Compute the FFT
    N_half = (N_fft/2)+1;                   % Half FFT size
    fft_s = fft(s,N_fft);
    mod_fft_s = abs(fft_s);
    arg_fft_s = 0;%unwrap(phase(fft_s));

    if (PLOT_ON)
        subplot(2,1,1)
            plot(f,20*log(mod_fft_s))
            title('|X_{DFT}(f)|^2')
            xlabel('Frequency (Hz)')
            ylabel('|X_{DFT}(f)|^2')
            grid on
            hold on

        subplot(2,1,2)
            plot(f,arg_fft_s*180/pi)
            title('arg(X_{DFT}(f))')
            xlabel('Frequency (Hz)')
            ylabel('arg(X_{DFT}(f))')
            grid on
    end    

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

    levels(:,u) = interp1((1:N_fft)',mod_fft_s,x,'spline');


    if (PLOT_ON)
        subplot(2,1,1)
            hold on
            plot(mod_argmax,20*log(mod_max),'r+')
            plot(interp_argmax,20*log(interp_max),'g+')
            plot((x-1)*fe/N_fft,20*log(levels(:,1)),'b+')
    end

    disp([num2str(round(100*u/N)),' % achieved'])
    
end

plot(FCR,levels')
legend('1','2','3','4','5','6','7','8','9')
grid on
