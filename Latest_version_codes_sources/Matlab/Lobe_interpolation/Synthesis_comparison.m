%% Projet long 2015
% Quentin Biache - Lili Zheng - Emilie Abia
% Lobe with different AM/FM modulation animation

% Comparison - generated lobe/theorical lobe

% This scripts generates some chirps at given frequencies, with given ACR
% and FCR, and compares it to the reconstructed version, using lobe
% interpolation.

%% Initialize workspace
clear all
close all
clc

load data.mat

%% Parameters definition

% Sampling parameters
fe = 44100;             % Sampling frequency

% FFT parameters
w_time = 0.023;         % Window duration (s)
N_padding = 8;          % Zero-padding factor
N_pts = round(fe*w_time);               % Total number of points
if (mod(N_pts,2) == 0)
    N_pts = N_pts + 1;
end

w_time = N_pts/fe;

N_fft = 2^(nextpow2(N_pts)+N_padding);  % FFT size
t = (0:(N_pts-1))'/fe;                  % Time vector
f = (0:(N_fft-1))'*fe/N_fft;            % Frequency vector

w = window(@hanning,N_pts);             % Create window
w = w./sum(w);                          % Normalize window

AM = 0;
dB = 1;

%% Compute FFT

% Number of points on the main lobe
N_lobe = size(x_lobe,1);


%% Grid settings

% Number of points
N_ACR = 3;
N_FCR = 3;

% Total number of simulations
N = N_ACR*N_FCR;

% Range 
FCR_min = -10000;
FCR_max = 10000;
ACR_min = -3/0.023;
ACR_max = 3/0.023;

% Create uniform grid of ACR and FCR parameters
ACR_q = kron(linspace(ACR_min,ACR_max,N_ACR)',ones(N_FCR,1));
FCR_q = kron(ones(N_ACR,1),linspace(FCR_min,FCR_max,N_FCR)');

% Matrix to store lobe data
x_lobe_reconst = zeros(N_lobe,N);
mod_lobe_reconst = zeros(N_lobe,N);
arg_lobe_reconst = zeros(N_lobe,N);

figure
pause

for u = 1:N

    % Parameters for the model
    f0 = 1000;              % Instantaneous frequency at t = 0
    A0 = 5;                 % Signal level
    phi0 = 0;               % Initial phase

    % Adjust parameters with the definition
    alpha0 = log(10^(ACR_q(u)/20));
    beta0 = FCR_q(u)/pi;        
    omega0 = 2*pi*f0;
    lambda0 = log(A0);

    % Create signal
    s = exp(alpha0.*t).*exp(lambda0).*exp(1i*((beta0*t.^2)+(omega0.*t)+phi0));
    
    subplot(2,2,1)
        plot(t,real(s),'b',t,exp(alpha0.*t).*exp(lambda0),'r',t,-exp(alpha0.*t).*exp(lambda0),'r')
        title('Real part of the signal')
        xlabel('Time (s)')
        ylabel('Signal')
        grid on

    % Apply window
    s = s.*w;

    % Correct phase
    s_w = zeros(N_fft,1);
    s_w(1:(N_pts-1)/2) = s(((N_pts+1)/2)+1:end);
    s_w(N_fft-(N_pts-1)/2:end) = s(1:((N_pts+1)/2));
    
   
    % Compute the FFT
    N_half = (N_fft/2)+1;                   % Half FFT size
    fft_s = fft(s_w,N_fft);
    mod_fft_s = abs(fft_s);

    plot_index = find((f > 0.8*f0) & (f < 1.2*f0));
    
    subplot(2,2,2)
        if dB
            plot(f(plot_index),20*log10(mod_fft_s(plot_index)))
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

    % Perform bilinear interpolation to find the new points on the lobe
    [lobe_freq,lobe_mag,lobe_phase] = lobe_reconst(ACR_q(u),FCR_q(u),LUT_ACR,LUT_FCR,x_lobe,mag_lobe,arg_lobe,N_lobe);

    x_lobe_reconst(:,u) = lobe_freq;
    mod_lobe_reconst(:,u) = lobe_mag;
    arg_lobe_reconst(:,u) = lobe_phase;


    % Split the frequency axis into N_lobe points in the main lobe
    x = linspace(lower_zero_loc,upper_zero_loc,N_lobe)';

    % Reconstruct signal
    Y = zeros(N_pts,1);     % Complex spectrum
    
    % Estimate the lobe phase over the synthesis domain using the 9 points
    lobe_phase = interp1(1+((x-1)*(N_pts-1)./(N_fft-1)),arg_lobe_reconst(:,u),(1:N_pts)','spline',0);
    %lobe_phase = quad_interp(1+((x-1)*(N_pts-1)./(N_fft-1)),arg_lobe_reconst(:,u),(1:N_pts)');

    % De-normalize the phase
    lobe_phase = lobe_phase + (pi*w_time*(f0-LUT_f0)) + phi0;
    
    % Estimate the lobe magnitude over the synthesis domain using the 9 points
    lobe_mag = interp1(1+((x-1)*(N_pts-1)./(N_fft-1)),mod_lobe_reconst(:,u),(1:N_pts)','spline',0);
    %lobe_mag = quad_interp(1+((x-1)*(N_pts-1)./(N_fft-1)),mod_lobe_reconst(:,u),(1:N_pts)');
    Y = A0*lobe_mag.*exp(1i*lobe_phase);
    
    s_reconst = real(ifft(Y,N_pts));

    % Undo window phase correction
    s_reconst_w = zeros(N_pts,1);
    s_reconst_w(1:(N_pts+1)/2) = s_reconst(((N_pts+1)/2):end);
    s_reconst_w(1+((N_pts+1)/2):end) = s_reconst(1:((N_pts-1)/2));
    
    s_reconst = s_reconst_w;
       
    % Plot the reconstructed magnitude
    subplot(2,2,2)
        hold on
        if dB
            %plot((x_lobe_reconst(:,u)*fe)+f0,20*log10(abs(mod_lobe_reconst(:,u))),'r') 
            plot((x_lobe_reconst(:,u)*LUT_fe)+f0,20*log10(abs(A0*mod_lobe_reconst(:,u))),'k+') 
        else
            %plot(f_lobe_interp,(y_lobe_interp),'r')
            plot((x_lobe_reconst(:,u)*LUT_fe)+f0,abs(A0*mod_lobe_reconst(:,u)),'k+') 
        end
        hold off
          
    % Plot the reconstructed signal
    subplot(2,2,3)
        plot(t,real(s),'b',t,s_reconst,'r')
        xlabel('Time (s)')
        ylabel('Signal')
        legend('Original','Reconstructed')
        err = round(100*norm(real(s)-s_reconst)./norm(real(s)));
        title(['Reconstruction error = ',num2str(err), '%'])
        grid on    

    % Plot the phase
    subplot(2,2,4)
        x_phase = linspace(lower_zero_loc,upper_zero_loc,1000)';
        lobe_index = (lower_zero_loc:upper_zero_loc)';
        
        % Theorical phase
        plot((x_phase-1)*fe/N_fft,interp1(lobe_index,unwrap(phase(fft_s(lobe_index))),x_phase),'b')
        hold on
        plot((x_phase-1)*fe/N_fft,interp1(x,arg_lobe_reconst(:,u)+(pi*w_time*(f0-LUT_f0)) + phi0,x_phase,'spline'),'k')
        xlabel('Frequency (Hz)')
        ylabel('Phase (rad)')
        title(['Current coordinates on the (ACR,FCR) plane: (',num2str(ACR_q(u)),',',num2str(FCR_q(u)),')'])
        grid on
        hold off
    
    pause
end
    
close all
        