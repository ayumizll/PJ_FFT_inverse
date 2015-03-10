%% Projet long 2015
% Quentin Biache - Lili Zheng - Emilie Abia
% Script generates the look-up table for the main lobe on a uniform (or
% non-uniform) grid. The range can be adjusted by the user.

%% Initialize workspace
clear all
close all
clc

%% Parameters definition

% Sampling parameters
fe = 44100;             % Sampling frequency

% FFT parameters
w_time = 0.023;         % Window duration (s)
N_padding = 9;          % Zero-padding factor
N_pts = round(fe*w_time);               % Total number of points
if (mod(N_pts,2) == 0)
    N_pts = N_pts + 1;
end

N_fft = 2^(nextpow2(N_pts)+N_padding);  % FFT size
t = (0:(N_pts-1))'/fe;                  % Time vector
f = (0:(N_fft-1))'*fe/N_fft;            % Frequency vector

w = window(@hanning,N_pts);             % Create window
w = w./sum(w);                          % Normalize window

dB = 1;

%% Grid settings

% Number of points
N_ACR = 5;
N_FCR = 5;
p_max = N_ACR*N_FCR;

% Range 
FCR_min = -10000;
FCR_max = 10000;
ACR_min = -3/0.023;
ACR_max = 3/0.023;

% Uncomment to create an uniform grid (error analysis)
LUT_ACR = kron(linspace(ACR_min,ACR_max,N_ACR)',ones(N_FCR,1));
LUT_FCR = kron(ones(N_ACR,1),linspace(FCR_min,FCR_max,N_FCR)');

% Uncomment to create non-uniform grid (obtained from the error analysis)
% load non_uniform_grid.mat
% LUT_ACR = X_c;
% LUT_FCR = Y_c;
% p_max = length(X_c);

%% Generate chirps

% Number of points on the main lobe
N_lobe = 9;

% Matrix to store lobe data
levels = zeros(N_lobe,N_ACR*N_FCR);
levels_normalized = zeros(N_lobe,N_ACR*N_FCR);
levels_loc = zeros(N_lobe,N_ACR*N_FCR);
levels_phase = zeros(N_lobe,N_ACR*N_FCR);

x_lobe = zeros(N_lobe,N_ACR*N_FCR);
mag_lobe = zeros(N_lobe,N_ACR*N_FCR);
arg_lobe = zeros(N_lobe,N_ACR*N_FCR);

figure
pause

for p = 1:p_max
        
    % Parameters for the model
    f0 = 1000;              % Instantaneous frequency at t = 0
    A0 = 1;                 % Signal level
    phi0 = 0;               % Initial phase

    % Frequency plot domain
    f_plot = [0.8*f0;1.2*f0];        

    % Adjust parameters with the definition
    alpha0 = log(10^(LUT_ACR(p)/20));
    beta0 = LUT_FCR(p)/pi;        
    omega0 = 2*pi*f0;
    lambda0 = log(A0);

    % Create signal
    s = exp(alpha0.*t).*exp(lambda0).*exp(1i*((beta0*t.^2)+(omega0.*t)+phi0));

    subplot(2,2,1)
        plot(t,real(s),'b',t,exp(alpha0.*t).*exp(lambda0),'r',t,-exp(alpha0.*t).*exp(lambda0),'r')
        xlabel('Time (s)')
        ylabel('Re[s(t)]')
        title(['Real part - FCR = ',num2str(LUT_FCR(p)),' Hz/s - ACR = ',num2str(LUT_ACR(p)),' dB/s'])
        grid on

    % Apply window
    s = s.*w;

    % Correct phase
    s_w = zeros(N_fft,1);
    s_w(1:(N_pts-1)/2) = s(((N_pts+1)/2)+1:end);
    s_w(N_fft-(N_pts-1)/2:end) = s(1:((N_pts+1)/2));

    % Compute the FFT
    fft_s = fft(s_w,N_fft);
    mod_fft_s = abs(fft_s);
    arg_fft_s = 0;%unwrap(phase(fft_s));

    plot_index = find((f > f_plot(1)) & (f < f_plot(2)));

    % Plot the FFT magnitude
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

    % Search for the main lobe
    [mod_max,index] = max(mod_fft_s);
    mod_argmax = f(index);

    % Find a more accurate value of the main lobe peak
    interp_index = quad_argmax(index,mod_fft_s(index-1),mod_fft_s(index),mod_fft_s(index+1));
    interp_max = quad_max(mod_fft_s(index-1),mod_fft_s(index),mod_fft_s(index+1));
    interp_argmax = interp1(1:length(f),f,interp_index);

    % Initialize tests to find the zeros
    left_peak = mod_fft_s(1:end-2) > mod_fft_s(2:end-1);
    right_peak = mod_fft_s(3:end) > mod_fft_s(2:end-1);

    % Find location of zeros in the FFT
    mod_zeros_loc = find(left_peak & right_peak) + 1; 

    % Find the closest zero to the lobe peak
    [~,index_min] = min(abs(index-mod_zeros_loc)); 

    % Set the lower and the upper zero of the lobe
    if (mod_zeros_loc(index_min) > index)
        lower_zero_loc = mod_zeros_loc(index_min - 1); 
        upper_zero_loc = mod_zeros_loc(index_min); 
    else
        lower_zero_loc = mod_zeros_loc(index_min); 
        upper_zero_loc = mod_zeros_loc(index_min + 1); 
    end    

    % Index where the main lobe is located
    lobe_index = (lower_zero_loc:upper_zero_loc)';

    % Split the lobe in N_lobe points
    x = linspace(lower_zero_loc*1.001,upper_zero_loc*0.999,N_lobe)';
    x_phase = linspace(lower_zero_loc,upper_zero_loc,1000)';

    levels(:,p) = interp1((1:N_fft)',mod_fft_s,x);
    levels_phase(:,p) = interp1(lobe_index,unwrap(phase(fft_s(lobe_index))),x);

    x_lobe_interp = linspace(lower_zero_loc,upper_zero_loc,N_pts)';
    f_lobe_interp = (x_lobe_interp-1)*fe/N_fft;
    y_lobe_interp = quad_interp(x,levels(:,p),x_lobe_interp);            


    % Reconstruct signal
    Y = zeros(N_pts,1);

    % Estimate the lobe phase over the synthesis domain using the 9 points
    lobe_phase = interp1(1+((x-1)*(N_pts-1)./(N_fft-1)),levels_phase(:,p),(1:N_pts)','spline',0);
    %lobe_phase = quad_interp(1+((x-1)*(N_pts-1)./(N_fft-1)),squeeze(levels_phase(u,v,:)),(1:N_pts)');

    % Estimate the lobe magnitude over the synthesis domain using the 9 points
    Y = interp1(1+((x-1)*(N_pts-1)./(N_fft-1)),levels(:,p),(1:N_pts)','spline',0).*exp(1i*lobe_phase);
    %Y = quad_interp(1+((x-1)*(N_pts-1)./(N_fft-1)),squeeze(levels(u,v,:)),(1:N_pts)').*exp(1i*lobe_phase);
    s_reconst = real(ifft(Y,N_pts));

    s_reconst_w = zeros(N_pts,1);
    s_reconst_w(1:(N_pts+1)/2) = s_reconst(((N_pts+1)/2):end);
    s_reconst_w(1+((N_pts+1)/2):end) = s_reconst(1:((N_pts-1)/2));

    s_reconst = s_reconst_w;

    subplot(2,2,3)
        plot(t,real(s),'b',t,s_reconst,'r')
        xlabel('Time (s)')
        ylabel('Signal')
        legend('Original','Reconstructed')
        err = round(100*norm(real(s)-s_reconst)./norm(real(s)));
        title(['Reconstruction error = ',num2str(err), '%'])
        grid on

    subplot(2,2,2)
        hold on
        if dB
            plot((x-1)*fe/N_fft,20*log10(levels(:,p)),'r+')
            plot(f_lobe_interp,20*log10(abs(y_lobe_interp)),'r')
            plot((0:(N_pts-1))'*fe/N_pts,20*log10(abs(Y)),'k+')
            axis([f_plot(1) f_plot(2) -80 5])
        else
            plot((x-1)*fe/N_fft,levels(:,p),'r+')
            plot(f_lobe_interp,(y_lobe_interp),'r')
        end
        hold off

    % Plot the phase
    subplot(2,2,4)
        plot((x_phase-1)*fe/N_fft,interp1(lobe_index,unwrap(phase(fft_s(lobe_index))),x_phase),'b',(x-1)*fe/N_fft,levels_phase(:,p),'r+')
        hold on
        plot((x_phase-1)*fe/N_fft,interp1(x,levels_phase(:,p),x_phase,'spline'),'k')
        xlabel('Frequency (Hz)')
        ylabel('Phase (rad)')
        title('Lobe phase and reconstructed lobe phase')
        grid on
        hold off

    levels_loc(:,p) = (x-1)*fe/N_fft;

    % Store references
    x_lobe(:,p) = (levels_loc(:,p)-f0)/fe;  % Frequency normalization
    mag_lobe(:,p) = levels(:,p);            % Store magnitude
    arg_lobe(:,p) = levels_phase(:,p);      % Store phase

    pause(0.1)
end

% Unwrap along the lines ACR = constant to get rid of phase jumps (see
% report)
for t = 1:N_lobe
    y = arg_lobe(t,:);
    y = y';
    w = reshape(y,N_FCR,N_FCR);
    q = unwrap(w);
    arg_lobe(t,:) = q(:)';
end

% Parameters of the synthesized lobes (must be stored also)
LUT_f0 = f0;
LUT_fe = fe;

save('data.mat','x_lobe','mag_lobe','arg_lobe','LUT_ACR','LUT_FCR','LUT_f0','LUT_fe')    
disp('"data.mat" generated successfully')

close all

    
    