%% Projet long 2015
% Quentin Biache - Lili Zheng - Emilie Abia

% This function fits a two-terms sine model to the sampled values of the
% main lobe, and returns a plot of the coefficients.

close all
clear all
clc

load data.mat

% Number of simulations (ie, number of values for the parameter)
N_simu = size(y_lobe,2);

a1 = zeros(N_simu,1);
b1 = zeros(N_simu,1);
c1 = zeros(N_simu,1);
a2 = zeros(N_simu,1);
b2 = zeros(N_simu,1);
c2 = zeros(N_simu,1);
r2 = zeros(N_simu,1);

for i = 1:N_simu
     [a,goodness] = fit(x_lobe(:,i),y_lobe(:,i),'sin2');
     a1(i,1) = a.a1;
     b1(i,1) = a.b1;
     c1(i,1) = a.c1;
     a2(i,1) = a.a2;
     b2(i,1) = a.b2;
     c2(i,1) = a.c2;  
     r2(i,1) = goodness.rsquare;
end

c1 = unwrap(c1);
c2 = unwrap(c2);

figure
    subplot(3,2,1)
    plot(parameter,a1,'b')
    grid on
    hold on
    title('a1')
    xlabel('Frequency change rate (Hz/s)')

    subplot(3,2,2)
    plot(parameter,a2,'b')
    grid on
    hold on
    title('a2')

    subplot(3,2,3)
    plot(parameter,b1,'b')
    grid on
    hold on
    title('b1')

    subplot(3,2,4)
    plot(parameter,b2,'b')
    grid on
    hold on
    title('b2')

    subplot(3,2,5)
    plot(parameter,c1,'b')
    grid on
    hold on
    title('c1')
    xlabel('Frequency change rate (Hz/s)')

    subplot(3,2,6)
    plot(parameter,c2,'b')
    grid on
    hold on
    title('c2')
    xlabel('Frequency change rate (Hz/s)')    

% Plot the R^2 coefficient to evaluate fit accuracy    
figure
    plot(parameter,r2)
    grid on
    xlabel('Parameter value')
    ylabel('R^2')
    axis([min(parameter) max(parameter) 0.8 1.2])
    
    
    
