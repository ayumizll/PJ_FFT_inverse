%% Projet long 2015
% Non-uniform grid generator

% Generate a non-uniform grid of control points. The grid density is based
% on the reconstruction error (so the script error_plot.m must be launched
% first)

%% Initialize worskpace
close all
clear all
clc

%% Settings

% Min and maximal number of points in a cluster
N_min = 1;
N_max = 40;

% Number of subdivisions
N_ACR = 5;
N_FCR = 5;

%% Create grid

% Initialize coordinates
X_c = [-3/0.023;+3/0.023;-3/0.023;+3/0.023;0;0;-3/0.023;3/0.023];
Y_c = [+10000;10000;-10000;-10000;10000;-10000;0;0];

% Re-scale error data
load error.mat
a = min(err(:));
b = max(err(:));
e = round(N_min + ((N_max-N_min).*(err-a)./(b-a)));

% Create primary grid
[ACR_g,FCR_g] = meshgrid(linspace(min(ACR_q(:)),max(ACR_q(:)),N_ACR)',linspace(max(FCR_q(:)),min(FCR_q(:)),N_FCR)');

ACR_g_step = ACR_g(1,2)-ACR_g(1,1);
FCR_g_step = FCR_g(1,1)-FCR_g(2,1);

for u = 1:(N_ACR-1)
    for v = 1:(N_FCR-1)
        % Number of points required for this region
        N_pts = round(griddata(ACR_q,FCR_q,e,ACR_g(v,u)+(0.5*ACR_g_step),FCR_g(v,u)-(0.5*FCR_g_step)));
        
        FCR_min = FCR_g(v+1,u);
        FCR_max = FCR_g(v,u);
        ACR_min = ACR_g(v,u);
        ACR_max = ACR_g(v,u+1);
        
        for w = 1:N_pts
            X_c = [X_c; ACR_min + (ACR_max-ACR_min).*rand];
            Y_c = [Y_c; FCR_min + (FCR_max-FCR_min).*rand];
        end
        
    end
end        

figure
    plot(X_c,Y_c,'k.','MarkerSize',15)
    xlabel('Amplitude change rate (dB/s)')
    ylabel('Frequency change rate (Hz/s)')
    title(['Nombre total de points = ',num2str(length(X_c))])
    grid on

% Interpolate the error on a uniform grid (for the 3D plot)
[ACR_grid,FCR_grid] = meshgrid(linspace(-3,3,100)'/0.023,linspace(10000,-10000,100)');
err_grid = griddata(ACR_q,FCR_q,err,ACR_grid,FCR_grid,'cubic');         
resampled_err = griddata(ACR_q,FCR_q,err,X_c,Y_c,'cubic');         

figure
    surf(ACR_grid,FCR_grid,err_grid)
    hold on
    plot3(X_c,Y_c,resampled_err,'rd','MarkerSize',8,'MarkerFaceColor','r')
    xlabel('Amplitude change rate (dB/s)')
    ylabel('Frequency change rate (Hz/s)')
    zlabel('Normalized error')
    title('Non-uniform sampling of the (ACR,FCR) plane')
    colorbar    

save('non_uniform_grid.mat','X_c','Y_c')  
        