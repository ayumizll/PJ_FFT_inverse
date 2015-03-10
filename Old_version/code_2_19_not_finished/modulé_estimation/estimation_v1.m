clc
clear all
close all

p=2
simga2 = 2.5;

x0 = 5;
x1 = 6;
x2 = 7;

y0 = -0.4605;
y1 = -0.4015;
y2 = -0.4057;

p0 = -0.08318;
p1 = -0.1697;
p2 = -0.2088;

w0_tild=calc_w0(x1,y0,y1,y2)
lambda_tild = est_value(x0,x1,x2,y0,y1,y2,w0_tild);
phi_tild = est_value(x0,x1,x2,p0,p1,p2,w0_tild);
% 
alpha_tild2 = est_ACR(x0,x1,x2,p0,p1,p2,w0_tild,p)
beta_tild2 = est_FCR(x0,x1,x2,y0,y1,y2,p0,p1,p2,w0_tild,p)
[w0_tild2,lambda_tild2,phi_tild2]=est_correction(w0_tild,lambda_tild,phi_tild,alpha_tild2,beta_tild2,p)
% 
