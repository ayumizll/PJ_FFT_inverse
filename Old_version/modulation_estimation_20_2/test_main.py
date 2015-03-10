# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 15:28:29 2015

@author: Lili_ZHENG
"""
from func import *
from numpy import *
from matplotlib.pyplot import *
from scipy import signal

#   trame size
N = 1024

#   ampling frequency
Fs = 44100

#   test amplitude
a = 10

#   test amplitude modulation
mu = 20

#   test phase
phi = pi/2.0

#   test frequency
omega = 500

#   test frequency modulation
psi = 2000


# time axis
T = N/float(Fs)
t = arange(N)/float(Fs)-T/2.0


#   test signal
s = a*exp(mu*t+1j*(phi+omega*t+0.5*psi*t**2));
N = size(s)

#   used window and its first and second derivative
w = signal.blackmanharris(N)
wd = win_d1(N,Fs)
wdd = win_d2(N,Fs)
tw = (t)*w
twd = (t)*wd

#   compute the spectre
Xw = fft.fft(w*s)
Xwd = fft.fft(wd*s)
Xtw = fft.fft(tw*s)
Xwdd = fft.fft(wdd*s)
Xtwd = fft.fft(twd*s)

#   find the pic
tmp = abs(Xw[1:-1])
m = array(where(tmp==max(tmp)))+1
m = m[0,:]

#   frequency
delta_omega = -imag (Xwd[m] / Xw[m])
base_omega = (m)*2*pi*Fs/N;
omega1 = base_omega + delta_omega;

#   amplitude modulation
mu1 = -real (Xwd[m] / Xw[m])

#   time offeset
delta_t=real(Xtw[m]/Xw[m])

#   frequency modulation
psi1 = (imag(Xwdd[m]/Xw[m]) - imag((Xwd[m]/Xw[m])**2)) / (real((Xtw[m]*Xwd[m])/(Xw[m])**2) - real(Xtwd[m]/Xw[m]))

#   amplitude et phase
p = Gamma_est(N, Fs, delta_omega, mu, psi, w)
phi1 = angle (Xw[m] / p);
a1 = abs (Xw[m] / p);

print 'frequency original and estimated:',omega,omega1
print 'amplitude modulation original and estimated:',mu,mu1
print 'frequency modulation original and estimated:',psi,psi1
print 'amplitude original and estimated:',a,a1
print 'phase original and estimated:',phi,phi1


s1 = a1*exp(mu1*t+1j*(phi1+omega1*t+0.5*psi1*t**2));

plot(t,s.real)
hold(True)
plot(t,s1.real)
legend(['Original signal','Estimated signal'])












