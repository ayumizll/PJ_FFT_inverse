# -*- coding: utf-8 -*-
"""
Created on Tue Feb 24 20:16:21 2015

@author: Quentin Biache, Emilie Abia, Lili_ZHENG
This script simulate the estimation of 5 parameters of mudulated sinus model
such as a*exp(mu*t+1j*(phi+omega*t+0.5*psi*t^2))
"""

#==============================================================================
# import section
#==============================================================================
from estimation_derivative import *
from numpy import *
from matplotlib.pyplot import *
from scipy import signal
from scipy.io import wavfile 


#################initialize reference signal###################################
# read the wav file
filename="sound/Wave_1.wav"
[Fs,x1] = wavfile.read(filename)
filename="sound/Wave_5.wav"
[Fs,x2] = wavfile.read(filename)
filename="sound/Wave_9.wav"
[Fs,x3] = wavfile.read(filename)
x = (x1*1.0+x3*1.0+x2*1.0)
x = x/32767.0


###################definition of parameters####################################

#   signal size
L = size(x)

#   duration of each tram 
dureefenetre = 23e-3     

#   Window size: M=fs*dureefenetre (odd)
M = round(dureefenetre*Fs)              
if (M%2==0):M+=1
    
#   Ratio (hop size)/(window size) 
H_ratio=0.25

#   Window                      
w = signal.blackmanharris(M)

#   size of fft
N = 4096                          

#   threshold for chose of peaks in spectral domain (dB)
th = -30

#   Hop size
H = int(H_ratio*M)
  

#   half analysis window size
hN = (M-1)/2

#   time axis
t = arange(-hN,hN+1)/float(Fs) 

#   add modulation to reference signal
tmp = x[1050:M+1050]*exp(-100*t+100*1j*t**2)

#   performe estimation of parameters
[a1,omega1,phi1,mu1,psi1]  = estimator1(tmp,Fs,w,N,th)


#   recorver the signal with 5 estimated parameters
L = size(a1)
s = zeros(M)
for i in xrange(L):
    s +=(a1[i]*exp(mu1[i]*t+1j*(phi1[i]+omega1[i]*t+0.5*psi1[i]*t**2))).real

#==============================================================================
# output section
#==============================================================================
subplot(122)
title('comparison')
plot(t,tmp.real)
hold(True)
plot(t,s.real,'r')
xlabel('time')
ylabel('amplitude')
legend(['orignal signal','estimated signal'])
print 'mu estimated',mu1
print 'omega estimated',omega1
print 'psi estimated',psi1
print 'phi estimated',phi1
print 'a estimated',a1

























