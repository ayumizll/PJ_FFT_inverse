# -*- coding: utf-8 -*-
"""
Created on Fri Feb 13 12:57:46 2015

@author: Emily Abia, Quentin Biache,Lili_ZHENG
This script used to simulate function sine_generator. Give the initial frequency, phase and
amplitude (they can be a combinaison of servral sines), this function can synthetize the mainlobe of this sine in spectral domain, the 
calculate it's fft-1 to recover the sine signal.
"""

from numpy import *
from matplotlib.pyplot import * 
from scipy.io import wavfile 
from scipy import signal
from sinemodel import sine_generator



#   duration of reference signal   
T = 0.01
#   sampling frequency
Fs = 44100
#   time axis
L = T*Fs 
t = arange(L)/float(Fs)

# amplitude of reference signal
a0 = array([0.3,0.5,0.7])
# frequency of reference signal
f0 = array([1000.0,500,2000.0])
# phase of reference signal
phi0=array([1.0,3.0,0.5])

#==============================================================================
#   generate references sinusoidal
#==============================================================================
x = zeros(L)
# generation of reference signal
for i in xrange(size(a0)):
    x += (a0[i]*exp(1j*(2.0*pi*f0[i]*t+phi0[i]))).real


# parameter            
dureefenetre = 23e-3                #each slice longer
#M=fs*dureefenetre
M = dureefenetre*Fs                 #window size
if (M%2==0):M+=1
H_ratio=0.25                      # Ratio (hop size)/(window size)
w = signal.blackmanharris(M)      # window for slice
N = 1024
Ns= 1024

#==============================================================================
#   synthetize sine 
#   a0,f0,phi0 are parameters of sinemodel a0*exp(1j*(2.0*pi*f0*t+phi0)),
#   L: signal size
#   Fs: sampling frequency
#   w: window used -> blackman-harris
#   N: FFT_size
#   H_ratio: Ratio (hop size)/(window size)
#==============================================================================
#   En cas de non-normalise
if a0.any >1:
    a0=a0/32767.0

y = sine_generator(a0,f0,phi0,L,Fs,w,N,N,H_ratio)


if a0.any >1:
    y=y*32767.0 

#==============================================================================
# output
#==============================================================================
figure(facecolor='w')
plot(t,x)
plot(t,y,'rx')
grid(True)
legend(['reference signal','synthetized signal'])
xlabel('time (s)')
ylabel('amplitude')
print sqrt(sum(abs(x-y)**2)/sum(x**2))