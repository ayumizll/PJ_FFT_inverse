# -*- coding: utf-8 -*-
"""
Created on Fri Feb 13 12:57:46 2015

@author: Emily Abia, Quentin Biache, Lili_ZHENG
"""

from numpy import *
from matplotlib.pyplot import * 
from scipy.io import wavfile 
from scipy import signal
from sinemodel import *

# read the wav file
filename="sound/Wave_5.wav"
[fs,x] = wavfile.read(filename)

# parameter            
L = size(x)
dureefenetre = 23e-3                #duration of each window
#M=fs*dureefenetre
M = dureefenetre*fs                 #window size
if (M%2==0):M+=1
H_ratio=0.25                      # Ratio (hop size)/(window size)
w = signal.blackmanharris(M)      # window for tram
N = 1024                          # fft size
Ns= 1024                          # window size for synthesis 
seuil = -50


#t = linspace(0,1,16000);
#x = sin(2*pi*1000*t);
y = sinemodel(x,w,N,Ns,H_ratio,seuil)
#==============================================================================
#   Analysis/synthesis of a sound using the sinusoidal model
#   x: input sound, w: analysis window (odd size), N: FFT size,
#   seuil: threshold in negative dB, y: output sound
#   Ns: FFT size for synthesis (even)
#   Ratio: (hop size)/(window size) 
#==============================================================================
wavfile.write('test.wav',fs,y)

figure()
plot(x,'+')
plot(y,'rx')
legend(['original signal','sythesis signal'])
print sum(abs(x-y))
