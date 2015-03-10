# -*- coding: utf-8 -*-

"""
Created on Fri Feb 13 12:57:46 2015

@author: Emily Abia, Quentin Biache, Lili_ZHENG

This script perform the simulation of synthesis by modulated sinusoidal model. 
"""

#==============================================================================
# Import section
# Import used library
#==============================================================================
from numpy import *
from matplotlib.pyplot import * 
from scipy.io import wavfile 
from scipy import signal
from sinemodel import *
from estimation_func import *




#==============================================================================
# Test section 
#==============================================================================

# read the wav file
filename="Sound/AM.wav"
[fs,x] = wavfile.read(filename)
 
# signal size       
L = size(x)

# duration of each tram
dureefenetre = 23e-3

# window size (odd)
M = dureefenetre*fs
if (M%2==0):M+=1

# Ratio (hop size)/(window size)
H_ratio=0.25                  

# analysis window
w = signal.blackmanharris(M)

# fft size
N = 1024

# synthesis window
Ns= 1024

# threshold for peaks detection, in negative dB
seuil = -80

# perform synthesis
y = sinemodel(x,fs,w,N,Ns,H_ratio,seuil)

# records the test result
wavfile.write('test.wav',fs,y)




#==============================================================================
# output section
#==============================================================================
figure()
t = arange(L)/float(fs)
title('modulated sinusoidal synthesis test')
plot(t,x)
hold(True)
plot(t,y,'r')
ylabel('amplitude')
xlabel('time (s)')
legend(['test signal','synthesized signal'])

# compute nomalized error
error = sqrt((sum((x-y)**2.0))/(sum((x)**2.0)))
print 'total nomalized error: ',error
