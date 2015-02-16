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
dureefenetre = 23e-3                #each slice longer
#M=fs*dureefenetre
M = dureefenetre*fs                 #window size
if (M%2==0):M+=1
H_ratio=0.25                      # Ratio (hop size)/(window size)
w = signal.blackmanharris(M)      # window for slice
N = 1024
Ns= 1024
seuil = -50


#t = linspace(0,1,16000);
#x = sin(2*pi*1000*t);
y = sinemodel(x,w,N,Ns,H_ratio,seuil)
wavfile.write('test.wav',fs,y)

figure()
plot(x,'+')
plot(y,'rx')
legend(['original signal','sythesis signal'])
print sum(abs(x-y))
