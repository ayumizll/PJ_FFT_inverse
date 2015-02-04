# -*- coding: utf-8 -*-
"""
Created on Mon Feb  2 21:26:53 2015

@author: Quentin Biache, Emily Abia, Lili Zheng
"""
import numpy as np
import matplotlib.pyplot as plt 
from scipy.io import wavfile 
#import stft_v1 as stft
import stft_v3 as stft

# read the wav file
filename="sound/alpha.wav"
[fs,data] = wavfile.read(filename)

#t = np.linspace(0,2,2000)
#data = np.sin(2*3.14*10*t)*15000

# parameter
#N=2048                               #fft size  
#if (N%2==1):N+=1              
dureefenetre = 20e-3                #each slice longer
M=fs*dureefenetre                   #window size
if (M%2==0):M+=1
H_ratio=0.25                         #Ratio (hop size)/(window size)
w1 = np.hanning(M)                  # window for slice
#w2 = np.ones(M)                     # window for ifft
w2 = np.hanning(M)
#w2 = np.hamming(M)

out=stft.stft(data,fs,H_ratio,w1,w2)
#==============================================================================
# Function stft used to do analysis and synthesis
# data : input data
# fs : sampling frequency
# H : Ratio (hop size)/(window size)
# w1 : window ---> slice of original signal
# w2 : window ---> slice of restored signal
#==============================================================================