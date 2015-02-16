# -*- coding: utf-8 -*-
"""
Created on Tue Feb  3 23:40:21 2015

@author: Quentin Biache, Emily Abia, Lili Zheng
"""

import numpy as np
import matplotlib.pyplot as plt 
from scipy.io import wavfile
import chose_pic as pic


def stft(data,fs,H_ratio,w1,w2):
#==============================================================================
# Function stft used to do analysis and synthesis v-3
# data : input data
# fs : sampling frequency
# H : Ratio (hop size)/(window size)
# w1 : window ---> slice of original signal
# w2 : window ---> slice of restored signal
#==============================================================================    
    s_win = np.size(w1) #window size
    n1 = s_win*H_ratio  #Hop size
    n2 = n1
    N = 2**(nextpow2(s_win)+1) #fft size ---> zero padding
    
    #initialisation
    L = np.size(data)
    fftbuffer = np.zeros(N)
    out = np.zeros(L,dtype=np.int16)
    # zero-padding and normalization
    data_in = np.append(np.append(np.zeros(s_win,dtype=np.int16),data),np.zeros(s_win-L%n1,dtype=np.int16))/np.float(np.max(np.abs(data)))  
    data_out = np.zeros(np.size(data_in))
    
    #indice
    pin = 0
    pout= 0
    pend = np.size(data_in)-s_win
    
    #algo
    while pin<pend:
        #analysis
        grain = data_in[pin:pin+s_win]*w1
        fftbuffer[:] = 0                            #clear buffer
        fftbuffer[:s_win]=grain[:]                  #zeropadding
        f = np.fft.fft(np.fft.fftshift(fftbuffer))  #fft
        r = 20*np.log10(np.abs(f[:N/2+1])) #magnitude
        phi = np.unwrap(np.angle(f[:N/2+1])) #phase
        #synthesis
        ft = np.append(10**(r/20)*np.exp(1j*phi),(10**(r[-1:0:-1]/20))*np.exp(-1j*phi[-1:0:-1]))
        fftbuffer=np.fft.ifftshift(np.real(np.fft.ifft(ft)))
        grain = fftbuffer[:s_win]*w2 #ifft
        data_out[pout:pout+s_win]=data_out[pout:pout+s_win]+grain #overlap-add
        pin+=n1
        pout+=n2
    data_out = data_out/np.max(np.abs(data_out))
    out[:]=data_out[s_win:s_win+L]*np.max(np.abs(data)) #output
    
    #result
    wavfile.write('sound/test.wav',fs,out)
    plt.figure(facecolor=[1,1,1])
    plt.plot(np.array(range(L),dtype=np.float)/fs,data,'x')
    plt.hold(True)
    plt.plot(np.array(range(L),dtype=np.float)/fs,out,'r+')
    plt.xlabel('time')
    plt.ylabel('magnitude')
    plt.legend(['original signal','restored signal'])   
    plt.xlim(0,np.float(L)/fs)
    err = np.sum(np.abs(data-out))
    print 'total err is',err

    return out
    
def nextpow2(i):
    if (i == 0): return -float ('inf') 
    elif (i <= 1): return 0.0
    else:  
        n = 1
    while n < i: n *= 2
    return np.log(n)/np.log(2)