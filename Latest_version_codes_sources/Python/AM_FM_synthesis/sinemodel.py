# -*- coding: utf-8 -*-
"""
Created on Thu Feb 12 13:44:30 2015

@author: Emily Abia, Quentin Biache, Lili_ZHENG
This scripte used to test systhesis with modulated sinusoidal model
"""

#==============================================================================
# Import section
#==============================================================================
from numpy import *
from matplotlib.pyplot import *
from scipy import signal
from estimation_func import *


#==============================================================================
# functions secton
#==============================================================================

def sinemodel(x,Fs,w,N,Ns,H_ratio,th):
#==============================================================================
#   Analysis/synthesis of a sound using the sinusoidal model
#   x: input sound, w: analysis window (odd size), N: FFT size,
#   th: threshold in negative dB, y: output sound
#   Ns: FFT size for synthesis (even)
#   H_ratio: Ratio (hop size)/(window size)
#==============================================================================

#   data nomalization
    max_x = 32767.0
    x = x/max_x    
    
#   analysis window size    
    M = size(w)

#   Hop size
    H = int(H_ratio*M)
    
#   size of positive spectrum
    N2 = N/2+1
    
#   length of data
    soundlength = size(x)
    
#   half of synthesis window size
    hNs = Ns/2
    
#   half analysis window size
    hM = (M-1)/2

#   Initialization of sound pointer    
    pin = max(hNs, hM)

#   reduce effect of limite
    x = concatenate((zeros(M),x,zeros(M)))
    
#   last sample to start a frame
    pend = size(x) - pin -1
    
#   Initialization of buffer for FFT
    fftbuffer = zeros(N)
    
#   Initialization of output vector
    y = zeros(size(x))
    out = zeros(soundlength,dtype=int16)

#   normalisation of windows    
    w = w/float(sum(w))
    
    sw = zeros(Ns)
    
#   overlap index
    ovidx = arange(2*H-1)+Ns/2-H+1

#   overlapping window 
    ow = signal.triang(2*H-1)
    sw[ovidx] = ow[:2*H-1]
    
#   synthesis window
    bh = signal.blackmanharris(Ns)
    bh = bh/float(sum(bh))
    sw[ovidx] = sw[ovidx]/bh[ovidx]
        
#   time axis
    t = arange(-hM,hM+1)/float(Fs) 
        
    while pin <= pend:

#==============================================================================
#   analysis
#==============================================================================
#       tram shift
        xw = x[pin-hM:pin+hM+1]
        
        
###############################################################################        
#       perform modulation estimation for each tram
#       a1:amplitude
#       omega1:frequency in radian
#       phi1:initial phase
#       mu1:modulation of amplitude in s^-1
#       psi1:modulation of frequency in Hz/s 
###############################################################################
        [a1,omega1,phi1,mu1,psi1]  = estimator1(xw,Fs,w,N,th)
#        [a1,omega1,phi1,mu1,psi1]  = estimator2(xw,Fs,th)  

#       generation of sinusoidal from estimated parameters        
        L = size(a1)
        s = zeros(M)
        for i in xrange(L):
            if (abs(mu1[i])>100000 or abs(psi1[i])>100000):
                pass
            else:            
                s +=(a1[i]*exp(mu1[i]*t+1j*(phi1[i]+omega1[i]*t+0.5*psi1[i]*t**2))).real
                
#       add window to each tram               
        xw = s*w

#       zero-phase window in buffer
        fftbuffer[:]=0
        fftbuffer[:hM+1] = xw[hM:]
        fftbuffer[-hM:] = xw[:hM]
        

#       compute FFT 
        X = fft.fft(fftbuffer)
        mX = 20*log10(abs(X[:N2]))
        pX = unwrap(angle(X[:N2]))
        
#       find peaks
        ploc = array(where((mX[1:-1]>mX[2:])& (mX[1:-1]>th) & (mX[1:-1]>mX[:-2])))+1
        
        plocs = (ploc)*float(Ns)/N-1.0
        pmag = mX[ploc]
        pphase = pX[ploc]
        
#       generate spec sines 
        Y = genspecsines(plocs[0,:], pmag[0,:], pphase[0,:], Ns)
        
#       restore signal from synthesized mainlobe
        yw = fft.fftshift(real(fft.ifft(Y))) 

#       overlap add        
        y[pin-hNs:pin+hNs]+=sw*yw
        
#       avance hop-size        
        pin += H
        
    out[:]=y[M:-M]*max_x      
    return out
    
    
    
    
    
    
    
    
    