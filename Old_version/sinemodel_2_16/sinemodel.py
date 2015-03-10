# -*- coding: utf-8 -*-
"""
Created on Thu Feb 12 13:44:30 2015

@author: Emily Abia, Quentin Biache, Lili_ZHENG
"""

from numpy import *
from matplotlib.pyplot import *
from scipy import signal
from peakinterp import *
from genspecsines import *

def sinemodel(x,w,N,Ns,H_ratio,t):
#==============================================================================
#   Analysis/synthesis of a sound using the sinusoidal model
#   x: input sound, w: analysis window (odd size), N: FFT size,
#   t: threshold in negative dB, y: output sound
#   Ns: FFT size for synthesis (even)
#==============================================================================
    max_x = 32767.0#float(max(x))
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
    
    while pin <= pend:
#==============================================================================
#   analysis
#==============================================================================
#   add window to our signal
        xw = x[pin-hM:pin+hM+1]*w
        
#   zero-phase window in buffer
        fftbuffer[:]=0
        fftbuffer[:hM+1] = xw[hM:]
        fftbuffer[-hM:] = xw[:hM]
        
#       compute FFT 
        X = fft.fft(fftbuffer)
        mX = 20*log10(abs(X[:N2]))
        pX = unwrap(angle(X[:N2]))
        
#       find peaks
        ploc = array(where((mX[1:-1]>mX[2:])& (mX[1:-1]>t) & (mX[1:-1]>mX[:-2])))+1
          
#       refine peak value
        [ploc,pmag,pphase]=peakinterp(mX,pX,ploc)  
        
#==============================================================================
#       Synthesis  
#==============================================================================
        
#        adapt peak location to synthesis FFT
        plocs = (ploc)*float(Ns)/N-1.0
        
#       generate spec sines 
        Y = genspecsines(plocs[0,:], pmag[0,:], pphase[0,:], Ns)

#       compute time domain sines
        yw = fft.fftshift(real(fft.ifft(Y)))  
        
#       over-lap  add
        y[pin-hNs:pin+hNs]+=sw*yw
        
        if pin==max(hNs, hM)+10*H:
            #test code
            figure
            plot(20*log10(abs(Y)),'r')
            plot(20*log10(abs(X)))
            legend(['spectrum analysis','spectrum synthesis'])
            pass
#       avance hop-size        
        pin += H
    out[:]=y[M:-M]*max_x        
    return out

#    y = y*(max_x/max(abs(y)))
#    out[:]=y[:] #output
    
    
    
    
    
    
    
    
    
    