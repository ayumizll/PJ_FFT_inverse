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

def sine_generator(a0,f0,phi0,L,Fs,w,N,Ns,H_ratio):
#==============================================================================
#   generate sinusoidal with model a0*exp(1j*(2.0*pi*f0*t+phi0)).real
#   a0: amplitude of sine
#   f0: frequency of sine
#   phi0: phase of sine
#   L: signal size
#   w: analysis window (odd size), N: FFT size,
#   t: threshold in negative dB, y: output sound
#   Ns: FFT size for synthesis (even)
#==============================================================================
    max_x = 32767.0#float(max(x))
    a0 = a0/max_x    
    
#   analysis window size    
    M = size(w)

#   Hop size
    H = int(H_ratio*M)
       
#   half of synthesis window size
    hNs = Ns/2
    
#   half analysis window size
    hM = (M-1)/2

#   Initialization of sound pointer    
    pin = max(hNs, hM)
    
#   last sample to start a frame
    pend = L+2*M - pin -1
    
#   Initialization of output vector
    y = zeros(L+2*M)
    out = zeros(L,dtype=int16)

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
    
    
    ploc,pmag,pphase = zeros(size(f0)),zeros(size(f0)),zeros(size(f0))
    ploc[:] = f0*N/float(Fs)
    pmag[:] = 20*log10(abs(a0/2.0))
    pphase[:] = phi0+2*pi*f0*(hNs-M)/float(Fs)
    
    while pin <= pend:                
#==============================================================================
#       Synthesis  
#==============================================================================
        
#        adapt peak location to synthesis FFT
        plocs = (ploc)*float(Ns)/N-1.0
        
#       generate spec sines 
        Y = genspecsines(plocs, pmag, pphase, Ns)
        pphase = pphase + 2.0*pi*f0*H/float(Fs)

#       compute time domain sines
        yw = fft.fftshift(real(fft.ifft(Y)))  
        
#       over-lap  add
        y[pin-hNs:pin+hNs]+=sw*yw
#        print pin
        if pin==max(hNs, hM):
            #test code
            figure(facecolor='w')
            title('synthetized mainlobe')
            plot(arange(size(Y))/float(size(Y)),20*log10(abs(Y)))
            ylabel('amplitude (dB)')
            xlabel('normalized frequency')
            pass
#       avance hop-size        
        pin += H
    out[:]=y[M:-M]*max_x        
    return out

    
    
    
    
    
    
    
    
    
    