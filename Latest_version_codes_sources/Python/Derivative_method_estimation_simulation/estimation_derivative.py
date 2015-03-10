# -*- coding: utf-8 -*-
"""
Created on Tue Feb 24 20:17:39 2015

@author: Lili_ZHENG
"""

from numpy import *
from matplotlib.pyplot import *
from scipy import signal
import time



###############################################################################
#   parameters estimator of modulated sinusoidal model (derivative method)
#   reference: GENERALIZATION OF THE DERIVATIVE ANALYSIS METHOD TO NON-STATIONARY SINUSOIDAL MODELING
#              Conference on Digital Audio Effects (DAFx-08), Espoo, Finland, September 1-4, 2008
#   input:
#       s: input signal
#       Fs: sampling frequency
#       w: analysis window
#       N_fft: fft size
#       th: threshold for peaks detection, in negative dB
#   output:
#       a1:amplitude
#       omega1:frequency in radian
#       phi1:initial phase
#       mu1:modulation of amplitude in s^-1
#       psi1:modulation of frequency in Hz/s 
###############################################################################


def estimator1(s,Fs,w,N_fft,th):
#   calculate the first and second derivative of signal
    d1 = gradient(s)*float(Fs)
    d2 = gradient(d1)*float(Fs)
    
#   size of the signal and size of the fft
    N = size(s,0)
#    N_fft = int(2**(nextpow2(N)))

#   Time axis
    hN = (N-1)/2
    t = arange(-hN,hN+1)/float(Fs)
    
#   calculate the spectrum of the signal and it's first derivative
    X = fft.fft((w*s),N_fft)
    X1 = fft.fft((w*d1),N_fft)
    X2 = fft.fft((w*d2),N_fft)
    
    #   find the pic
    mX = 20*log10(abs(X[1:N_fft/2]))
    m = array(where((mX[1:-1]>mX[2:])& (mX[1:-1]>th) & (mX[1:-1]>mX[:-2])))+1 
    m = m[0,:]
    
#   plot the spectrum 
    figure(facecolor='w')
    subplot(121)
    title('spectrum')
    plot(arange(1,N_fft/2)/float(N),mX)
    xlabel('normalized frequency')
    ylabel('magnitude')

#   estimater  omega
    r1 = X1[m]/X[m]
    omega1 = imag(r1)
    
#   initialize the output 
    L = size(omega1)
    mu1 = zeros(L)
    psi1 = zeros(L)
    phi1 = zeros(L)
    a1 = zeros(L)
    
#   estimate other parameters    
    for i in xrange(L):
        X =sum(w*s*exp(-1j*omega1[i]*t))
        X1 =sum(w*d1*exp(-1j*omega1[i]*t))
        X2 =sum(w*d2*exp(-1j*omega1[i]*t))     
          
        r1 = X1 / X
        r2 = X2 / X
        
        mu1[i] = real(r1)
        psi1[i] = imag(r2) - 2.0*mu1[i]*omega1[i]
        p = Gamma_est(N, Fs, 0, mu1[i], psi1[i], w)
        phi1[i] = angle(X/p)
        a1[i] = 2.0*abs(X/p)
    return [a1,omega1,phi1,mu1,psi1]    
    


###############################################################################
#   Gamma function used to perform estimation
#   reference: GENERALIZATION OF THE DERIVATIVE ANALYSIS METHOD TO NON-STATIONARY SINUSOIDAL MODELING
#              Conference on Digital Audio Effects (DAFx-08), Espoo, Finland, September 1-4, 2008
###############################################################################
def Gamma_est(N,Fs,Delta,mu,psi,w):
    L = size(Delta)
    Mw = zeros([L,size(w)])
    Mw[:]=w[newaxis,:]
#   time axis    
    hN = (N-1)/2
    t = arange(-hN,hN+1)/float(Fs)
    
    result = transpose (sum (Mw*exp(mu*t+1j*(Delta*t+psi*(t**2)/2.0)),1))

#   in case of problems with large mu value... 
    idx_nan = array(where(isnan(result)))
    idx_inf = array(where(isinf(result)))
    idx_nan,idx_inf=idx_nan[0,:],idx_inf[0,:]
    idx = concatenate((idx_nan,idx_inf))

#   reset to default: no correction    
    result[idx] = 1
    return result
    
    
#==============================================================================
#  Time consume 
#==============================================================================
def tic():
    #Homemade version of matlab tic and toc functions
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc():
    if 'startTime_for_tictoc' in globals():
        print "Elapsed time is " + str(time.time() - startTime_for_tictoc) + " seconds."
    else:
        print "Toc: start time not set"