# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 15:31:33 2015

@author: Lili_ZHENG
"""

from numpy import *
from matplotlib.pyplot import *
from scipy import signal


def win_d1(N,Fs):
#   time axis
    T = N/float(Fs)
    t = arange(N)/float(Fs)-T/2.0
    

#   coef of blackman_harris window
    a1 = 0.48829;
    a2 = 0.14128;
    a3 = 0.01168;
    
    wd = -a1*(2*pi/T)*sin(2*pi*t/T)-a2*(4*pi/T)*sin(4*pi*t/T)-a3*(6*pi/T)*sin(6*pi*t/T)
    
    return wd
    
def win_d2(N,Fs):
#   time axis
    T = N/float(Fs)
    t = arange(N)/float(Fs)-T/2.0

#   coef of blackman_harris window
    a1 = 0.48829;
    a2 = 0.14128;
    a3 = 0.01168;
    
    wdd = -a1*(2*pi/T)**2*cos(2*pi*t/T)-a2*(4*pi/T)**2*cos(4*pi*t/T)-a3*(6*pi/T)**2*cos(6*pi*t/T)
    
    return wdd
    
def Gamma_est(N,Fs,Delta,mu,psi,w):
    L = size(Delta)
    Mw = zeros([L,size(w)])
    Mw[:]=w[newaxis,:]
#   time axis    
    T = N/float(Fs)
    t = arange(N)/float(Fs)-T/2.0
    
    result = transpose (sum (Mw*exp(mu*t+1j*(Delta*t+psi*(t**2)/2.0)),1))
    
#   in case of problems with large mu value... 
    idx = concatenate((where(isnan(result)), where(isinf(result)))) 

#   reset to default: no correction    
    result[idx] = 1
    return result
    
    
    
    