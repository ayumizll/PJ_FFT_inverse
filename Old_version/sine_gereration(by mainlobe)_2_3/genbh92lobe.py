# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 22:08:29 2015

@author: Emily Abia, Quentin Biache, Lili_ZHENG
"""
from numpy import *
from matplotlib.pyplot import *

#==============================================================================
# x: bin positions to compute (real values)
# y: transform values
# M: Window size
# N: fft size
#==============================================================================
def genbh92lobe(x,M,N):
    
    x = array(x)
  
#   frequency sampling 
    f = 2.0*x*pi/N
    df = 2.0*pi/N
    r = float(N)/M
    
#   initialize window 
    y = zeros(size(x))
    
#   window constants 
    consts = array([0.35875, 0.48829, 0.14128, 0.01168])

#    sum Dirichlet kernels 
#    for m in xrange(4):
#        y = y + 0.5*consts[m]*(D(f-df*m,N)+D(f+df*m,N))
    y = 0.5*consts[0]*(D(f-df*0*r,M)+D(f+df*0*r,M)) \
      + 0.5*consts[1]*(D(f-df*1*r,M)+D(f+df*1*r,M)) \
      + 0.5*consts[2]*(D(f-df*2*r,M)+D(f+df*2*r,M)) \
      + 0.5*consts[3]*(D(f-df*3*r,M)+D(f+df*3*r,M))
        
#   normalize
    y = (y)/((M)*consts[0])
    return abs(y)
    
def D(x,M):
#   Calculate rectangular window transform (Dirichlet kernel) 
    y = sin(M*x/2.0)/sin(x/2.0)
#    y = y*exp(-1j*x*(N-1)/2.0)
    
#   avoid NaN if x==0
    y[where(y!=y)]=M
    return y
    
#n=77
#t = arange(n)-(n-1)/2
#y = genbh92lobe(t,128)
#plot(t,20*log10(abs(y)))
