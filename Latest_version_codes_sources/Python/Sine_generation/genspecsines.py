# -*- coding: utf-8 -*-
"""
Created on Thu Feb 12 09:15:47 2015

@author: Emily Abia, Quentin Biache, Lili_ZHENG
"""

from numpy import *
from matplotlib.pyplot import *
from genbh92lobe import *

def genspecsines(ploc, pmag, pphase, N):
#==============================================================================
#   Compute a spectrum from a series of sine values
#   iploc, ipmag, ipphase: sine locations, magnitudes and phases
#   N: size of complex spectrum
#   Y: generated complex spectrum of sines  
#==============================================================================

# initialize output spectrum
    Y = zeros(N,dtype=complex)
    
# size of positive freq. spectrum
    hN = N/2+1
    
# generate all sines lobes
    for i in xrange(size(ploc)):
        
#    location of peak
        loc = ploc[i]
        
#    avoid frequencies out of range
        if ((loc<=0) or (loc>=hN-2)):continue
        
        binremainder = round(loc) - loc

#   main lobe bin to read        
        lb = arange(9)-4+binremainder
        
#   lobe magnitude of the complex exponential
        lmag = genbh92lobe(lb,N,N)*(10**(pmag[i]/20.0))

#   spectrum bins to fill
        b = arange(9)-4+round(loc)+1

        for m in xrange(9):            
#           peak lobe croses DC bin
            if (b[m]<0):Y[-b[m]] += lmag[m]*exp(-1j*pphase[i])

#           peak lobe croses Nyquist bin
            elif (b[m]>hN-1):Y[2*(hN-1)-b[m]] +=  lmag[m]*exp(-1j*pphase[i])

#           peak lobe in positive freq. range
            else:
                Y[b[m]]+=(lmag[m]*exp(1j*pphase[i])+lmag[m]*exp(-1j*pphase[i])*(b[m]==0 or b[m]==hN-1))

#   fill the rest of the spectrum
    Y[hN:] = conjugate(Y[hN-2:0:-1])
    return Y










