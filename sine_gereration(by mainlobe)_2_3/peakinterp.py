# -*- coding: utf-8 -*-
"""
Created on Thu Feb 12 13:14:01 2015

@author: Emily Abia, Quentin Biache, Lili_ZHENG
"""
from numpy import *
from matplotlib import *
from scipy.interpolate import interp1d


def peakinterp(mX,pX,ploc):
#==============================================================================
#   interpolation of spectral peak 
#   mX: magnitude spectrum, pX: phase spectrum, ploc: locations of peaks
#   iploc, ipmag, ipphase: interpolated values
#   note that ploc values are assumed to be between 2 and length(mX)-1    
#==============================================================================

#   magnitude of peak bin
    val = mX[ploc]

#   magnitude of bin at lef    
    lval = mX[ploc-1]

#    magnitude of bin at right
    rval = mX[ploc+1]

#   center of parapola
    iploc = ploc+0.5*(lval-rval)/(lval-2*val+rval)
    
#   magitude of peaks
    ipmag = val-0.25*(lval-rval)*(iploc-ploc)

#   phase of peaks
    f = interp1d(arange(size(pX)), pX,kind='linear')
    ipphase = f(iploc)  
    
    return [iploc,ipmag,ipphase]
    