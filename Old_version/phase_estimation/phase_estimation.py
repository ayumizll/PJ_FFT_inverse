# -*- coding: utf-8 -*-
# =============================================================================
# Instantaneous phase estimation
# @author: Quentin
# Version 1.0
# =============================================================================

# Import NumPy library
import numpy as np

# =============================================================================
# Function to find the closest power of 2 for a given real number.
# This is useful to add zero-padding when computing the FFT

def nextpow2(i):
    if (i == 0): return -float ('inf') 
    elif (i <= 1): return 0.0
    else:  
        n = 1
        while n < i: n *= 2
        return np.log(n)/np.log(2)


def princ_arg(phi):
    temp = phi    
    if (phi > np.pi):
        while (temp > np.pi):
            temp -= 2*np.pi
    elif (phi <= -np.pi):
        while (temp < -np.pi):
            temp += 2*np.pi
    return temp
    

def quadratic_max(u,v,w):    
    temp = pow(u,2) + (16*pow(v,2)) + pow(w,2)
    temp += (-8*u*v) + (-8*v*w) + (-2*u*w)
    temp = temp/(u-(2*v)+w)
    temp *= (-1/8)
    return temp  

    
def quadratic_argmax(location,u,v,w):    
    temp = u-w
    temp = temp/(u-(2*v)+w)
    temp *= 1/2
    temp += location
    return temp    
# =============================================================================

def phase_estimation(signal,fe,w_length,ov_f):
#    Window size
    w = round(w_length*fe)
    
#    Hop step
    delta = round(w/ov_f)

#    Add zeros at the beginning and at the end to avoid side-effect
    s = np.concatenate((np.zeros(w),signal,np.zeros(w)))
        
#    Number of points for the signal
    N_pts = np.size(s,0)
    
#    Number of points for the FFT (with 0-padding)
    N_fft = int(pow(2,nextpow2(w)+5))
    
#    Maximal position for the moving window
    k_max = int(np.floor((N_pts-w)/delta))

#    Initialize output
    phase = np.zeros(k_max+1)
    phase_peak = np.zeros(k_max+1)
    
#    Phase estimated using the measured frequency
    phi_t = 0    
    
    for k in range(k_max+1):
        s_temp = s[(k*delta):((k*delta)+w)]
        
#       Apply window
        s_temp = s_temp*np.hanning(w)    
        
#       Compute FFT    
        s_temp_fft = np.fft.fft(s_temp,N_fft)
        s_temp_mod = np.abs(s_temp_fft)
        s_temp_arg = np.angle(s_temp_fft)
        
#        Estimate the main frequency (quadratic interpolation)
        i = np.argmax(s_temp_mod[:int(N_fft/2)])
        p = s_temp_mod[i-1]
        q = s_temp_mod[i]
        r = s_temp_mod[i+1]
        i0 = quadratic_argmax(i,p,q,r)
        f0 = i0*fe/N_fft

#        Store the phase of the peak in magnitude
        phase_peak[k] = s_temp_arg[i]   # Todo: interpolate       
        
#        Unwrap phase
        phi_t = phase[k-1] + (2*np.pi*f0*w_length/ov_f)
        if (k >= 1): phase[k] = phi_t + princ_arg(phase_peak[k]-phi_t)
#        if (k > 1): phase[k] = phase[k-1] + (2*np.pi*f0*w_length/ov_f)
        else: phase[k] = 0
    
#    Remove borders
#    s_processed = s_processed[w:(N_pts-w)]/1.5 
    
    return phase
# =============================================================================    