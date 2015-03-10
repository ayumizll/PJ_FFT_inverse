# -*- coding: utf-8 -*-
# =============================================================================
# Study the main lobe shape in case of FM modulations
# @author: Quentin
# Version 2.0
# This function plots the FFT of sine waves with different FM parameters
# Minor modifications
# =============================================================================

# =============================================================================
# Import plot library
import matplotlib.pyplot as plt
plt.close('all')

# Import NumPy library
import numpy as np

import scipy.signal as sig
# =============================================================================



# =============================================================================
# Functions section


# Function to find the closest power of 2 for a given real number.
# This is useful to add zero-padding when computing the FFT
def nextpow2(i):
    if (i == 0): return -float ('inf') 
    elif (i <= 1): return 0.0
    else:  
        n = 1
        while n < i: n *= 2
        return np.log(n)/np.log(2)
# =============================================================================


# -----------------------------------------------------------------------------
# Settings

# Window length in s
window_length = 0.046

# Sampling frequency
f_e = 44100
  
# Modulation time
t_mod = window_length

# Modulation rate parameters
q_min = -20
q_max = 20 

# Number of plots
N_plots = 11

# Initial frequency
f_0 = 1000

# 0-padding constant
padding = 5

# -----------------------------------------------------------------------------

# Create time vector
N_pts = int(round(window_length*f_e))
if ((N_pts % 2) == 0): N_pts += 1
t = np.array(range(N_pts))/float(f_e)

# FFT parameters
a_win = sig.blackmanharris(N_pts)
N_fft_a = int(2**(nextpow2(N_pts)+padding))

a_win = a_win/np.sum(a_win)

q = np.linspace(q_min,q_max,N_plots)

array_max = np.zeros(np.size(q))
array_argmax = np.zeros(np.size(q))


for k in range(np.size(q)):
    
    # Create signal
    phase = 2*np.pi*f_0*(t + (0.5*0.01*q[k]*(t**2)/t_mod))
    s = np.sin(phase)
    s_w = s*a_win
    
    # Compute FFT
    fft_s = np.fft.fft(s_w,N_fft_a)
    f = f_e*np.array(range(N_fft_a))/float(N_fft_a)
    
    array_max[k] = np.max(20*np.log10(abs(fft_s[:N_fft_a/2])))
    array_argmax[k] = np.argmax(abs(fft_s[:N_fft_a/2]))
    
    zoom = np.where((f > 700) & (f < 1300))
      
    # Plot some results
    plt.figure('Results - 1')
    plt.subplot(2,2,1)
    plt.plot(f[zoom],20*np.log10(abs(fft_s[zoom])))
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('FFT Magnitude (dB)')
    plt.title('Main lobe shape for different linear FM ('+str(q_min)+'% to + '+str(q_max)+'%)')
    plt.grid(True)
    plt.hold(True)
    plt.plot(array_argmax[k]*f_e/N_fft_a,array_max[k],'rx')
    plt.subplot(2,2,3)
    plt.plot(f[zoom],np.unwrap(np.angle(fft_s[zoom])*5.)/5.)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('FFT Phase')    
    plt.grid(True)


print('FFT magnitudes maximal values = ' + str(10**(array_max/20)))

plt.subplot(2,2,2)
plt.plot(q/t_mod,array_max)
plt.xlabel('Modulation rate (% of f_0 per seconds)')
plt.ylabel('FFT Magnitude peak value (dB)')
plt.title('Main lobe parameters for different linear FM ('+str(q_min)+'% to + '+str(q_max)+'%)')
plt.grid(True)
plt.subplot(2,2,4)
plt.plot(q/t_mod,(array_argmax*f_e/N_fft_a)-f_0)
plt.xlabel('Modulation rate (% of f_0 per seconds)')
plt.ylabel('Peak location relative to f_0 (Hz)')
plt.grid(True)

# 16-bits normalization and write the signal in a wave file
#wavfile.write('Chirp.wav',f_e,np.asarray(s*32767,dtype=np.int16))