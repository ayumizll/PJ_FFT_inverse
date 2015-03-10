# -*- coding: utf-8 -*-
# =============================================================================
# Script to call the transposition  module
# @author: Quentin Biache, Emilie Abia, Lili ZHENG
# Version 2.0
# =============================================================================

# =============================================================================
# Import plot library
import matplotlib.pyplot as plt
plt.close('all')

# Import NumPy library
import numpy as np

# Import wave library
from scipy.io import wavfile

# Import function
from sine_model import *

import scipy.signal as sig


# =============================================================================

# Settings
filename = 'Sons/Lili'    # Name without extension
window_length = 0.023   # Window's length (in ms)
overlap_factor = 4.0      # Overlap_factor = window_size/hop_step
N_semitones = 3.0

# Import data
[f_e, s_int16] = wavfile.read(filename + '.wav')
#s_int16 = s_int16[:,0]

# Convert in floating point data
s = s_int16/32767.0

# We work only with the first component (mono)
if np.ndim(s) >= 2: s = s[:,0]
N_pts = np.size(s,0)

# Create time vector
t = np.array(range(N_pts))/float(f_e)

N_w = int(round(window_length*f_e))
if ((N_w % 2) == 0): N_w += 1

a_win = sig.blackmanharris(N_w)

N_fft_a = int(2**(nextpow2(N_w)))

# Call the function
s_processed = sine_model(s,a_win,overlap_factor,N_fft_a,1024,-80,2**(N_semitones/12.0))

# 16-bits normalization and write the signal in a wave file
wavfile.write(filename + '_processed.wav',f_e,np.asarray(s_processed*32767,dtype=np.int16))

  
# Plot some results
plt.figure('Results - 1')
plt.plot(t,s,t,s_processed)
plt.xlabel('Time (s)')
plt.ylabel('Signal')
plt.title('Comparison - original and resynthesized signal')
plt.legend(['Original signal','Processed signal'])
plt.grid()
