# -*- coding: utf-8 -*-
# =============================================================================
# Script to call test2
# @author: Quentin
# Version 1.0
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
from phase_estimation import phase_estimation as phase_estimation

# =============================================================================

# Settings
filename = 'Wave_1'    # Name without extension
window_length = 0.023   # Window's length (in ms)
overlap_factor = 4      # Overlap_factor = window_size/hop_step

# Import data
[fe, s_int16] = wavfile.read(filename + '.wav')

# Convert in floating point data
s = s_int16/32767.0

# We work only with the first component (mono)
if np.ndim(s) >= 2: s = s[:,0]
N_pts = np.size(s,0)

N_pts = 88200
t_real = np.array(range(N_pts))/float(fe)
s = np.sin((2*np.pi*440*t_real)+(np.pi*500*t_real*t_real))

# Call the function
phase = phase_estimation(s,fe,window_length,overlap_factor)
# Create time vector
#t = np.array(range(N_pts))/float(fe)
t = (np.array(range(np.size(phase,0)))*window_length/overlap_factor) + (window_length/2)

#t_real = np.array(range(N_pts))/float(fe)
plt.figure()
plt.plot(t,phase)
plt.grid()
plt.hold(True)
plt.plot(t_real,(2*np.pi*440*t_real)+(np.pi*500*t_real*t_real),'r')
plt.xlabel('Time (s)')
plt.ylabel('Phase(t)')