# =============================================================================
# Test of the generate_sinus function
# @author: Quentin Biache, Emilie Abia, Lili ZHENG
# Version 3.0
# =============================================================================

# =============================================================================
# Import plot library
import matplotlib.pyplot as plt

# Import NumPy library
import numpy as np

# Import function
from generate_sinus import generate_sin as generate_sin

# =============================================================================

# Settings
A0_range    = [0.1,0.9]
f0_range    = [100.,1000.]
phi0_range  = [-np.pi,np.pi]
duration    = 2
N           = 10
fe          = 44100.

# Call the function
M = generate_sin(A0_range,f0_range,phi0_range,duration,N,fe)

# Plot one of the sine wave
plt.plot(M[:1001,0],M[:1001,1])
plt.xlabel('Time (s)')
plt.ylabel('Signal')
plt.title('First generated sine wave')
plt.grid()