# =============================================================================
# Generate_sinus function
# @author: Quentin Biache, Emilie Abia, Lili ZHENG
# Version 3.0


# Arguments:
# - A0_range: list containing the lower and the upper value for the amplitude
# - f0_range: list containing the lower and the upper value for the frequency
# - phi0_range: list containing the lower and the upper value for the phase
# - duration: length (in seconds) of the sine waves
# - N: number of sinusoids to synthesize
# - fe: sampling frequency    
# =============================================================================

# Import NumPy library
import numpy as np

# Import wave library
from scipy.io import wavfile

# Import Random number generator library
import random


def generate_sin(A0_range,f0_range,phi0_range,duration,N,fe):
# -----------------------------------------------------------------------------
#    Generate the data matrix
# -----------------------------------------------------------------------------

#    Number of points required
    N_pts = int(round(duration*fe))

#   Ouput matrix
    M = np.zeros((N_pts,N+1))
    
#    Time, Frequency, Amplitude and Phase matrices    
    T = np.zeros((N_pts,N))
    A = np.zeros((N_pts,N)) 
    F = np.zeros((N_pts,N))
    Phi = np.zeros((N_pts,N))
        
#   Time vector
    t = np.linspace(0,duration,round(duration*fe))
    T[:] = t[:,np.newaxis]

#    Generate random amplitudes
    a = np.random.rand(1,N)*(A0_range[1]-A0_range[0])+A0_range[0]
    A[:] = a[np.newaxis,:]

#    Generate random frequencies
    f = np.random.rand(1,N)*(f0_range[1]-f0_range[0])+f0_range[0]
    F[:] = f[np.newaxis,:]
    
#    Generate random initial phase
    phi = np.random.rand(1,N)*(phi0_range[1]-phi0_range[0])+phi0_range[0]  
    Phi[:] = phi[np.newaxis,:]
    
#    Phase matrix
    phase = (2*np.pi*F*T)+Phi
    
#    Generate data matrix
    M[:,1:] = a*np.sin(phase)
    M[:,0] = t

# -----------------------------------------------------------------------------
#    Store each signal in a file
# -----------------------------------------------------------------------------
    for i in range(N):
        wavfile.write('Wave_' + str(i+1) + '.wav',fe,np.asarray(M[:,i+1]*32767,dtype=np.int16))


# -----------------------------------------------------------------------------
#    Create a table with the caracteristics of each file
# -----------------------------------------------------------------------------
    file_table = open('Table.txt', 'w')
    file_table.write('File number \tAmplitude \t\tFrequency \t\tPhase \n')
    
    for i in range(N):
        file_table.write(str(i+1)+'\t\t'+str(a[0,i])+'\t\t'+str(f[0,i])+'\t\t'+str(phi[0,i])+'\n')                
    
    file_table.close()

    


    return M   