# -*- coding: utf-8 -*-
# =============================================================================
# Sine model analysis and synthesis with transposition
# @author: Quentin
# Version 4.0
# =============================================================================





# =============================================================================
# Import section

# Import NumPy library
import numpy as np

# Import Signal for windows
import scipy.signal as sig

import matplotlib.pyplot as plt
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
        
        
# Function to estimate the maximum using a quadratic interpolation
# u = f(p-1)
# v = f(p)
# w = f(p+1)
# quadratic_max(u,v,w) = max f(x) for x in [p-1,p+1]
def quadratic_max(u,v,w):    
    temp = pow(u,2) + (16.0*pow(v,2)) + pow(w,2)
    temp += (-8*u*v) + (-8*v*w) + (-2*u*w)
    temp = temp/(u-(2*v)+w)
    temp *= (-1./8.)
    return temp  

# Function to estimate the argmax using a quadratic interpolation
# u = f(p-1)
# v = f(p)
# w = f(p+1)
# quadratic_argmax(location,u,v,w) = argmax f(x) for x in [location-1,location+1]
def quadratic_argmax(location,u,v,w):    
    temp = u-w
    temp = temp/(u-(2*v)+w)
    temp *= 1./2.
    temp += location
    return temp          

# Dirichlet kernel function
def D_kernel(x,N_w):
    y = np.sin(N_w*x/2)/np.sin(x/2)

#    Regularize value at zero    
    y[np.isnan(y)] = float(N_w)
    return y

# Function to synthesize the FFT of Blackman-Harris window
def blackman_harris_92db(x,N_w,N_fft):
    
#    Convert index into angular velocity
    omega = 2*np.pi*x/N_fft

#    Angular velocity step
    d_omega = 2*np.pi/N_fft
    
    y = np.zeros((np.size(x)))    
    
#    Define constants of the Blackman window    
    a = np.array([0.35875, 0.48829, 0.14128, 0.01168])
    
#   Sum the frequency-shifted kernels    
    for k in range(4):
        y += 0.5*a[k]*(D_kernel(omega-(k*d_omega*N_fft/N_w),N_w)+D_kernel(omega+(k*d_omega*N_fft/N_w),N_w))
    
#    Normalize values    
    y = y/(N_w*a[0])    
    
    return y

def sine_synthesis(p_loc,p_mag,p_phase,N_fft,N_w):
    
#    Initialize output spectrum
    Y = np.zeros((N_fft),dtype = complex)

#    Number of points in the positive spectrum range (f = 0.5 is also included)
    N_h = 1+int(0.5*N_fft)

#    Iterate on the number of peaks    
    for k in range(np.size(p_loc)):

#           Current pic location                        
            loc = p_loc[k]                        
                        
            if ((loc <= 1.0) or (loc >= N_h-1.0)):
                continue

#           Center peak location at 0
            zero_loc = round(loc)-loc

#           Points where the main lobe must be evaluated 
            x = zero_loc + np.array(range(-4,5))
            
#           Generate lobe with the required magnitude           
            l_mag = blackman_harris_92db(x,N_fft,N_w)*(10**((p_mag[k])/20))
            
#           Frequency index which have to be filled (we take the closest one)
            f_index = round(loc) + np.array(range(-4,5))

#           Add the whole main lobe to the spectrum
            for l in np.array(range(9)):
#               Change the phase if there are negative components in the lobe
#               (Hermitian symmetry: X(-f) = X*(+f)) 
                if (f_index[l] < 0):
                    Y[-f_index[l]] += l_mag[l]*np.exp(-1.0j*p_phase[k])
                        
                elif (f_index[l] > (N_h-1.0)):
                    Y[N_fft-f_index[l]] += l_mag[l]*np.exp(-1.0j*p_phase[k])
                    
                else:                    
                    Y[f_index[l]] += l_mag[l]*np.exp(1.0j*p_phase[k])
                    Y[f_index[l]] += l_mag[l]*np.exp(-1.0j*p_phase[k])*((f_index[l] == 0.0) or (f_index[l] == (N_h-1.0)))
            
            Y[N_h:] = np.conj(Y[N_h-2:0:-1])            
            
    return Y
# =============================================================================



# =============================================================================
# Main code

# Function that performs the analysis and synthesis of the signal
# (No processing is performed up to now) using the additive sine model
# Requiered arguments:
# - s           = signal to analyse
# - a_win       = analysis window
# - N_fft_s     = synthesis FFT size
# - N_fft_a     = analysis FFT size
# - thres       = threshold for the peak detection (dB)

def sine_model(s,a_win,overlap_factor,N_fft_a,N_fft_s,thres,alpha):

#   Window size (analysis)
    N_a = np.size(a_win)

#   Window size (synthesis). We keep N_fft_synthesis and N_s identical, since
#   0-padding will not add any relevant information
    N_s = N_fft_s
    mid_fft_a = int((N_fft_a+1)/2)
    
#   Normalize analysis window
    a_win = a_win/np.sum(a_win)    

#   Hop size
    H = int(round(N_a/overlap_factor))
       
#   Total number of signal points
    N_pts = np.size(s,0)
    
#   Initialize output
    y = np.zeros(N_pts)
    
#   Build corrected Blackman-Harris window (synthesis)    
    s_win = np.zeros(N_fft_s)
    t_win = sig.triang((2*H)-1)

#   Add triangular window at the center    
    s_win[(N_s/2)+1-H+1:(N_s/2)+H+1] = t_win

    b_win = sig.blackmanharris(N_s)
    b_win = b_win/np.sum(b_win)    
    
    s_win[(N_s/2)+1-H+1:(N_s/2)+H+1] = s_win[(N_s/2)+1-H+1:(N_s/2)+H+1]/b_win[(N_s/2)+1-H+1:(N_s/2)+H+1]
        
#   Initialize the reading pointer
    p_in = int(max(0.5*(N_a-1),(0.5*N_s)-1))
    
#   Initialize the upper limit of the pointer
    p_end = N_pts-1-int(max(0.5*N_s,0.5*(N_a-1)))
    
#   Debug variable    
    k = 0
    debug_tab = np.zeros(10000)
    debug_tab_2 = np.zeros(10000)    
    print('H = ' + str(H))
    print('N_fft_a = ' + str(N_fft_a))
    
#   ---------------------------------------------------------------------------

#   Phase increment vector
    delta_phi = 2.0*(alpha-1.0)*np.pi*np.arange(mid_fft_a+1)*H/N_fft_a

    
#   ---------------------------------------------------------------------------

    while p_in < p_end:

#       Fill the buffer
        s_w = s[p_in-int(0.5*(N_a-1)):p_in+int(0.5*(N_a-1))+1]
        
#       Apply analysis window
        s_w = s_w*a_win         

#       Zero phase correction and 0-padding
        s_temp = np.zeros(N_fft_a)
        s_temp[0:((N_a-1)/2)] = s_w[((N_a+1)/2):]
        s_temp[N_fft_a-((N_a+1)/2):] = s_w[:((N_a+1)/2)]
        
#       Compute FFT     
        X = np.fft.fft(s_temp,N_fft_a)
        mod_X = 20*np.log10(np.abs(X[:mid_fft_a+1]))
        arg_X = (np.angle(X[:mid_fft_a+1]))
        
        lower_bound     = mod_X[1:mid_fft_a-1] > mod_X[0:mid_fft_a-2]
        max_estimation  = mod_X[1:mid_fft_a-1] > thres
        upper_bound     = mod_X[1:mid_fft_a-1] > mod_X[2:mid_fft_a]

#       Search for the indexes where the three conditions are satisfied
        peak_loc = np.where(lower_bound & max_estimation & upper_bound)        
        peak_loc = peak_loc[0] + 1                     
        
#       Interpolated values
        peak_mag = quadratic_max(mod_X[peak_loc-1],mod_X[peak_loc],mod_X[peak_loc+1])
        peak_loc_interp = quadratic_argmax(peak_loc,mod_X[peak_loc-1],mod_X[peak_loc],mod_X[peak_loc+1])
        
        if (k == 0):
            arg_X_shifted = np.unwrap(arg_X)
            arg_X_shifted_cumulated = np.zeros(np.size(arg_X))
#            arg_X_shifted = np.zeros(np.size(arg_X))
        else:
            arg_X_shifted_cumulated += delta_phi
            arg_X_shifted_cumulated = arg_X_shifted_cumulated % (2*np.pi)
            arg_X_shifted = arg_X
            arg_X_shifted += arg_X_shifted_cumulated            
            arg_X_shifted = np.unwrap(arg_X_shifted)
            
            

        peak_phase = np.interp(peak_loc_interp,np.arange(np.size(arg_X_shifted)),arg_X_shifted)

#        arg_X_shifted = arg_X_shifted % (2*np.pi)
#        debug_tab[k] = arg_X[10]
#        if (k == 0):
#            debug_tab_2[k] = arg_X[peak_loc[0]]
#        else:
##            debug_tab_2[k] = debug_tab[k]-debug_tab[k-1]
#            debug_tab_2[k] = arg_X[peak_loc[0]]
#        print[peak_loc]
#        print[peak_loc_interp]
#        print(np.size(arg_X_shifted))
#        print(np.size(arg_X))
#        print('Position du lobe initiale = ' + str(peak_loc))
#        print('Position du lobe interpolée = ' + str(peak_loc_interp))
#        print('Fréquence interpolée = ' + str(peak_loc_interp*44100.0/N_fft_a))
##        print(arg_X[peak_loc[0]])
#
#        print('Phase du lobe avant = ' + str(arg_X_shifted[peak_loc]))
#        print('Phase du lobe interpolée = ' + str(peak_phase))
#        print('Phase du lobe après = ' + str(arg_X_shifted[peak_loc+1]))        
##        print((delta_phi[peak_loc[0]]+0*np.pi)/(2*np.pi*0.023/4.))
#        print('------------------')


        
#       Fun ...(0 phase reconstruction, aka robotic voice)
#        peak_phase = np.zeros(np.size(arg_X))
        
#       Synthesize signal and remove zero-phase correction
        Y = sine_synthesis(peak_loc_interp*alpha*N_fft_s/N_fft_a,peak_mag,peak_phase,N_fft_s,N_fft_s)
        y_w = np.fft.fftshift(np.real(np.fft.ifft(Y)))

#       -----------------------------------------------------------------------
#       Section for debugging purposes
        if (k == 9): 
            plt.plot(s_w,'b')
            plt.title('Input signal versus output signal (time domain)')            
            plt.hold(True)
            plt.plot(y_w,'r')
            plt.grid()
            
            plt.figure()
            plt.plot(20*np.log10(X),'b')
            plt.title('Input signal versus output signal (frequency domain)')
            plt.hold(True)
            plt.plot(20*np.log10(Y),'r')
            plt.grid()

            plt.figure()
            plt.plot(arg_X_shifted,'b')
            plt.title('Input signal versus output signal (frequency domain)')
            plt.grid()                
            
        k += 1            
#       -----------------------------------------------------------------------

#       Apply synthesis window           
        y_w = y_w*s_win

#       Overlap-add                
        y[p_in-int(0.5*N_s)+1:p_in+int(0.5*N_s)+1] += y_w

#       Jump to the next position
        p_in += H
        
#    plt.figure()
#    plt.plot(np.arange(10000),lobe_phase,'r',np.arange(10000),lobe_phase_2,'b')        
#
        
        
        
    # -----------------------------------------------------------------------------
#    Create a table with the caracteristics of each file
# -----------------------------------------------------------------------------
    file_table = open('Data.txt', 'w')
    
    for i in range(10000):
        file_table.write(str(debug_tab_2[i])+'\n')                
    
    file_table.close()        
    
    
    plt.figure()
    plt.plot(debug_tab_2)
    plt.title('Debug_tab content')
    
    return y
# =============================================================================    