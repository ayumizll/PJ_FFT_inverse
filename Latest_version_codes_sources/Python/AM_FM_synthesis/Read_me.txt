Script main_test.py perform the simulation of synthesis by modulated sinusoid model.
This simulation did not synthesize the main lobe of each partial.
This script is used to verify the performance of the estimator.
The estimator function is in the file estimation_fun.py

function:
1. sinemodel(x,Fs,w,N,Ns,H_ratio,th)  

arguments:

##############################################################################
#   parameters estimator of modulated sinusoidal model (derivative method)
#   reference: GENERALIZATION OF THE DERIVATIVE ANALYSIS METHOD TO NON-STATIONARY SINUSOIDAL MODELING
#              Conference on Digital Audio Effects (DAFx-08), Espoo, Finland, September 1-4, 2008
#   input:
#       s: input signal
#       Fs: sampling frequency
#       w: analysis window
#       N_fft: fft size
#       th: threshold for peaks detection, in negative dB
#   output:
#       a1:amplitude
#       omega1:frequency in radian
#       phi1:initial phase
#       mu1:modulation of amplitude in Neper/s
#       psi1:modulation of frequency in Hz/s 
#############################################################################