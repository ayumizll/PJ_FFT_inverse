Execute test_main.py to simulate the estimation of 5 parameters
of modulated sinus model such as a*exp(mu*t+1j*(phi+omega*t+0.5*psi*t^2)).
The estimator function named estimator1() is in the file estimation_derivative.py.

functions:
1.estimator1(s,Fs,w,N_fft,th)

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