Execute test_mian.py file to simulate the generation of sine signals with synthesized mainlbe.
Once 3 parameters of sines are known, function sine_generator() in file sinemodel.py will generate sine corresponding to these parameters.

functions:

1. sine_generator(a0,f0,phi0,L,Fs,w,N,Ns,H_ratio)

arguments:
#============================================================================
#   generate sinusoidal with model a0*exp(1j*(2.0*pi*f0*t+phi0)).real
#   a0: amplitude of sine
#   f0: frequency of sine
#   phi0: phase of sine
#   L: signal size
#   w: analysis window (odd size), N: FFT size,
#   t: threshold in negative dB, y: output sound
#   Ns: FFT size for synthesis (even)
#   H_ratio: Ratio (hop size)/(window size)
#============================================================================

