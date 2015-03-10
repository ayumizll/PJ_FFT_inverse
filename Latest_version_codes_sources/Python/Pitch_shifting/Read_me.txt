Execute main.py file to simulate the effect of pitch shifting.

The main function used to do this process is sine_model() in file sine_model.py

functions:

1. sine_model(s,a_win,overlap_factor,N_fft_a,N_fft_s,thres,alpha)

arguments:
# - s           = signal to analyse
# - a_win       = analysis window
# - N_fft_s     = synthesis FFT size
# - N_fft_a     = analysis FFT size
# - thres       = threshold for the peak detection (dB)
# - alpha       = shift factor (2^(N_semitones/12.0))