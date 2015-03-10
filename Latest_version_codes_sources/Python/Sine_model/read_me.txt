Execute mian.py file to simulate the fft-1 synthesis (using stationary sinusoid model) which contain analyse and synthesis steps.
The main function sine_model() used to do this process is in file sine_model.py.

functions:

1.sine_model(s,a_win,overlap_factor,N_fft_a,N_fft_s,thres,robot)

# arguments:
# - s           = signal to analyse
# - a_win       = analysis window
# - N_fft_s     = synthesis FFT size
# - N_fft_a     = analysis FFT size
# - thres       = threshold for the peak detection (dB)
# - robot       = Activate robotic voice