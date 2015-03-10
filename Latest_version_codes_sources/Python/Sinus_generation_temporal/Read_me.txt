Execute main.py file to generate wave files containing stationnary sine waves with random parameters.

File generate_sinus.py contain the necessaries functions such as generate_sin().

functions:
1. generate_sin(A0_range,f0_range,phi0_range,duration,N,fe)

# Arguments:
# - A0_range: list containing the lower and the upper value for the amplitude
# - f0_range: list containing the lower and the upper value for the frequency
# - phi0_range: list containing the lower and the upper value for the phase
# - duration: length (in seconds) of the sine waves
# - N: number of sinusoids to synthesize
# - fe: sampling frequency 