Long project 2015
Quentin Biache - Lili Zheng - Emilie Abia

The codes in this folder were developped to perform the lobe synthesis using
a two-sines approximation.

The script "AM_FM_animation" plots a sampled main lobe with either amplitude
or frequency modulation. Each point of the lobe is stored in a .mat file.

The coefficients for the two-sines model can be computed using sine_interp.
Warning, "AM_FM_animation" must be run first. 

The function main_lobe_FM reconstructs a lobe using the coefficients for the
sine model.

Finally, "AM_FM_Reconstruction" reconstructs a signal in the temporal 
domain using the model developped previously.
