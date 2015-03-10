Long project 2015
Quentin Biache - Lili Zheng - Emilie Abia

The codes in this folder were developped to perform the lobe synthesis
using interpolation.

First of all, one must launch 'LUT_generator' to create a .mat file 
containing the reference lobes.

Then, an error plot can be performed in the (ACR,FCR) plane using 
error_plot.m

One can also generate chirps with custom parameters and compare them to the
reconstructed one using Synthesis_comparison.m

Finally, the error plot can be exploited to generate a non-uniform grid of
control points. For this, launch Non_uniform grid. The updated error plot
can be done simply be executing again Error_plot.m