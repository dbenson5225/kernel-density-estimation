# kernel-density-estimation
Codes to optimally interpolate/extrapolate continuous density functions (1-d) from data

These matlab codes are used to create figures in the paper "Nonparametric, data-based kernel interpolation for particle-tracking simulations and kernel density estimation" by D. A. Benson, D. Bolster, S. Pankavich, and M. J. Schmidt (submitted).

The main code to run in matlab is called "BTC_smoother_DB_2.m".  This main program calls the subroutine "FT_h.m" to estimate the global bandwidth from data.  The main code then iteratively interpolates the density using kernels that are scaled copies of the interploated density itself.  Details are in the preprint in ArXiv.  There are a number of options within the code to generate data from many known distributions and do the interpolations using either Gaussian kernels or the new data-based kernels.  One may also choose to read a file with random data.  One such file is provided from the paper called "BTC_50000", which are random particle arrival times. Any questions about options in the code may be submitted to dbenson@mines.edu.
