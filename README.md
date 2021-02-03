# kernel-density-estimation
Codes to optimally interpolate/extrapolate continuous density functions (1-d) from data

These matlab codes are used to create figures in the paper "Nonparametric, data-based kernel interpolation for particle-tracking simulations and kernel density estimation" by D. A. Benson, D. Bolster, S. Pankavich, and M. J. Schmidt (submitted).

The main code to run in matlab is called "BTC_smoother_DB_final.m".  This main program calls the subroutine "FT_h.m" to estimate the global bandwidth h_0 from data. Newly added (Feb. 2021) is an option to calculate the plug-in value of h_0 from Engel et al., 1994 by calling function "h_plug_in.m".  The main code first calculates the kernel density estimate (KDE) using Pedretti and Fernandez Garcia (2013) methods.  The main program then iteratively interpolates the density using kernels that are scaled copies of the interpolated density itself. Details are in the preprint by Benson et al. (2020) in ArXiv.  There are a number of options within the code to generate data from many known distributions and do the interpolations using either Gaussian kernels or the new data-based kernels.  One may also choose to read a file with random data.  One such file is provided from the paper called "BTC_50000", which are random particle arrival times. Any questions about options in the code may be submitted to dbenson@mines.edu.

Also included are "Fourier_test_eps.m", which estimates h_0 based on a Fourier transform of data.  This is done for an ensemble of data replicates and user-defined distributions.  The "h_v_n.m" code calculates the h_0 for various values of n from 10 to 10^6 and different assumptions discussed in the paper and plots results.

References
Benson, David A., Diogo Bolster, Stephen Pankavich, Michael J Schmidt (2020), Nonparametric, data-based kernel interpolation for particle-tracking simulations and kernel density estimation, ArXiv preprint, https://arxiv.org/abs/2010.06737.

Engel, J., Herrmann, E., & Gasser, T. (1994). An iterative bandwidth selector for ker- nel estimation of densities and their derivatives. Journal of Nonparametric Statistics, 4, 21–34. URL: https://doi.org/10.1080/10485259408832598. doi:10.1080/10485259408832598.

Pedretti, D., & Fernandez-Garcia, D. (2013). An automatic locally-adaptive method to estimate heavily-tailed breakthrough curves from particle distributions. Advances in Water Resources, 59, 52 – 65. URL: http://www.sciencedirect.com/science/article/pii/S0309170813000869. doi:https://doi.org/10.1016/j.advwatres.2013.05.006.
