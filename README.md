## Uniform-in-time rational approximation of the matrix exponential with real poles

These are the MATLAB codes accompanying the paper of the above title. Some of these files have additional requirements: 

* Advanpix Multiprecision Toolbox available at https://www.advanpix.com/ needs to be in the MATLAB path
* MATLAB Rational Krylov Toolbox available at https://guettel.com/rktoolbox/ needs to be in the MATLAB path
* The matrices $C,M,q$ available from https://sparse.tamu.edu/Guettel/TEM181302 need to be in the same folder

The figures in the paper are reproduced as follows:

* Figure 1: Run `drive_weight_vs_unweighted.m`
* Figure 2: Run `drive_contour_plot.m`
* Figure 3: Run `drive_three_approximants.m`
* Figure 4 (subplots `x=1,...,4`): Run `drive_time_uniform_error_x.m`
* Figure 5 (subplots `x=1,...,3`): Run `drive_illustrate_error_analysis_x.m`
* Figure 6: Run `drive_timing_comparison.m`
* Figure 7: Run `drive_tem.m`
