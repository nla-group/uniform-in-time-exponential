## Uniform-in-time rational approximation of the matrix exponential with real poles

These are the MATLAB codes accompanying the paper of the above title. Some of these files have additional requirements: 

* Advanpix Multiprecision Toolbox available at https://www.advanpix.com/ needs to be in the MATLAB path
* MATLAB Rational Krylov Toolbox available at https://guettel.com/rktoolbox/ needs to be in the MATLAB path
* The matrices $C,M,q$ available from https://sparse.tamu.edu/Guettel/TEM181302 need to be in the same folder

The figures in the paper are reproduced by the following scripts:

* Figure 1: `drive_weight_vs_unweighted.m`
* Figure 2: `drive_contour_plot.m`
* Figure 3: `drive_three_approximants.m`
* Figure 4 (subplots `x=1,...,4`): `drive_time_uniform_error_x.m`
* Figure 5 (subplots `x=1,...,3`): `drive_illustrate_error_analysis_x.m`
* Figure 6: `drive_timing_comparison.m`
* Figure 7: `drive_tem.m`

In addition, we provide the script `myaddpath.m`, which should be modified to include the paths with the installation directory of your RKFIT Toolbox and Advanpix Multiprecision Computing Toolbox. `myaddpath.m` is called by all of the above scripts.
