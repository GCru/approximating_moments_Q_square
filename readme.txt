
1. For Figure 1, 'use radius_of_convergence_sqrt'

2. 'imhoff.py' iis the imhoff routine translated from R  eith the help of ChatGPT

3. Figure 2 uses 'draw_linear_comb_chi_sqr' with Inhoff function

4. Use 'error_histograms_monte_carlo_vs_taylor' to complete Tables 1,2,3. It it run seperately for n=2,...10,100. It
    creates an intermediate 'monte_carlo_list_n.npy' file which gives eigenvalue set and expected value and variance.
    This is then used to creaate csv files 'taylor_errors_n' with all the info for required for the tables as well as
    a file 'taylor_max_eigenvalues_n' which has the eigenvalue sets where the errors for expected value and variance
    are the largest.

5.  To draw Figures 3 and 4 use 'histograms_of_Taylor' which reads data from the numpy files of point 4.

6. The script 'minimum_of_eigenvalue_expressions' double chacks Appendix E

7.  The script 'check_Taylor_approximations', 'examine_cov_matrix_behavious' and 'test_marchenko_pastur', 'linear_combinations_of_chi_square',
    'compare_chi_and_gamma', 'compare_sum_of_chi_squares' were used
    to check and aid understanding.

8.  To check which script is the quickest for drawing multivariate independent normal I used 'speed_check_multivariate_normal'