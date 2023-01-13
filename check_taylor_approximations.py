
import numpy, math
from datetime import date
from statistics import stdev, mean
from math import pi, log, sqrt
from scipy.stats import kstest, kurtosis, skew, skewtest, jarque_bera, wasserstein_distance, median_abs_deviation
from scipy.special import erf
from scipy import optimize

from momentchi2 import hbe

from bokeh.plotting import figure, output_file, show
from bokeh.models import Range1d, PrintfTickFormatter, NumeralTickFormatter, Title
from bokeh.layouts import row
from bokeh.models import ColumnDataSource

"""from fundskill_utilities.fundskill_utilities import change_to_another_month_end, \
	change_to_month_end

from latex_label_utility import LatexLabel

from fundskill_utilities.fundskill_shared_bokeh_utilities import setup_date_ticker_jump

from bokeh_constants import *"""


def create_covariance_matrix(sigma, rho, n):
	""" Create a covariance matrix
	:param sigma:
	:param rho:
	:param n:
	:return:
	"""
	V = numpy.ones((n, n)) * rho * sigma * sigma
	
	for idx in range(n):
		V[idx, idx] = sigma * sigma
	
	return V

def monte_carlo_moments_wVw_and_sqrt_wVw(V, iterations=1000):
	
	n=numpy.shape(V)[0]
	
	I = numpy.zeros((n, n))
	numpy.fill_diagonal(I, 1)
	wVw = []
	wVw_sqrt = []

	for count in range(iterations):
		
		w = numpy.random.default_rng().multivariate_normal(numpy.zeros(n), (1 / n) * I)
		
		wVw.append(numpy.dot(w, numpy.dot(V, w)))
		wVw_sqrt.append(wVw[-1] ** 0.5)
		
	mean_wVw_monte_carlo = numpy.mean(wVw)
	mean_sqrt_wVw_monte_carlo = numpy.mean(wVw_sqrt)
	
	var_wVw_monte_carlo = numpy.var(wVw)
	var_sqrt_wVw_monte_carlo =numpy.var(wVw_sqrt)
	
	taylor_support_min_monte_carlo = min(wVw)
	taylor_support_max_monte_carlo = max(wVw)
	return mean_wVw_monte_carlo, var_wVw_monte_carlo, mean_sqrt_wVw_monte_carlo, var_sqrt_wVw_monte_carlo, \
		taylor_support_min_monte_carlo, taylor_support_max_monte_carlo

def calculate_algebraic_moments_wVw(V):
	
	n = numpy.shape(V)[0]
	
	mean_wVw_algebraic = numpy.trace(V) / n
	
	var_wVw_algebraic = 2 * numpy.average(numpy.square(V))

	# or in term of rho and sigma
	# 2 * ((rho ** 2) * (n - 1) * n + (sigma ** 2) * n) / (n * n)
	
	return mean_wVw_algebraic, var_wVw_algebraic


def calculate_cumulant(k, eigenvalue_list ):
	
	n=len(eigenvalue_list)
	kappa = 2**(k-1)*math.factorial(k-1) / n**k

	sum_product = 0
	for eigenvalue in eigenvalue_list:
		sum_product = sum_product + eigenvalue**k
		
	return kappa*sum_product


def calculate_taylor_mean_sqrt_wVw(V):
	n = numpy.shape(V)[0]
	
	mean_sqrt_wVw_taylor = (numpy.trace(V) / n) ** 0.5 - (1.0 / 8.0) * (numpy.trace(V) / n) ** (
		-1.5) * 2 * numpy.average(
		numpy.square(V))
	
	return mean_sqrt_wVw_taylor


def calculate_taylor_2_var_sqrt_wVw(V):
	n = numpy.shape(V)[0]
	
	var_sqrt_wVw_taylor_2 = (numpy.linalg.norm(V)**2)/ (2*n*numpy.trace(V))
	
	return var_sqrt_wVw_taylor_2


def calculate_taylor_3_var_sqrt_vWv(V):
	
	eigenvalue_list, _ = numpy.linalg.eig(V)

	kappa_list=[]
	
	for idx in range(5):
		if idx == 0:
			kappa_list.append(0)
		else:
			kappa_list.append(calculate_cumulant(idx, eigenvalue_list))

	var_Q = kappa_list[2]
	
	var_Q_sqr = kappa_list[4]+4*kappa_list[3]*kappa_list[1] + 2*kappa_list[2]**2 + 4* kappa_list[2]*kappa_list[1]**2
	
	cov_Q_Q_sqr = kappa_list[3] + 2* kappa_list[2] *kappa_list[1]
	
	term1 = (9/16) * var_Q /kappa_list[1]
	
	term2 = (1/64)*var_Q_sqr/kappa_list[1]**3
	
	term3 = (3/16)*cov_Q_Q_sqr/kappa_list[1]**2
	
	return term1+term2-term3
	
	

def calculate_upper_bound_mean_sqrt_wVw(V):
	
	n = numpy.shape(V)[0]
	
	return   (numpy.trace(V) / n) ** 0.5


def calculate_taylor_suport_algebraic(V):
	
	n = numpy.shape(V)[0]
	return 0, 2 * numpy.trace(V)


if __name__ == '__main__':
	
	# examine how random funds weights work
	# see_how_random_fund_weights_work(100)
	# exit()
	
	# check if funds drawn from standard random normal gives normal return distributions, e.g. affine transformation
	# calculate_distribution_of_random_funds_returns()
	# exit()
	
	# check how marcenko pastur works
	# sigma =1
	# n=1000
	# T=1000
	# V= test_Marcenko_Pastur_covariance_matrix(sigma, n, T)
	# exit()
	
	# Set up covariance matrix of size n
	n = 20
	rho = 0# set to zero for all eigenvalues the same and set to one for only one positve eigenvalue
	sigma = 0.1**0.5
	V = create_covariance_matrix(sigma, rho, n)
	print(V)
	#V[1,1] =0.4
	# Calculate eigenvalues
	eigenvalue_list, _ = numpy.linalg.eig(V)
	print('Max eigenvalue:', max(eigenvalue_list), ' Min eigenvalue', min(eigenvalue_list))
	print(eigenvalue_list)
	
	# Calculate cumulant
	k=3
	print('Cumulant ',k, ': ', calculate_cumulant(k, eigenvalue_list))
	
	# Calculate mean and variance using monte carlo
	mean_wVw_monte_carlo, var_wVw_monte_carlo, mean_sqrt_wVw_monte_carlo, var_sqrt_wVw_monte_carlo,\
		taylor_support_min_monte_carlo, taylor_support_max_monte_carlo = monte_carlo_moments_wVw_and_sqrt_wVw(V)
	
	
	print(mean_wVw_monte_carlo, var_wVw_monte_carlo, mean_sqrt_wVw_monte_carlo, var_sqrt_wVw_monte_carlo)
	# Calculate algebraic mean and varince of wVw
	mean_wVw_algebraic, var_wVw_elgebraic = calculate_algebraic_moments_wVw(V)
	
	print(mean_wVw_algebraic, var_wVw_monte_carlo)
	
	mean_sqrt_wVw_taylor = calculate_taylor_mean_sqrt_wVw(V)
	print(mean_sqrt_wVw_taylor)
	
	print('Percentage error taylor vs monte carlo for expected value')
	print(100* (mean_sqrt_wVw_taylor-mean_sqrt_wVw_monte_carlo)/(mean_sqrt_wVw_monte_carlo),'%')
	
	
	
	var_sqrt_wVw_taylor_2 = calculate_taylor_2_var_sqrt_wVw((V))
	
	print(var_sqrt_wVw_taylor_2, var_sqrt_wVw_monte_carlo, numpy.trace(V)/(2*numpy.shape(V)[0]**2))
	
	print('Percentage error taylor 2 vs monte carlo for variance')
	print(100 * (var_sqrt_wVw_taylor_2 - var_sqrt_wVw_monte_carlo) / (var_sqrt_wVw_monte_carlo), '%')
	
	var_sqrt_wVw_taylor_3 = calculate_taylor_3_var_sqrt_vWv(V)
	print(var_sqrt_wVw_taylor_3)
	
	print('Percentage error taylor 2 vs monte carlo for variance')
	print(100 * (var_sqrt_wVw_taylor_3 - var_sqrt_wVw_monte_carlo) / (var_sqrt_wVw_monte_carlo), '%')
	
	print(numpy.trace(V)/(2*numpy.shape(V)[0]**2))
	exit()
	

	#draw_sum_of_chi_cdf(w)
	exit()