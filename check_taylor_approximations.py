
import numpy
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
	
	w, v = numpy.linalg.eig(V)
	print('Max eigenvalue:', max(w), ' Min eigenvalue', min(w))
	# print(w)
	# print(w)
	
	return V


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
	n = 100
	rho = 1.0
	sigma = 0.1
	V = create_covariance_matrix(sigma, rho, n)
	# V[1,1] =0.4
	
	I = numpy.zeros((n, n))
	numpy.fill_diagonal(I, 1)
	
	wVw = []
	wVw_sqrt = []
	adjusted_wVw = []
	fund_sigma_if_zero_sum = []
	for count in range(1000):
		w = numpy.random.default_rng().multivariate_normal(numpy.zeros(n), (1 / n) * I)
		
		wVw.append(numpy.dot(w, numpy.dot(V, w)))
		wVw_sqrt.append(wVw[-1] ** 0.5)
		
		adjusted_w = w + (0 - numpy.sum(w) / n)
		adjusted_wVw.append(numpy.dot(adjusted_w, numpy.dot(V, adjusted_w)))
		
		# print(numpy.sum(w), wVw[-1], numpy.sum(adjusted_w), adjusted_wVw[-1])
		
		adjusted_A_euclid = numpy.dot(adjusted_w, adjusted_w) ** 0.5
		fund_sigma_if_zero_sum.append((adjusted_A_euclid * sigma * (1 - rho) ** 0.5) ** 2)
	
	print()
	print('Mean of tracking error variance. Actual ', numpy.mean(wVw), 'Theoretical', numpy.trace(V) / n)
	
	# print('Mean when weights are adjusted', numpy.mean(adjusted_wVw), 'Theoretical mean',  numpy.mean(fund_sigma_if_zero_sum) )
	
	print('Variance of tracking error variance. Actual:', numpy.var(wVw), 'Theoretical',
		  2 * ((rho ** 2) * (n - 1) * n + (sigma ** 2) * n) / (n * n), 'or', 2 * numpy.average(numpy.square(V)))
	
	print()
	print('Expected value (mean) of tracking error. Actual: ', numpy.mean(wVw_sqrt), 'Theoretical upper bound:',
		  (numpy.trace(V) / n) ** 0.5)
	
	# print('Check by calculating mean tracking error variance from formula', numpy.mean(wVw_sqrt)**2+ numpy.var(wVw_sqrt) )
	
	taylor = (numpy.trace(V) / n) ** 0.5 - (1.0 / 8.0) * (numpy.trace(V) / n) ** (-1.5) * 2 * numpy.average(
		numpy.square(V))
	print('Expected value (mean) of tracking error. Taylor: ', taylor)
	print('Percentage error taylor', 100 * (taylor - numpy.mean(wVw_sqrt)) / numpy.mean(wVw_sqrt), '%',
		  'Percentage error upper bound',
		  100 * ((numpy.trace(V) / n) ** 0.5 - numpy.mean(wVw_sqrt)) / numpy.mean(wVw_sqrt), '%')
	
	print()
	print('Numerical support,', min(wVw), 'to', max(wVw), 'Valid support for Taylor series', 2 * numpy.trace(V) / n)
	
	w, v = numpy.linalg.eig(V)
	
	#draw_sum_of_chi_cdf(w)
	exit()