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
	#print(w)
	# print(w)
	
	return V


def calculate_active_returns(V, z):
	# Caluclate minimum posibble bencmark return
	
	# Use own decompositioning
	own_LU = True
	
	if own_LU is True:
		# Own decomposition
		L = numpy.linalg.cholesky(V)
		return numpy.dot(L, z)
	
	else:
		# Using numpy function
		n = z.shape[0]
		answer = numpy.random.default_rng().multivariate_normal(numpy.zeros(n), V)
		
		return answer


def calculate_distribution_of_random_funds_returns():
	
	# for any set of returns r the distribution of fund returns where funds are drawn from
	# standard normal distribution then we get a normal distribution
	# bevause it is an affine transormation
	
	
	# draw randomm funds from standard multivariate normal
	return_list = []
	
	I = numpy.identity(4)
	r = numpy.array([2, 1, 0.01, 0.01])
	for idx in range(10000):
		answer = numpy.random.default_rng().multivariate_normal(numpy.zeros(4), I)
		
		return_list.append(numpy.dot(r, answer))
	
	print('Mean:' , numpy.mean(return_list), 'expected is 0')
	print('Variance',  numpy.var(return_list), 'expected is ', numpy.dot(r, r ))
	
	print('Kurtosis', 3+kurtosis(return_list), 'which should be 3')
	print('Skew', skew(return_list), 'which is zero')
	
	return


def calculate_variance_spread_from_V(V, sigma, rho):
	""" Calculate the variance of wVw with w a random normal with variance = 1/n. This is similar to calculate
		XX with X a draw form the N(0,V). Under certain conditons (e.g. rho=0) it follows a cgi squared distribtion
		with mean n/n and variance 2/n. See
		https://stats.stackexchange.com/questions/364446/distribution-of-the-sample-variance-of-values-from-a-multivariate-normal-distrib
	:param V:
	:return:
	"""
	
	n = V.shape[0]
	fund_var_list_w = []
	
	for i in range(1000):
	
		z = numpy.random.default_rng().normal(0.0, 1, size=V.shape[0])
		# set sigma^2_w = 1/n
		w = z / sqrt(n)  # same as normal(0.0,1/sqrt(n), size=n)
		
		fund_var_list_w.append(numpy.dot(w, numpy.dot(V, w)))
		
	
	print('Mean fund volatility', mean(fund_var_list_w), 'with expected value Trace(V)/n:', n*sigma * sigma/n)
	print('Var of fund volatility', stdev(fund_var_list_w)**2)
	print('Var of fund should be:', (2*sigma**4/n) + 2*(n-1)*(rho**2)/n)
	
	return


def see_how_random_fund_weights_work(n):
	"""
	
	:param n: Number of active positons
	:return:
	"""
	
	n = n
	
	stdev = 0.1
	I = numpy.zeros((n, n))
	numpy.fill_diagonal(I, 1)
	
	sum_list = []
	for count in range(100):
		w = numpy.random.default_rng().multivariate_normal(numpy.zeros(n), (stdev ** 2) / n * I)
		
		A_euclid = numpy.dot(w, w) ** 0.5
		adjust_w = w + (0 - numpy.sum(w)) / n
		A_euclid_adjust = numpy.dot(adjust_w, adjust_w) ** 0.5
		
		print(' The holdings add to', numpy.sum(w), 'with Euclidean active share', A_euclid,
			  'The adjusted holdings should add to zero', numpy.sum(adjust_w), ' with Euclidean active share',
			  A_euclid_adjust)
		
		print('EUclidean share difference', A_euclid - A_euclid_adjust)
		sum_list.append(numpy.sum(w))

	
	return
	

def examine_standardized_fund_return(V):
	""" Calculate the variance of wVw with w a random normal with variance = 1/n. This is similar to calculate
		XX with X a draw form the N(0,V). Under certain conditons (e.g. rho=0) it follows a cgi squared distribtion
		with mean n/n and variance 2/n. See
		https://stats.stackexchange.com/questions/364446/distribution-of-the-sample-variance-of-values-from-a-multivariate-normal-distrib
	:param V:
	:return:
	"""
	
	n= V.shape[0]
	fund_return_list_w = []
	fund_var_list_w = []
	standard_fund_return_w = []
	
	fund_var_5_perc_list =[]
	fund_var_95_perc_list = []
	
	standard_fund_mean_list =[]
	standard_fund_var_list =[]
	standard_fund_skew_list = []
	standard_fund_kurtosis_list = []
	
	# first calculate confidence intervals by doing many return draws
	draw_total = 1000
	for idx in range(draw_total):
		#draw the market return
		market_return = numpy.random.default_rng().multivariate_normal(numpy.zeros(n), 1.0 * V)
		# market_return = numpy.random.default_rng().multivariate_normal(numpy.zeros(n), numpy.identity(n))
		
		# create a distribution of fund returns
		for i in range(100):
			#market_return = numpy.random.default_rng().multivariate_normal(numpy.zeros(n), 1.0 * V)
			z = numpy.random.default_rng().normal(0.0, 1, size=V.shape[0])
			# set sigma^2_w = 1/n
			w = z / sqrt(n)  # same as normal(0.0,1/sqrt(n), size=n)
			
			fund_return_list_w.append(numpy.dot(w, market_return))
			fund_var_list_w.append(numpy.dot(w, numpy.dot(V, w)))
			
			standard_fund_return_w.append(fund_return_list_w[-1] / sqrt(fund_var_list_w[-1]))
		
		standard_fund_mean_list.append(numpy.mean(standard_fund_return_w))
		standard_fund_var_list.append(numpy.var(standard_fund_return_w))
		standard_fund_skew_list.append(skew(standard_fund_return_w))
		standard_fund_kurtosis_list.append(3 + kurtosis(standard_fund_return_w))
		
		print('Draw number', idx)
	print('varrr', numpy.var(standard_fund_return_w))
	fund_var_list_w.sort()
	# This should be the same independent of market returns
	fund_var_5_perc = fund_var_list_w[int(0.05*len(fund_var_list_w))]
	fund_var_95_perc = fund_var_list_w[int(0.95 * len(fund_var_list_w))]
	
	# outside bounds on standard fund var list would be
	# print(mean(fund_var_list_w), fund_var_list_w[50],fund_var_list_w[950])
	epsilon = ((mean(fund_var_list_w) - fund_var_5_perc ** 0.5) + (
			fund_var_95_perc ** 0.5 - mean(fund_var_list_w))) / 2
	print('Variance of fund variance ', numpy.var(fund_var_list_w))
	
	print('Epsilon', epsilon)
	
	standard_fund_var_list.sort()
	
	print('Mean of standard var list should be 1:', mean(standard_fund_var_list))
	# This would difffer according to market returns
	standard_fund_var_5_perc = standard_fund_var_list[int(0.05 * len(standard_fund_var_list))]
	standard_fund_var_95_perc = standard_fund_var_list[int(0.95 * len(standard_fund_var_list))]
	
	# the variance of fund variance should be related to epsilon via chi square distribution if \rho =0
	
	
	
		
	error_around_standard_variance = [(1 - epsilon) ** 2, (1 + epsilon) ** 2]
	
	print('Error bounds for standard fund return variance', error_around_standard_variance)
	
	print('Variance of standard fund returns', standard_fund_var_5_perc, standard_fund_var_95_perc)
	print()
	print('Standard fund mean over many return draws', mean(standard_fund_mean_list), 'vs expected 0')
	print('Standard fund variance over many return draws', mean(standard_fund_var_list), ' vs expected 1')
	print('Standard fund skew over many return draws', mean(standard_fund_skew_list), 'vs expected 0')
	print('Standard fund kurtosis over many return draws',mean(standard_fund_kurtosis_list), 'vs expected 3')
	
	#standard_fund_variance_list
	return

	
if __name__ == '__main__':
	
	# examine how random funds weights work
	#see_how_random_fund_weights_work(100)
	#exit()
	
	
	# check if funds drawn from standard random normal gives normal return distributions, e.g. affine transformation
	#calculate_distribution_of_random_funds_returns()
	#exit()
	
	# Set up covariance matrix of size n
	n = 100
	rho = 0.3
	sigma = 1.0
	V = create_covariance_matrix(sigma, rho, n)
	#V[1,1] =0.4
	
	
	
	
	I = numpy.zeros((n,n))
	numpy.fill_diagonal(I,1)
	
	wVw = []
	wVw_sqr = []
	adjusted_wVw= []
	fund_sigma_if_zero_sum = []
	for count in range(10000):
		
		w = numpy.random.default_rng().multivariate_normal(numpy.zeros(n), (1/n) * I)
		
		wVw.append(numpy.dot(w, numpy.dot(V, w)))
		wVw_sqr.append(wVw[-1]**0.5)
		
		adjusted_w = w +(0- numpy.sum(w)/n)
		adjusted_wVw.append(numpy.dot(adjusted_w, numpy.dot(V, adjusted_w)))
		
		#print(numpy.sum(w), wVw[-1], numpy.sum(adjusted_w), adjusted_wVw[-1])
	
		
		adjusted_A_euclid = numpy.dot(adjusted_w, adjusted_w) ** 0.5
		fund_sigma_if_zero_sum.append((adjusted_A_euclid * sigma * (1-rho)**0.5)**2)
	
	print('Mean of tracking error variance. Actual ', numpy.mean(wVw), 'Theoretical', numpy.trace(V)/n)
	
	print('Mean when weights are adjusted', numpy.mean(adjusted_wVw), 'Theoretical mean',  numpy.mean(fund_sigma_if_zero_sum) )
	
	print('Variance of tracking error variance. Actual:', numpy.var(wVw), 'Theoretical', 2*((rho**2)*(n-1)*n + (sigma**2)*n)/(n*n), 2*numpy.average(numpy.square(V)))
	
	print('Support,', min(wVw), max(wVw), 2*numpy.trace(V)/n )
	
	print('Mean of tracking error. Actual: ', numpy.mean(wVw_sqr), 'Theoretical:', (numpy.trace(V) / n)**0.5)
	print('Variance of tracking error. Actual:', numpy.var(wVw_sqr))\
	
	print('Check by calculating tracking error variance from formula', numpy.mean(wVw_sqr)**2+ numpy.var(wVw_sqr) )
	
	taylor=(numpy.trace(V)/n)**0.5 - (1.0/8.0)*(numpy.trace(V)/n)**(-1.5)*2*numpy.average(numpy.square(V))
	print('Taylor',taylor)
	print('Percentage error', 100*(taylor-numpy.mean(wVw_sqr))/numpy.mean(wVw_sqr),'%')
	
	# should give value close to 0.95, actually 0.94908
	w, v = numpy.linalg.eig(V)
	#print(sum(w)/n,2*numpy.trace(V)/n)
	#print(w)
	w=[i.real for i in w]
	x = hbe(coeff=w, x=2*sum(w))
	print('cut off point of distribution',x)
	
	exit()
	#draw histogram of vWv
	from make_histogram_procedure import make_histogram
	
	frequency_list, bin_edges_list = make_histogram(wVw_sqr)
	
	left_edges = bin_edges_list[:-1]
	right_edges = bin_edges_list[1:]
	data = {'frequency': frequency_list, 'left_edges': left_edges, 'right_edges': right_edges}
	source = ColumnDataSource(data=data)
	
	# Set up plot
	plot0 = figure(plot_width=int(500), plot_height=500)
	r0 = plot0.quad(
		bottom=0, top='frequency', left='left_edges', right='right_edges', source=source,
		fill_color='lightgrey', line_color='black')  # legend='Underperforming monkey portfolios')
	
	show(plot0)
	
	exit()
	
	print((numpy.std(wVw))**2, 2*numpy.linalg.norm(V)**2/(n*n), numpy.std(fund_sigma_if_zero_sum)**2)
	
	print(numpy.mean(adjust_wVw), numpy.std(adjust_wVw)**2)
	
	#print(numpy.mean(w), numpy.std(w), (1/n)**0.5)
	
	
	# what is expected if sum is zero
	print('sum is zero', )
	exit()
	
	
	
	
	# check spread of variances calculated from covariance matrix
	#calculate_variance_spread_from_V(V, sigma, rho)
	#exit()
	examine_standardized_fund_return(V)
	exit()
	
	
	# Note that by making matrix larger does not improve accuracy as the size of V grow, so even
	# if the sum of z's approaches zero as n increases
	
	
	# create covariance matrix V
	
	market_return = numpy.random.default_rng().multivariate_normal(numpy.zeros(n), 1 * V)
	# print('dotproduct', numpy.dot(market_return, market_return))
	
	cross_sect_with_z_list = []
	cross_sect_with_w_list = []
	fund_var_list_z = []
	fund_var_list_w = []
	standardized_fund_return_list = []
	
	expected_value_with_z_list = []
	expected_value_with_w_list = []
	cross_z_list = []
	for i in range(10):
		# set sigma^2_w = 1/n
		
		z = numpy.random.default_rng().normal(0.0, 1, size=n)
		
		print()
		# Normalise
		# z = z - numpy.mean(z)
		print("Norm:  z*z", numpy.dot(z, z), 'should equal', n)
		
		# create cross z list
		# for j in range(n):
		#	for k in range(j+1,n):
		#		cross_z_list.append(z[j]*z[k])
		
		w = z / sqrt(n)  # same as normal(0.0,1/sqrt(n), size=n)
		print("Norm w*w", numpy.dot(w, w), 'should equal', 1)
		
		fund_var_list_z.append(numpy.dot(z, numpy.dot(V, z)))
		print('Fund risk with z', fund_var_list_z[-1])
		
		fund_var_list_w.append(numpy.dot(w, numpy.dot(V, w)))
		print('Fund risk with w', fund_var_list_w[-1])
		
		# market_return = numpy.random.default_rng().normal(0,sigma,n)
		# market_return = calculate_active_returns(V, w)
		
		standardized_fund_return_list.append(numpy.dot(w, market_return) / sqrt(fund_var_list_w[-1]))
		
		second_term_with_z = sigma * sigma * n * numpy.dot(z, z) / n
		print("Trace term with z", second_term_with_z)
		
		second_term_with_w = sigma * sigma * n * numpy.dot(w, w) / n
		print("Trace term with w", second_term_with_w)
		
		expected_value_with_z_list.append(abs(fund_var_list_z[-1] - second_term_with_z))
		print("Expected value with z", expected_value_with_z_list[-1])
		
		expected_value_with_w_list.append(abs(fund_var_list_w[-1] - second_term_with_w))
		print("Expected value with w", expected_value_with_w_list[-1])
		
		r_active_returns_with_z = calculate_active_returns(V, z)
		# tecnically the active retuns should add to zero
		r_active_returns_with_z = r_active_returns_with_z - numpy.mean(r_active_returns_with_z)
		
		cross_section_with_z = numpy.dot(r_active_returns_with_z, r_active_returns_with_z)
		cross_sect_with_z_list.append(cross_section_with_z)
		
		r_active_returns_with_w = calculate_active_returns(V, w)
		# tecnically the active retuns should add to zero
		r_active_returns_with_w = r_active_returns_with_w - numpy.mean(r_active_returns_with_w)
		
		cross_section_with_w = numpy.dot(r_active_returns_with_w, r_active_returns_with_w)
		cross_sect_with_w_list.append(cross_section_with_w)
	
	print()
	# print('St dev', stdev(cross_sect_list), mean(cross_sect_list))
	
	rhs_bound_z = log(n) * log(n) * sqrt(n)
	print('Mean expected value with z list', mean(expected_value_with_z_list), 'K value ',
		  mean(expected_value_with_z_list) / rhs_bound_z)
	
	print('Mean fund volatility', mean(fund_var_list_z), 'var with z', stdev(fund_var_list_z))
	print('Trace with z', n * sigma * sigma)
	print('Fund volatility equals cross sect with z ', sigma * sigma * (1 - rho) * n)
	print('Cross section calculated with z ', mean(cross_sect_with_z_list))
	
	print()
	rhs_bound_w = log(n) * log(n) / sqrt(n)
	error_vs_trace = [abs(fund_var_list_w[idx] - sigma * sigma) for idx in range(len(fund_var_list_w))]
	
	print('Mean expected value with w list', mean(expected_value_with_w_list), 'K value ',
		  mean(expected_value_with_w_list) / rhs_bound_w)
	print('Mean fund volatility', mean(fund_var_list_w),
		  'var with w', stdev(fund_var_list_w))
	print('Trace with w', sigma * sigma)
	print('Fund volatility equals cross sect with w ', sigma * sigma * (1 - rho))
	print('Cross section calculated with w ', mean(cross_sect_with_w_list))
	
	print('Mean_error', mean(error_vs_trace))
	
	print("norm matrix", sigma * sigma + n * rho * sigma * sigma)
	# print('Mean cross z list', mean(cross_z_list), len(cross_z_list), n*(n-1))
	
	print('Stdev van standarzidez fund', stdev(standardized_fund_return_list))
	
	# draw histogram of fund_stdev_list
	from shared_with_active_risk_paper import make_histogram
	
	real_stdev_list = [sqrt(fund_var_list_w[idx]) for idx in range(len(fund_var_list_w))]
	frequency_list, bin_edges_list = make_histogram(real_stdev_list)
	print(stdev(real_stdev_list), median_abs_deviation(real_stdev_list))
	
	left_edges = bin_edges_list[:-1]
	right_edges = bin_edges_list[1:]
	data = {'frequency': frequency_list, 'left_edges': left_edges, 'right_edges': right_edges}
	source = ColumnDataSource(data=data)
	
	# Set up plot
	plot0 = figure(plot_width=int(500), plot_height=500)
	r0 = plot0.quad(
		bottom=0, top='frequency', left='left_edges', right='right_edges', source=source,
		fill_color='lightgrey', line_color='black')  # legend='Underperforming monkey portfolios')
	
	show(plot0)
	
	# draw histogram of inverse
	from shared_with_active_risk_paper import make_histogram
	
	frequency_list, bin_edges_list = make_histogram(standardized_fund_return_list)
	
	left_edges = bin_edges_list[:-1]
	right_edges = bin_edges_list[1:]
	data = {'frequency': frequency_list, 'left_edges': left_edges, 'right_edges': right_edges}
	source = ColumnDataSource(data=data)
	
	# Set up plot
	plot1 = figure(plot_width=int(500), plot_height=500)
	r0 = plot1.quad(
		bottom=0, top='frequency', left='left_edges', right='right_edges', source=source,
		fill_color='lightgrey', line_color='black')  # legend='Underperforming monkey portfolios')
	my_plot = row(plot0, plot1)
# show(plot1)

