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

def test_Marcenko_Pastur_covariance_matrix(sigma,n,T):
	"""
	:param sigma: standard deviation
	:param n: number of assets
	:param T: number of time intervals
	:return:
	"""
	
	R = numpy.random.default_rng().normal(loc=0, scale = sigma, size=(n,T))
	V = numpy.matmul(R, numpy.transpose(R))/T
	
	w, v = numpy.linalg.eig(V)
	print('Max eigenvalue:', max(w), ' Min eigenvalue', min(w))
	print('Theoretical max eigenvalue', (1 + (n / T) ** 0.5) ** 2)
	
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


def draw_sum_of_chi_cdf(w):
	
	
	
	grain=100
	
	
	plot0 = figure(plot_width=int(500), plot_height=500)
	
	n = 100
	rho = 0.94
	sigma = 1.0
	V = create_covariance_matrix(sigma, rho, n)
	w, v = numpy.linalg.eig(V)
	w = [i.real / n for i in w]
	
	w=[0.005,0.005,0.005, 0.995]
	w=[0.01,0.01,0.01,0.97]
	#w=[0.2,0.2,0.2,0.2,0.2]
	w=[0.40,0.2,0.3,0.1]
	print('*****', hbe(coeff=w,x=1))
	
	#w=[0.499,0.001,0.499,0.001]
	print(sum(w))
	x_axis = [i*4* sum(w) / grain for i in range(1, grain + 1)]
	
	y_axis = [hbe(coeff=w, x=item) for item in x_axis]
	plot0.line(x_axis,y_axis, line_width=2)
	
	rho = 0.75
	V = create_covariance_matrix(sigma, rho, n)
	w, v = numpy.linalg.eig(V)
	w = [i.real / n for i in w]
	w = [0.106, 0.106, 0.106, 0.682]
	w=[0.25,0.25,0.25,0.25]
	#w=[0.02,0.02,0.01,0.95]
	print(sum(w))
	x_axis = [i *4* sum(w) / grain for i in range(1, grain + 1)]
	
	y_axis = [hbe(coeff=w, x=item) for item in x_axis]
	plot0.line(x_axis, y_axis, line_width=2, line_color="red" )
	
	show(plot0)
	
	return





	
if __name__ == '__main__':
	
	# examine how random funds weights work
	#see_how_random_fund_weights_work(100)
	#exit()
	
	
	# check if funds drawn from standard random normal gives normal return distributions, e.g. affine transformation
	#calculate_distribution_of_random_funds_returns()
	#exit()
	
	# check how marcenko pastur works
	#sigma =1
	#n=1000
	#T=1000
	#V= test_Marcenko_Pastur_covariance_matrix(sigma, n, T)
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
	wVw_sqrt = []
	adjusted_wVw= []
	fund_sigma_if_zero_sum = []
	for count in range(1000):
		
		w = numpy.random.default_rng().multivariate_normal(numpy.zeros(n), (1/n) * I)
		
		wVw.append(numpy.dot(w, numpy.dot(V, w)))
		wVw_sqrt.append(wVw[-1]**0.5)
		
		adjusted_w = w +(0- numpy.sum(w)/n)
		adjusted_wVw.append(numpy.dot(adjusted_w, numpy.dot(V, adjusted_w)))
		
		#print(numpy.sum(w), wVw[-1], numpy.sum(adjusted_w), adjusted_wVw[-1])
	
		
		adjusted_A_euclid = numpy.dot(adjusted_w, adjusted_w) ** 0.5
		fund_sigma_if_zero_sum.append((adjusted_A_euclid * sigma * (1-rho)**0.5)**2)
	
	print()
	print('Mean of tracking error variance. Actual ', numpy.mean(wVw), 'Theoretical', numpy.trace(V)/n)
	
	#print('Mean when weights are adjusted', numpy.mean(adjusted_wVw), 'Theoretical mean',  numpy.mean(fund_sigma_if_zero_sum) )
	
	print('Variance of tracking error variance. Actual:', numpy.var(wVw), 'Theoretical', 2*((rho**2)*(n-1)*n + (sigma**2)*n)/(n*n), 'or', 2*numpy.average(numpy.square(V)))
	
	print()
	print('Expected value (mean) of tracking error. Actual: ', numpy.mean(wVw_sqrt), 'Theoretical upper bound:', (numpy.trace(V) / n)**0.5)

	# print('Check by calculating mean tracking error variance from formula', numpy.mean(wVw_sqrt)**2+ numpy.var(wVw_sqrt) )
	
	taylor=(numpy.trace(V)/n)**0.5 - (1.0/8.0)*(numpy.trace(V)/n)**(-1.5)*2*numpy.average(numpy.square(V))
	print('Expected value (mean) of tracking error. Taylor: ', taylor)
	print('Percentage error taylor', 100*(taylor-numpy.mean(wVw_sqrt))/numpy.mean(wVw_sqrt),'%',
		  'Percentage error upper bound',100*((numpy.trace(V)/n)**0.5-numpy.mean(wVw_sqrt))/numpy.mean(wVw_sqrt), '%')
	
	print()
	print('Numerical support,', min(wVw), 'to', max(wVw), 'Valid support for Taylor series', 2*numpy.trace(V)/n )
	
	w, v = numpy.linalg.eig(V)
	
	draw_sum_of_chi_cdf(w)
	exit()
	#print(sum(w)/n,2*numpy.trace(V)/n)
	#print(w)
	w=[i.real/n for i in w]
	x = hbe(coeff=w, x=2*sum(w))
	print('Probability that wVw >', 2*sum(w),' for wVw distribution', x)
	
	x = hbe(coeff=w, x=0.5*sum(w))
	print('Probability that wVw >', 0.001, ' for wVw distribution', x)
	
	print()
	print('Variance of tracking error: Estimated from Taylor', numpy.trace(V)/n-taylor**2, 'and stdev', (numpy.trace(V)/n-taylor**2)**0.5)
	print('Variance of tracking error: Actual:', numpy.var(wVw_sqrt), 'and stdev', numpy.var(wVw_sqrt)**0.5 )
	
	print('Check actual variance of tracking error from formula:',  numpy.mean(wVw) - numpy.mean(wVw_sqrt)**2)
	
	mu_x= numpy.trace(V)/n
	sigma_sq = 2*numpy.average(numpy.square(V))
	print((sigma_sq *0.25/mu_x)-(sigma_sq**2/(4*16*mu_x**3)))
 
	#exit()
	#draw histogram of vWv
	from make_histogram_procedure import make_histogram
	
	frequency_list, bin_edges_list = make_histogram(wVw_sqrt)
	
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
	
	
	


