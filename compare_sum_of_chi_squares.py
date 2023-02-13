import numpy
from datetime import date
from statistics import stdev, mean
from math import pi, log, sqrt
from scipy.stats import chi2, gamma, kstest, kurtosis, skew, skewtest, jarque_bera, wasserstein_distance, \
	median_abs_deviation
from scipy.special import erf
from scipy import optimize

from momentchi2 import hbe

if __name__ == '__main__':
	

	
	c = 1
	n = 10
	
	print('chi_sqaure', chi2.cdf(x=10, df=n, scale=c))
	
	print('gamma', gamma.cdf(x=10, a=n / 2, scale=2 * c))
	
	w_max = numpy.array([0.0, 0.0, 1.0])
	w_min = numpy.array([1/3, 1/3, 1/3])
	
	w_inbetween = numpy.array([0.60, 0.2, 0.2 ])

	size = 10000
	result_min = numpy.empty(size)
	result_max = numpy.empty(size)
	result_inbetween = numpy.empty(size)
	
	for idx in range(size):
		chi_square_list = numpy.random.chisquare(df=1,size=3)
		
		result_min[idx] = numpy.dot(w_min,chi_square_list)
		result_max[idx] = numpy.dot(w_max, chi_square_list)
		result_inbetween[idx] = numpy.dot(w_inbetween, chi_square_list)
	
		
	
	result_min =result_min[result_min<0.05]
	result_max = result_max[result_max < 0.05]
	result_inbetween = result_inbetween[result_inbetween < 0.05]
	
	print('Density', 100* len(result_min)/size, 100*len(result_inbetween)/size, 100*len(result_max)/size)
	
	
	
	
