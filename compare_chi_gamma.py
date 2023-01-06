import numpy
from datetime import date
from statistics import stdev, mean
from math import pi, log, sqrt
from scipy.stats import chi2, gamma, kstest, kurtosis, skew, skewtest, jarque_bera, wasserstein_distance, median_abs_deviation
from scipy.special import erf
from scipy import optimize

from momentchi2 import hbe

if __name__ == '__main__':

	c=1
	n=10
	
	print('chi_sqaure', chi2.cdf(x=10,df=n, scale=c))
	
	print('gamma', gamma.cdf(x=10, a=n/2,scale=2*c))