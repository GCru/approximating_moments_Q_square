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
	
	lambda_list =[0,0,0]
	for i1 in range(0,10):
		
		for i2 in range(0,10):
			
			for i3 in range(0,10):
				
				lambda_list[0] = (i1+1) * 0.0001
				lambda_list[1] = (i2) * 0.0001
				lambda_list[2] = (i3) * 0.0001
				
				lambda_squares =0
				lambda_cubes =0
				for j in range(0, len(lambda_list)):
					lambda_squares = lambda_squares + lambda_list[j]**2
					lambda_cubes = lambda_cubes + lambda_list[j] ** 3
					
				result = lambda_squares**3/ lambda_cubes**2
				
				print(lambda_list, result)
				if result<1:
					input()
					
				