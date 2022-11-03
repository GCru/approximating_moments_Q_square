import numpy
from datetime import date
from statistics import stdev, mean
from math import pi, log, sqrt
from scipy.stats import kstest, kurtosis, skew, skewtest, jarque_bera, wasserstein_distance, median_abs_deviation
from scipy.special import erf
from scipy import optimize




def create_random_n_by_T_matrix(n,T):

	return numpy.random.default_rng().normal(scale=2, size=(n,T))


def create_cov_matrix(X,n ):

	return numpy.cov(X)

	return numpy.matmul(X,numpy.transpose(X) ) / T

if __name__ == '__main__':
	
	#v= numpy.random.normal(size=(10000))
	#w= numpy.random.normal(size=(10000))
	#v_square = numpy.dot(v,v)
	#w_times_v =numpy.dot(w,v)
	#print(v_square, w_times_v)
	#exit()
	
	n=200
	T=25
	
	X = create_random_n_by_T_matrix(n,T)
	
	Y= create_cov_matrix(X,n)
	print(Y.shape)
	
	
	
	l , _ = numpy.linalg.eig(Y)
	print('Maximum eigenvalue',max(l), 'vs Marchenko-Pastur', 4* (1+(n/T)**0.5)**2)
	
	D = numpy.zeros((n, n))
	numpy.fill_diagonal(D, 1)
	
	
	
	print('Overall sum', numpy.sum(Y))
	
	
	print('Diagonal sum', numpy.sum(numpy.diag(Y)))
	print(numpy.linalg.norm(Y-D, ord=2))
	
	
	print('Average of all entries sqaured', numpy.average(Y)-1/n )