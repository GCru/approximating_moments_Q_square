import numpy
from sklearn import linear_model
from datetime import date
from statistics import stdev, mean
from math import pi, log, sqrt
from scipy.stats import kstest, kurtosis, skew, skewtest, jarque_bera, wasserstein_distance, median_abs_deviation
from scipy.special import erf
from scipy import optimize


def create_return_matrix():
	"We have three assets n=1,2,3 and four time steps t=1,2,3,4, , each R_j is one time step"
	
	R_1 = numpy.array([0.02,0.01, -0.01])
	R_2 = numpy.array([0.03, 0.02, 0.02])
	R_3 = numpy.array([0.04, 0.01, 0.01])
	R_4 = numpy.array([0.02, 0.03, 0.00])
	print(R_1)

	#Now put all in a numpy array
	
	R=numpy.array([R_1,R_2,R_3,R_4]).T
	
	print(R)
	print(R[0,0])
	
	return R

def create_factor_matrix():
	# at each time step t=1,2,3,4 there are two factors k=1,2
	
	f_1 = numpy.array([0.0067,0.01])
	f_2 = numpy.array([0.0233, 0.01])
	f_3 =numpy.array([0.02,0.02])
	f_4 = numpy.array([0.0167,0.02])
	
	f = numpy.array([f_1, f_2, f_3, f_4]).T
	
	print(f)
	
	return f
	
	
	
	

if __name__ == '__main__':
	#  see https://stackoverflow.com/questions/11479064/multiple-linear-regression-in-python
	
	R=create_return_matrix()
	
	f=create_factor_matrix()
	
	# Do multiple linear regression for each asset
	A = numpy.zeros(3)
	B = numpy.zeros((3,2))
	E = numpy.zeros((3,4))
	
	for i in range(0,3):
		# Count over assets

		F = f.T

		F = numpy.c_[numpy.ones(F.shape[0]),F]  # add a constant (the A in the model
		print(F)
	

		b,e = numpy.linalg.lstsq(F, R[i], rcond=None)[0:2]
		
		A[i]=b[0]
		B[i,0]=b[1]
		B[i,1]=b[2]

		print(b)
	
		print(e)
	
		#check answer
		r=[]
		e=[]
		for t in range(0,4):
			r.append(b[0]*F[t,0]+b[1]*F[t,1]+b[2]*F[t,2])
			e.append(r[-1]-R[i,t])
			E[i,t]= e[-1]
		
		print(r)
		print(sum(e))
		
	V_sample=numpy.cov(R, bias=False)
	
	V_factor = numpy.matmul(numpy.matmul(B,numpy.cov(f,bias=False)),B.T)+numpy.cov(E,bias=False)
	print(V_sample)
	
	print(V_factor)
	