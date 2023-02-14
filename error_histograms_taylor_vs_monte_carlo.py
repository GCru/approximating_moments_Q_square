import numpy, math

from drs import drs


def monte_carlo_simulations_lin_comb_chi(eigenvalues, iterations=10000):
	
	n = len(eigenvalues)
	
	I = numpy.zeros((n, n))
	numpy.fill_diagonal(I, 1)

	chi_list = numpy.empty([iterations])
	
	for count in range(iterations):
		z = numpy.random.default_rng().multivariate_normal(numpy.zeros(n), (1 / n) * I)
		
		chi_list[count] = (numpy.dot(eigenvalues, numpy.square(z)))**0.5
		
	mean_lin_comb_chi_monte_carlo = numpy.mean(chi_list)
	var_lin_comb_chi_monte_carlo = numpy.var(chi_list)
	
	return mean_lin_comb_chi_monte_carlo, var_lin_comb_chi_monte_carlo
	

def calculate_cumulant(k, eigenvalue_list):
	
	n = len(eigenvalue_list)
	kappa = (2 ** (k - 1)) * math.factorial(k - 1) / n ** k
	
	sum_product = 0
	for eigenvalue in eigenvalue_list:
		sum_product = sum_product + eigenvalue ** k
	
	return kappa * sum_product


def calculate_taylor_3_mean_sqrt(eigenvalues):
	
	u_Q_sqrt = calculate_cumulant(1, eigenvalues)**0.5
	
	mean_sqrt_taylor = u_Q_sqrt  - (1/8) * (1 / u_Q_sqrt)**3 * calculate_cumulant(2, eigenvalues)
	
	return mean_sqrt_taylor


def calculate_taylor_2_var_sqrt(eigenvalues):
	
	var_sqrt_wVw_taylor_2 = calculate_cumulant(2, eigenvalues) / (4*calculate_cumulant(1, eigenvalues))
	
	return var_sqrt_wVw_taylor_2


def calculate_taylor_3_var_sqrt(eigenvalues):
	
	kappa_list = []
	
	for idx in range(5):
		if idx == 0:
			kappa_list.append(0)
		else:
			kappa_list.append(calculate_cumulant(idx, eigenvalues))
	
	var_Q = kappa_list[2]
	
	var_Q_sqr = kappa_list[4] + 4 * kappa_list[3] * kappa_list[1] + 2 * kappa_list[2] ** 2 + 4 * kappa_list[2] * \
				kappa_list[1] ** 2
	
	cov_Q_Q_sqr = kappa_list[3] + 2 * kappa_list[2] * kappa_list[1]
	
	term1 = (9 / 16) * var_Q / kappa_list[1]
	
	term2 = (1 / 64) * var_Q_sqr / kappa_list[1] ** 3
	
	term3 = (3 / 16) * cov_Q_Q_sqr / kappa_list[1] ** 2
	
	return term1 + term2 - term3


if __name__ == '__main__':
	
	
	
	mean_taylor_2_errors = []
	mean_taylor_3_errors =[]
	var_taylor_2_errors = []
	var_taylor_3_errors = []
	
	sum_monte_carlo_variance =0
	for counter in range(50):
		n = 3
		eigenvalues = drs(n,1)
		#eigenvalues=[0.1,0.9]
		print()
		print(eigenvalues)
		#n=3
		#eigenvalues=eigenvalues+[0]
		#print(eigenvalues)
	
		#n=2
		#eigenvalues=[0.5,0.5]
		#eigenvalues=[1,0,0,0,0,0,0,0,0,0]
		#eigenvalues=[0.4, 0.4/9, 0.4/9, 0.4/9, 0.4/9,0.4/9, 0.4/9,0.4/9, 0.4/9,0.4/9]
		#eigenvalues=[0.44, 0.44, 0.02/8, 0.02/8, 0.02/8, 0.02/8, 0.02/8, 0.02/8, 0.02/8, 0.02/8,]
		#eigenvalues=[0.1, 0.1, 0.1, 0.1,0.1,0.1,0.1,0.1,0.1,0.1]
		#eigenvalues = [0.1]*10
		mean_lin_comb_chi_monte_carlo, var_lin_comb_chi_monte_carlo = monte_carlo_simulations_lin_comb_chi(eigenvalues)
		print('Monte carlo mean', mean_lin_comb_chi_monte_carlo, 'expec6ted mean', 0.1**0.5*(1-1/40))
		print('Monte carlo var', var_lin_comb_chi_monte_carlo, 'expec6ted var', 0.1  * (1/20))
		
		mean_sqrt_taylor_2 = calculate_cumulant(1, eigenvalues) ** 0.5
		mean_taylor_2_errors.append(100 * (mean_sqrt_taylor_2 - mean_lin_comb_chi_monte_carlo) / mean_lin_comb_chi_monte_carlo)
		print('Mean taylor 2 error: ', mean_taylor_2_errors[-1], '%')
		
		mean_sqrt_taylor_3 = calculate_taylor_3_mean_sqrt(eigenvalues)
		mean_taylor_3_errors.append(100* (mean_sqrt_taylor_3-mean_lin_comb_chi_monte_carlo)/mean_lin_comb_chi_monte_carlo)
		print('Mean taylor 3 error: ',mean_taylor_3_errors[-1], '%' )
		
		
		sum_monte_carlo_variance =sum_monte_carlo_variance+ var_lin_comb_chi_monte_carlo
		
		var_sqrt_taylor_2 = calculate_taylor_2_var_sqrt(eigenvalues)
		#print('Var Taylor 2: ', var_sqrt_taylor_2)
		var_taylor_2_errors.append(100 * (var_sqrt_taylor_2 - var_lin_comb_chi_monte_carlo) / var_lin_comb_chi_monte_carlo)
		print('Var Taylor 2 error: ',var_taylor_2_errors[-1], '%')
		
		var_sqrt_taylor_3 = calculate_taylor_3_var_sqrt(eigenvalues)
		#print('Taylor 3 variance: ', var_sqrt_taylor_3)
		var_taylor_3_errors.append(100 * (var_sqrt_taylor_3 - var_lin_comb_chi_monte_carlo) / var_lin_comb_chi_monte_carlo)
		print('Var Taylor 3  error: ',  var_taylor_3_errors[-1], '%')
	print()
	print(sum_monte_carlo_variance/100)
	print('Final result')
	print('Mean taylor 2 errors: ',numpy.mean(mean_taylor_2_errors), numpy.std(mean_taylor_2_errors), numpy.min(mean_taylor_2_errors),
		  numpy.max(mean_taylor_2_errors))
	print('Mean taylor 3 errors: ', numpy.mean(mean_taylor_3_errors), numpy.std(mean_taylor_3_errors),
		  numpy.min(mean_taylor_3_errors),
		  numpy.max(mean_taylor_3_errors))
	print('Var taylor 2 errors: ', numpy.mean(var_taylor_2_errors), numpy.std(var_taylor_2_errors),
		  numpy.min(var_taylor_2_errors), numpy.max(var_taylor_2_errors))
	print('Var taylor 3 errors: ', numpy.mean(var_taylor_3_errors), numpy.std(var_taylor_3_errors),
		  numpy.min(var_taylor_3_errors), numpy.max(var_taylor_3_errors))
	
