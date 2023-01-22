import numpy, math

from drs import drs



def monte_carlo_simulations_lin_comb_chi(eigenvalues, iterations=1000):
	
	n = len(eigenvalues)
	
	I = numpy.zeros((n, n))
	numpy.fill_diagonal(I, 1)
	#wVw = []
	chi_list = numpy.empty([iterations])
	
	#three_term_taylor_list = []
	
	for count in range(iterations):
		z = numpy.random.default_rng().multivariate_normal(numpy.zeros(n), (1 / n) * I)
		
		chi_list[count] = (numpy.dot(eigenvalues, numpy.square(z)))**0.5
		#wVw.append(numpy.dot(w, numpy.dot(V, w)))
		#wVw_sqrt.append(wVw[-1] ** 0.5)
		
		#mu_q = numpy.trace(V) / n
		#three_term_taylor_list.append(((3 / 4) / mu_q ** 0.5) * wVw[-1] - ((1 / 8) / mu_q ** 1.5) * wVw[-1] * wVw[-1])
	
		
	
	
	#print('##############', numpy.var(three_term_taylor_list))
	
	#mean_wVw_monte_carlo = numpy.mean(wVw)
	mean_lin_comb_chi_monte_carlo = numpy.mean(chi_list)
	
	#var_wVw_monte_carlo = numpy.var(wVw)
	var_lin_comb_chi_monte_carlo = numpy.var(chi_list)
	
	#taylor_support_min_monte_carlo = min(wVw)
	#taylor_support_max_monte_carlo = max(wVw)
	return mean_lin_comb_chi_monte_carlo, var_lin_comb_chi_monte_carlo
	
	#mean_wVw_monte_carlo, var_wVw_monte_carlo, mean_sqrt_wVw_monte_carlo, var_sqrt_wVw_monte_carlo, \
	#	   taylor_support_min_monte_carlo, taylor_support_max_monte_carlo

def calculate_cumulant(k, eigenvalue_list):
	n = len(eigenvalue_list)
	kappa = (2 ** (k - 1)) * math.factorial(k - 1) / n ** k
	
	sum_product = 0
	for eigenvalue in eigenvalue_list:
		sum_product = sum_product + eigenvalue ** k
	
	return kappa * sum_product


def calculate_taylor2_mean_sqrt(eigenvalues):
	
	n = len(eigenvalues)
	eigenvalues_sum=sum(eigenvalues)
	
	u_Q_sqrt = calculate_cumulant(1, eigenvalues)**0.5
	
	mean_sqrt_taylor = u_Q_sqrt  - (1/8) * (1 / u_Q_sqrt)**3 * calculate_cumulant(2, eigenvalues)
	
	return mean_sqrt_taylor


def calculate_taylor_2_var_sqrt(eigenvalues):
	
	n = len(eigenvalues)
	
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
	
	
	
	mean_taylor1_errors = []
	mean_taylor2_errors =[]
	var_taylor2_errors = []
	var_taylor3_errors = []
	for counter in range(100):
		n = 100
		eigenvalues = drs(n,1)
		#print(eigenvalues)
	
		#n=2
		#eigenvalues=[0.04, 0.64, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04,]
		#eigenvalues=[0.4, 0.4, 0.2/8, 0.2/8, 0.2/8, 0.2/8, 0.2/8, 0.2/8, 0.2/8, 0.2/8,]
		mean_lin_comb_chi_monte_carlo, var_lin_comb_chi_monte_carlo = monte_carlo_simulations_lin_comb_chi(eigenvalues)
		print('Monte carlo mean', mean_lin_comb_chi_monte_carlo)
		
		mean_sqrt_taylor1 = calculate_cumulant(1, eigenvalues) ** 0.5
		print(mean_sqrt_taylor1)
		mean_taylor1_errors.append(100 * (mean_sqrt_taylor1 - mean_lin_comb_chi_monte_carlo) / mean_lin_comb_chi_monte_carlo)
		print('Mean error: ', mean_taylor1_errors[-1], '%')
		
		
		mean_sqrt_taylor2 = calculate_taylor2_mean_sqrt(eigenvalues)
		
		print(mean_sqrt_taylor2)
		mean_taylor2_errors.append(100* (mean_sqrt_taylor2-mean_lin_comb_chi_monte_carlo)/mean_lin_comb_chi_monte_carlo)
		print('Mean error: ',mean_taylor2_errors[-1], '%' )
		
		
		print('Monte carlo variance: ', var_lin_comb_chi_monte_carlo)
		
		var_sqrt_taylor_2 = calculate_taylor_2_var_sqrt(eigenvalues)
		print('Taylor 2 variance: ', var_sqrt_taylor_2)
		var_taylor2_errors.append(100 * (var_sqrt_taylor_2 - var_lin_comb_chi_monte_carlo) / var_lin_comb_chi_monte_carlo)
		print('Var Taylor 2  error: ',var_taylor2_errors[-1], '%')
		
		var_sqrt_taylor_3 = calculate_taylor_3_var_sqrt(eigenvalues)
		print('Taylor 3 variance: ', var_sqrt_taylor_3)
		var_taylor3_errors.append(100 * (var_sqrt_taylor_3 - var_lin_comb_chi_monte_carlo) / var_lin_comb_chi_monte_carlo)
		print('Var Taylor 3  error: ',
			  var_taylor3_errors[-1], '%')
	print()
	print('Final result')
	print('mean taylor 1 errors: ',numpy.mean(mean_taylor1_errors), numpy.std(mean_taylor1_errors), numpy.min(mean_taylor1_errors),
		  numpy.max(mean_taylor1_errors))
	print('mean taylor 2 errors: ', numpy.mean(mean_taylor2_errors), numpy.std(mean_taylor2_errors),
		  numpy.min(mean_taylor2_errors),
		  numpy.max(mean_taylor2_errors))
	print('var taylor 2 errors: ', numpy.mean(var_taylor2_errors), numpy.std(var_taylor2_errors),
		  numpy.min(var_taylor2_errors), numpy.max(var_taylor2_errors))
	print('var taylor 3 errors: ', numpy.mean(var_taylor3_errors), numpy.std(var_taylor3_errors),
		  numpy.min(var_taylor3_errors), numpy.max(var_taylor3_errors))