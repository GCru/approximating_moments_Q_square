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


def calculate_taylor_mean_sqrt_wVw(eigenvalues):
	
	n = len(eigenvalues)
	eigenvalues_sum=sum(eigenvalues)
	
	mean_sqrt_wVw_taylor = (numpy.trace(V) / n) ** 0.5 - (1.0 / 8.0) * (numpy.trace(V) / n) ** (
		-1.5) * 2 * numpy.average(
		numpy.square(V))
	
	return mean_sqrt_wVw_taylor


def calculate_taylor_2_var_sqrt_wVw(V):
	n = numpy.shape(V)[0]
	
	var_sqrt_wVw_taylor_2 = (numpy.linalg.norm(V) ** 2) / (2 * n * numpy.trace(V))
	
	return var_sqrt_wVw_taylor_2


def calculate_taylor_3_var_sqrt_vWv(V):
	eigenvalue_list, _ = numpy.linalg.eig(V)
	
	kappa_list = []
	
	for idx in range(5):
		if idx == 0:
			kappa_list.append(0)
		else:
			kappa_list.append(calculate_cumulant(idx, eigenvalue_list))
	
	var_Q = kappa_list[2]
	
	var_Q_sqr = kappa_list[4] + 4 * kappa_list[3] * kappa_list[1] + 2 * kappa_list[2] ** 2 + 4 * kappa_list[2] * \
				kappa_list[1] ** 2
	
	cov_Q_Q_sqr = kappa_list[3] + 2 * kappa_list[2] * kappa_list[1]
	
	term1 = (9 / 16) * var_Q / kappa_list[1]
	
	term2 = (1 / 64) * var_Q_sqr / kappa_list[1] ** 3
	
	term3 = (3 / 16) * cov_Q_Q_sqr / kappa_list[1] ** 2
	
	return term1 + term2 - term3


if __name__ == '__main__':
	
	
	n = 10
	eigenvalues = drs(2,1)

	eigenvalues=[0.1]*10
	mean_lin_comb_chi_monte_carlo, var_lin_comb_chi_monte_carlo = monte_carlo_simulations_lin_comb_chi(eigenvalues)
	
	print(mean_lin_comb_chi_monte_carlo, var_lin_comb_chi_monte_carlo)
	
	