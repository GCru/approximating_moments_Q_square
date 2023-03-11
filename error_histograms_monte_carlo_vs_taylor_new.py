import numpy, math
import csv

from drs import drs


def update_monte_carlo_file(n, fname, total_eigenvalue_draws, iterations_per_eigenvalue=100000):
	
	from pathlib import Path
	# fname='monte_carlo-'+str(n)+'.npy'
	fname = Path('test.npy')
	if fname.is_file():
		print('Exists')
		with open(fname, 'rb') as fp:
			num_lines=0
			while True:
				try:
					item = numpy.load(fp)
					print(item)
					num_lines= num_lines+1
				except:
					print("EoF", num_lines)
					break
			
			#num_lines = sum(1 for line in fp)
			#print('Total lines:', num_lines)  # 8
			#for line in fp:
			#	print(line)
		total_remaining_eigenvalue_draws = total_eigenvalue_draws - num_lines
	else:
		print('Not found')
		num_lines=0
		total_remaining_eigenvalue_draws = total_eigenvalue_draws
	
	for count_eigenvalue_draws in range(total_remaining_eigenvalue_draws):
		eigenvalues = drs(n, 1)

		I = numpy.zeros((n, n))
		numpy.fill_diagonal(I, 1)
	
		chi_list = numpy.empty([iterations])
	
		for count in range(iterations_per_eigenvalue):
			z = numpy.random.default_rng().multivariate_normal(numpy.zeros(n), (1 / n) * I)
			# print('stuck', count)
			chi_list[count] = (numpy.dot(eigenvalues, numpy.square(z))) ** 0.5
	
		mean_lin_comb_chi_monte_carlo = numpy.mean(chi_list)
		var_lin_comb_chi_monte_carlo = numpy.var(chi_list)
	
		print('Eigenvalue_set: ', num_lines + count_eigenvalue_draws+1 )
	
		output_array = numpy.concatenate((eigenvalues, numpy.array([mean_lin_comb_chi_monte_carlo, var_lin_comb_chi_monte_carlo])))
		with open('test.npy', 'ab') as f:
			numpy.save(f, output_array)
			print(output_array)
	
	return


def calculate_cumulant(k, eigenvalue_list):
	n = len(eigenvalue_list)
	kappa = (2 ** (k - 1)) * math.factorial(k - 1) / n ** k
	
	sum_product = 0
	for eigenvalue in eigenvalue_list:
		sum_product = sum_product + eigenvalue ** k
	
	return kappa * sum_product


def calculate_taylor_3_mean_sqrt(eigenvalues):
	u_Q_sqrt = calculate_cumulant(1, eigenvalues) ** 0.5
	
	mean_sqrt_taylor = u_Q_sqrt - (1 / 8) * (1 / u_Q_sqrt) ** 3 * calculate_cumulant(2, eigenvalues)
	
	return mean_sqrt_taylor


def calculate_taylor_2_var_sqrt(eigenvalues):
	var_sqrt_wVw_taylor_2 = calculate_cumulant(2, eigenvalues) / (4 * calculate_cumulant(1, eigenvalues))
	
	# print(var_sqrt_wVw_taylor_2)
	# input('aha')
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


def calculate_extreme_mean_and_var(eigenvalue_sum, n):
	""" Calculate analytical values for mean and var when all eigenvalues are euqal
	:param eigenvaluesum:
	:param n:
	:return:
	"""
	
	mu_Q = eigenvalue_sum / n
	mean_extreme = ((2 * mu_Q / n) ** 0.5) * math.gamma((n + 1) / 2) / math.gamma(n / 2)
	# print(mean_extreme, mean_extreme / mu_Q ** 0.5)
	# exit()
	var_extreme = (2 ** 0.5 * math.gamma((n + 1) / 2) / math.gamma(n / 2)) ** 2
	
	var_extreme = mu_Q * (1 - var_extreme / n)
	
	# print(var_extreme, var_extreme / mu_Q)
	
	return mean_extreme, var_extreme


def calculate_extremes(n):
	
	if n>102:
		print('The value of n too large to calculate')
		mean_two_term_extreme_error = 0
		mean_three_term_extreme_error= 0
		var_two_term_extreme_error = 0
		var_three_tern_extreme_error =0
		return
	else:
		mean_extreme, var_extreme = calculate_extreme_mean_and_var(eigenvalue_sum, n)
		mu_Q = eigenvalue_sum / n
		answer = mu_Q ** 0.5
		mean_two_term_extreme_error = 100 * (answer - mean_extreme) / mean_extreme
		
		answer = (1 - 1 / (4 * n)) * mu_Q ** 0.5
		mean_three_term_extreme_error = 100 * (answer - mean_extreme) / mean_extreme
		
		answer = (1 - 1 / (4 * n) + 1 / (2 * n * n)) * mu_Q ** 0.5
		mean_four_term_extreme_error = 100 * (answer - mean_extreme) / mean_extreme
	
		answer = mu_Q / (2 * n)
		var_two_term_extreme_error = 100 * (answer - var_extreme) / var_extreme
		
		answer = (mu_Q / n) * (1 / 2 - 7 / (8 * n) + 3 / (4 * n ** 2))
		var_three_term_extreme_error = 100 * (answer - var_extreme) / var_extreme
	
		print('Mean two term extreme error', mean_two_term_extreme_error)
		print('Mean three term extreme error', mean_three_term_extreme_error)
		print('Mean four term extreme error', mean_four_term_extreme_error)
		print('Var two term extreme error', var_two_term_extreme_error)
		print('Var three term extreme error', var_three_term_extreme_error)
	
	return mean_two_term_extreme_error, mean_three_term_extreme_error, var_two_term_extreme_error, var_three_term_extreme_error


if __name__ == '__main__':
	
	results_list = []
	
	sum_monte_carlo_variance = 0
	
	eigenvalue_sum = 1
	
	n=7 # number of eigenvalues
	
	mean_two_term_extreme_error, mean_three_term_extreme_error, var_two_term_extreme_error, \
		var_three_term_extreme_error = calculate_extremes(n)

	# create file with monte carlo expected value and variance
	
	# first see if file exists and check how many items
	#fname='monte_carlo-'+str(n)+'.npy'
	fname='test.npy'
	update_monte_carlo_file(n, fname, total_eigenvalue_draws=20,iterations_per_eigenvalue=10)
	
	
	with open(fname, 'rb') as fp:
		num_lines = 0
		while True:
			try:
				item = numpy.load(fp)
				mean_lin_comb_chi_monte_carlo = item[-2]
				var_lin_comb_chi_monte_carlo = item[-1]
				#print(mean_lin_comb_chi_monte_carlo, var_lin_comb_chi_monte_carlo)
				eigenvalues = item[list(range(0,n))]
				#print(eigenvalues, sum(eigenvalues))
				#print(item)
				num_lines = num_lines + 1
			except:
				print("EoF", num_lines)
				break
	
	exit()
	

	
	max_mean_taylor_3_errors = 0
	max_var_taylor_3_errors = 0
	
	mean_taylor_2_errors = []
	mean_taylor_3_errors = []
	
	var_taylor_2_errors = []
	var_taylor_3_errors = []
	
	for counter in range(10000):
		
		print(counter)
		
		eigenvalues = drs(n, eigenvalue_sum)
		# eigenvalues=[0.1,0.9]
		print()
		print(eigenvalues)
		
		mean_lin_comb_chi_monte_carlo, var_lin_comb_chi_monte_carlo = monte_carlo_simulations_lin_comb_chi(
			eigenvalues)
		print('Monte carlo mean', mean_lin_comb_chi_monte_carlo, 'expec6ted mean', 0.1 ** 0.5 * (1 - 1 / 40))
		print('Monte carlo var', var_lin_comb_chi_monte_carlo, 'expec6ted var', 0.1 * (1 / 20))
		
		mean_sqrt_taylor_2 = calculate_cumulant(1, eigenvalues) ** 0.5
		mean_taylor_2_errors.append(
			100 * (mean_sqrt_taylor_2 - mean_lin_comb_chi_monte_carlo) / mean_lin_comb_chi_monte_carlo)
		print('Mean taylor 2 error: ', mean_taylor_2_errors[-1], '%')
		
		mean_sqrt_taylor_3 = calculate_taylor_3_mean_sqrt(eigenvalues)
		mean_taylor_3_errors.append(
			100 * (mean_sqrt_taylor_3 - mean_lin_comb_chi_monte_carlo) / mean_lin_comb_chi_monte_carlo)
		print('Mean taylor 3 error: ', mean_taylor_3_errors[-1], '%')
		
		if abs(mean_taylor_3_errors[-1]) > max_mean_taylor_3_errors:
			max_mean_taylor_3_errors = abs(mean_taylor_3_errors[-1])
			max_mean_taylor_3_eigenvalues = eigenvalues
		
		sum_monte_carlo_variance = sum_monte_carlo_variance + var_lin_comb_chi_monte_carlo
		
		var_sqrt_taylor_2 = calculate_taylor_2_var_sqrt(eigenvalues)
		# print('Var Taylor 2: ', var_sqrt_taylor_2)
		var_taylor_2_errors.append(
			100 * (var_sqrt_taylor_2 - var_lin_comb_chi_monte_carlo) / var_lin_comb_chi_monte_carlo)
		print('Var Taylor 2 error: ', var_taylor_2_errors[-1], '%')
		
		var_sqrt_taylor_3 = calculate_taylor_3_var_sqrt(eigenvalues)
		# print('Taylor 3 variance: ', var_sqrt_taylor_3)
		var_taylor_3_errors.append(
			100 * (var_sqrt_taylor_3 - var_lin_comb_chi_monte_carlo) / var_lin_comb_chi_monte_carlo)
		print('Var Taylor 3  error: ', var_taylor_3_errors[-1], '%')
		
		if abs(var_taylor_3_errors[-1]) > max_var_taylor_3_errors:
			max_var_taylor_3_errors = abs(var_taylor_3_errors[-1])
			max_var_taylor_3_eigenvalues = eigenvalues
		
		print()
		print(sum_monte_carlo_variance / 100)
		print('Final result')
		print('Mean taylor 2 errors: ', numpy.mean(mean_taylor_2_errors), numpy.std(mean_taylor_2_errors),
			  numpy.min(mean_taylor_2_errors),
			  numpy.max(mean_taylor_2_errors))
		print('Mean taylor 3 errors: ', numpy.mean(mean_taylor_3_errors), numpy.std(mean_taylor_3_errors),
			  numpy.min(mean_taylor_3_errors),
			  numpy.max(mean_taylor_3_errors))
		print('Var taylor 2 errors: ', numpy.mean(var_taylor_2_errors), numpy.std(var_taylor_2_errors),
			  numpy.min(var_taylor_2_errors), numpy.max(var_taylor_2_errors))
		print('Var taylor 3 errors: ', numpy.mean(var_taylor_3_errors), numpy.std(var_taylor_3_errors),
			  numpy.min(var_taylor_3_errors), numpy.max(var_taylor_3_errors))
		
		lower_bound_mean_taylor_2_errors = numpy.min(mean_taylor_2_errors)
		upper_bound_mean_taylor_2_errors = numpy.max(mean_taylor_2_errors)
		
		lower_bound_mean_taylor_3_errors = numpy.min(mean_taylor_3_errors)
		upper_bound_mean_taylor_3_errors = numpy.max(mean_taylor_3_errors)
		
		lower_bound_var_taylor_2_errors = numpy.min(var_taylor_2_errors)
		upper_bound_var_taylor_2_errors = numpy.max(var_taylor_2_errors)
		
		lower_bound_var_taylor_3_errors = numpy.min(var_taylor_3_errors)
		upper_bound_var_taylor_3_errors = numpy.max(var_taylor_3_errors)
		
		if mean_extreme != 0:
			mu_Q = eigenvalue_sum / n
			answer = mu_Q ** 0.5
			mean_two_term_extreme_error = 100 * (answer - mean_extreme) / mean_extreme
			
			answer = (1 - 1 / (4 * n)) * mu_Q ** 0.5
			mean_three_term_extreme_error = 100 * (answer - mean_extreme) / mean_extreme
			
			answer = mu_Q / (2 * n)
			var_two_term_extreme_error = 100 * (answer - var_extreme) / var_extreme
			
			answer = (mu_Q / n) * (1 / 2 - 7 / (8 * n) + 3 / (4 * n ** 2))
			var_three_term_extreme_error = 100 * (answer - var_extreme) / var_extreme
		else:
			mean_two_term_extreme_error = 0
			mean_sqrt_taylor_3 = 0
			var_two_term_extreme_error = 0
			var_three_term_extreme_error = 0
		
		dict_item = {'n': n,
					 'mean_mean_taylor_2_errors': numpy.mean(mean_taylor_2_errors),
					 'mean_two_term_extreme_error': mean_two_term_extreme_error,
					 'lower_bound_mean_taylor_2_errors': lower_bound_mean_taylor_2_errors,
					 'upper_bound_mean_taylor_2_errors': upper_bound_mean_taylor_2_errors,
					 'mean_mean_taylor_3_errors': numpy.mean(mean_taylor_3_errors),
					 'mean_three_term_extreme_error': mean_three_term_extreme_error,
					 'lower_bound_mean_taylor_3_errors': lower_bound_mean_taylor_3_errors,
					 'upper_bound_mean_taylor_3_errors': upper_bound_mean_taylor_3_errors,
					 'mean_var_taylor_2_errors': numpy.mean(var_taylor_2_errors),
					 'var_two_term_extreme_error': var_two_term_extreme_error,
					 'lower_bound_var_taylor_2_errors': lower_bound_var_taylor_2_errors,
					 'upper_bound_var_taylor_2_errors': upper_bound_var_taylor_2_errors,
					 'mean_var_taylor_3_errors': numpy.mean(var_taylor_3_errors),
					 'var_three_term_extreme_error': var_three_term_extreme_error,
					 'lower_bound_var_taylor_3_errors': lower_bound_var_taylor_3_errors,
					 'upper_bound_var_taylor_3_errors': upper_bound_var_taylor_3_errors}
		
		results_list.append(dict_item)
	
	field_names = ['n', 'mean_mean_taylor_2_errors', 'mean_two_term_extreme_error', 'lower_bound_mean_taylor_2_errors',
				   'upper_bound_mean_taylor_2_errors',
				   'mean_mean_taylor_3_errors', 'mean_three_term_extreme_error', 'lower_bound_mean_taylor_3_errors',
				   'upper_bound_mean_taylor_3_errors',
				   'mean_var_taylor_2_errors', 'var_two_term_extreme_error', 'lower_bound_var_taylor_2_errors',
				   'upper_bound_var_taylor_2_errors',
				   'mean_var_taylor_3_errors', 'var_three_term_extreme_error', 'lower_bound_var_taylor_3_errors',
				   'upper_bound_var_taylor_3_errors']
	
	with open('taylor_errors_new_' + str(n) + '.csv', 'w', ) as csvfile:
		writer = csv.DictWriter(csvfile, fieldnames=field_names, lineterminator='\n')
		writer.writeheader()
		writer.writerows(results_list)
	
	print(max_mean_taylor_3_eigenvalues)
	print(max_var_taylor_3_eigenvalues)
	# open file in write mode
	
	max_mean_taylor_3_eigenvalues = numpy.around(max_mean_taylor_3_eigenvalues, decimals=3)
	max_var_taylor_3_eigenvalues = numpy.around(max_var_taylor_3_eigenvalues, decimals=3)
	
	results = [max_mean_taylor_3_eigenvalues, max_var_taylor_3_eigenvalues]
	numpy.savetxt('taylor_max_eigenvalues_' + str(n) + '.csv', results)

# read
# result = np.loadtxt("results.txt")
# print(result.tolist())