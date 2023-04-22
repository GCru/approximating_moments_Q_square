import numpy
import drs


def update_monte_carlo_file(n, fname, total_eigenvalue_draws, iterations_per_eigenvalue=10000):

	if fname.is_file():
		print('Appending to an existing file')
		with open(fname, 'rb') as fp:
			num_lines = 0
			while True:
				try:
					item = numpy.load(fp)
					print(item)
					num_lines = num_lines + 1
				except:
					print("EoF", num_lines)
					break
		
		# num_lines = sum(1 for line in fp)
		# print('Total lines:', num_lines)  # 8
		# for line in fp:
		#	print(line)
		total_remaining_eigenvalue_draws = total_eigenvalue_draws - num_lines
	else:
		print('Starting a new file')
		num_lines = 0
		total_remaining_eigenvalue_draws = total_eigenvalue_draws
	
	for count_eigenvalue_draws in range(total_remaining_eigenvalue_draws):
		eigenvalues = drs(n, 1)
		
		I = numpy.zeros((n, n))
		numpy.fill_diagonal(I, 1)
		
		chi_list = numpy.empty([iterations_per_eigenvalue])
		
		for count in range(iterations_per_eigenvalue):
			z = numpy.random.default_rng().multivariate_normal(numpy.zeros(n), (1 / n) * I)
			# print('stuck', count)
			chi_list[count] = (numpy.dot(eigenvalues, numpy.square(z))) ** 0.5
		
		mean_lin_comb_chi_monte_carlo = numpy.mean(chi_list)
		var_lin_comb_chi_monte_carlo = numpy.var(chi_list)
		
		print('Eigenvalue set number: ', num_lines + count_eigenvalue_draws + 1)
		print('Eigenvalue set:', eigenvalues)
		print('Mean:', mean_lin_comb_chi_monte_carlo, 'Var:', var_lin_comb_chi_monte_carlo)
		
		output_array = numpy.concatenate(
			(eigenvalues, numpy.array([mean_lin_comb_chi_monte_carlo, var_lin_comb_chi_monte_carlo])))
		with open(fname, 'ab') as f:
			numpy.save(f, output_array)
	
	return


def update_monte_carlo_file(n, fname, total_eigenvalue_draws, iterations_per_eigenvalue=10000):
	if fname.is_file():
		print('Appending to an existing file')
		with open(fname, 'rb') as fp:
			num_lines = 0
			while True:
				try:
					item = numpy.load(fp)
					print(item)
					num_lines = num_lines + 1
				except:
					print("EoF", num_lines)
					break
		
		# num_lines = sum(1 for line in fp)
		# print('Total lines:', num_lines)  # 8
		# for line in fp:
		#	print(line)
		total_remaining_eigenvalue_draws = total_eigenvalue_draws - num_lines
	else:
		print('Starting a new file')
		num_lines = 0
		total_remaining_eigenvalue_draws = total_eigenvalue_draws
	
	for count_eigenvalue_draws in range(total_remaining_eigenvalue_draws):
		eigenvalues = drs(n, 1)
		
		I = numpy.zeros((n, n))
		numpy.fill_diagonal(I, 1)
		
		chi_list = numpy.empty([iterations_per_eigenvalue])
		
		for count in range(iterations_per_eigenvalue):
			z = numpy.random.default_rng().multivariate_normal(numpy.zeros(n), (1 / n) * I)
			# print('stuck', count)
			chi_list[count] = (numpy.dot(eigenvalues, numpy.square(z))) ** 0.5
		
		mean_lin_comb_chi_monte_carlo = numpy.mean(chi_list)
		var_lin_comb_chi_monte_carlo = numpy.var(chi_list)
		
		print('Eigenvalue set number: ', num_lines + count_eigenvalue_draws + 1)
		print('Eigenvalue set:', eigenvalues)
		print('Mean:', mean_lin_comb_chi_monte_carlo, 'Var:', var_lin_comb_chi_monte_carlo)
		
		output_array = numpy.concatenate(
			(eigenvalues, numpy.array([mean_lin_comb_chi_monte_carlo, var_lin_comb_chi_monte_carlo])))
		with open(fname, 'ab') as f:
			numpy.save(f, output_array)
	
	return

if __name__ == '__main__':
	
	n=10
	I = numpy.zeros((n, n))
	numpy.fill_diagonal(I, 1)
	print('Start mulitvariate normal')
	for count in range(100):
		z = numpy.random.default_rng().multivariate_normal(numpy.zeros(n), (1 / n) * I)
	print('End multivariate normal', numpy.mean(z))
	
	print('Start singlevariate normal')
	for count in range(1):
		z = numpy.random.default_rng().normal(size=n)
	z=z/n
	print('End multivariate normal', z)
	
	