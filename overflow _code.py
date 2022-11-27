
# check spread of variances calculated from covariance matrix
# calculate_variance_spread_from_V(V, sigma, rho)
# exit()
examine_standardized_fund_return(V)
exit()


# Note that by making matrix larger does not improve accuracy as the size of V grow, so even
# if the sum of z's approaches zero as n increases


# create covariance matrix V

market_return = numpy.random.default_rng().multivariate_normal(numpy.zeros(n), 1 * V)
# print('dotproduct', numpy.dot(market_return, market_return))

cross_sect_with_z_list = []
cross_sect_with_w_list = []
fund_var_list_z = []
fund_var_list_w = []
standardized_fund_return_list = []

expected_value_with_z_list = []
expected_value_with_w_list = []
cross_z_list = []
for i in range(10):
	# set sigma^2_w = 1/n
	
	z = numpy.random.default_rng().normal(0.0, 1, size=n)
	
	print()
	# Normalise
	# z = z - numpy.mean(z)
	print("Norm:  z*z", numpy.dot(z, z), 'should equal', n)
	
	# create cross z list
	# for j in range(n):
	#	for k in range(j+1,n):
	#		cross_z_list.append(z[j]*z[k])
	
	w = z / sqrt(n)  # same as normal(0.0,1/sqrt(n), size=n)
	print("Norm w*w", numpy.dot(w, w), 'should equal', 1)
	
	fund_var_list_z.append(numpy.dot(z, numpy.dot(V, z)))
	print('Fund risk with z', fund_var_list_z[-1])
	
	fund_var_list_w.append(numpy.dot(w, numpy.dot(V, w)))
	print('Fund risk with w', fund_var_list_w[-1])
	
	# market_return = numpy.random.default_rng().normal(0,sigma,n)
	# market_return = calculate_active_returns(V, w)
	
	standardized_fund_return_list.append(numpy.dot(w, market_return) / sqrt(fund_var_list_w[-1]))
	
	second_term_with_z = sigma * sigma * n * numpy.dot(z, z) / n
	print("Trace term with z", second_term_with_z)
	
	second_term_with_w = sigma * sigma * n * numpy.dot(w, w) / n
	print("Trace term with w", second_term_with_w)
	
	expected_value_with_z_list.append(abs(fund_var_list_z[-1] - second_term_with_z))
	print("Expected value with z", expected_value_with_z_list[-1])
	
	expected_value_with_w_list.append(abs(fund_var_list_w[-1] - second_term_with_w))
	print("Expected value with w", expected_value_with_w_list[-1])
	
	r_active_returns_with_z = calculate_active_returns(V, z)
	# tecnically the active retuns should add to zero
	r_active_returns_with_z = r_active_returns_with_z - numpy.mean(r_active_returns_with_z)
	
	cross_section_with_z = numpy.dot(r_active_returns_with_z, r_active_returns_with_z)
	cross_sect_with_z_list.append(cross_section_with_z)
	
	r_active_returns_with_w = calculate_active_returns(V, w)
	# tecnically the active retuns should add to zero
	r_active_returns_with_w = r_active_returns_with_w - numpy.mean(r_active_returns_with_w)
	
	cross_section_with_w = numpy.dot(r_active_returns_with_w, r_active_returns_with_w)
	cross_sect_with_w_list.append(cross_section_with_w)

print()
# print('St dev', stdev(cross_sect_list), mean(cross_sect_list))

rhs_bound_z = log(n) * log(n) * sqrt(n)
print('Mean expected value with z list', mean(expected_value_with_z_list), 'K value ',
	  mean(expected_value_with_z_list) / rhs_bound_z)

print('Mean fund volatility', mean(fund_var_list_z), 'var with z', stdev(fund_var_list_z))
print('Trace with z', n * sigma * sigma)
print('Fund volatility equals cross sect with z ', sigma * sigma * (1 - rho) * n)
print('Cross section calculated with z ', mean(cross_sect_with_z_list))

print()
rhs_bound_w = log(n) * log(n) / sqrt(n)
error_vs_trace = [abs(fund_var_list_w[idx] - sigma * sigma) for idx in range(len(fund_var_list_w))]

print('Mean expected value with w list', mean(expected_value_with_w_list), 'K value ',
	  mean(expected_value_with_w_list) / rhs_bound_w)
print('Mean fund volatility', mean(fund_var_list_w),
	  'var with w', stdev(fund_var_list_w))
print('Trace with w', sigma * sigma)
print('Fund volatility equals cross sect with w ', sigma * sigma * (1 - rho))
print('Cross section calculated with w ', mean(cross_sect_with_w_list))

print('Mean_error', mean(error_vs_trace))

print("norm matrix", sigma * sigma + n * rho * sigma * sigma)
# print('Mean cross z list', mean(cross_z_list), len(cross_z_list), n*(n-1))

print('Stdev van standarzidez fund', stdev(standardized_fund_return_list))

# draw histogram of fund_stdev_list
from shared_with_active_risk_paper import make_histogram

real_stdev_list = [sqrt(fund_var_list_w[idx]) for idx in range(len(fund_var_list_w))]
frequency_list, bin_edges_list = make_histogram(real_stdev_list)
print(stdev(real_stdev_list), median_abs_deviation(real_stdev_list))

left_edges = bin_edges_list[:-1]
right_edges = bin_edges_list[1:]
data = {'frequency': frequency_list, 'left_edges': left_edges, 'right_edges': right_edges}
source = ColumnDataSource(data=data)

# Set up plot
plot0 = figure(plot_width=int(500), plot_height=500)
r0 = plot0.quad(
	bottom=0, top='frequency', left='left_edges', right='right_edges', source=source,
	fill_color='lightgrey', line_color='black')  # legend='Underperforming monkey portfolios')

show(plot0)

# draw histogram of inverse
from shared_with_active_risk_paper import make_histogram

frequency_list, bin_edges_list = make_histogram(standardized_fund_return_list)

left_edges = bin_edges_list[:-1]
right_edges = bin_edges_list[1:]
data = {'frequency': frequency_list, 'left_edges': left_edges, 'right_edges': right_edges}
source = ColumnDataSource(data=data)

# Set up plot
plot1 = figure(plot_width=int(500), plot_height=500)
r0 = plot1.quad(
	bottom=0, top='frequency', left='left_edges', right='right_edges', source=source,
	fill_color='lightgrey', line_color='black')  # legend='Underperforming monkey portfolios')
my_plot = row(plot0, plot1)
# show(plot1)