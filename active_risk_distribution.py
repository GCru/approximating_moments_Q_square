





if __name__ == '__main__':
	
	date_range = set_up_in_sample_month_end_list()
	
	# date_range =[start_date]
	
	# Set up benchmarks of the funds to consider
	benchmark_list = ['Russell Mid Cap TR USD', 'Russell Mid Cap TR USD', 'Russell Mid Cap Growth TR US', ]
	benchmark_list = ['Russell 3000 TR USD']
	benchmark_list = ['Russell 2000 TR USD']
	benchmark_list = ['Russell 1000 TR USD', 'Russell 1000 Growth TR USD', 'Russell 1000 Value TR USD', ]
	benchmark_list = ['Russell 2000 Growth TR USD', 'Russell 2000 Value TR USD', ]
	
	benchmark_list = ['S&P 500 TR USD']
	
	for benchmark_name in benchmark_list:
		results_file_name = 'random_10_10000_funds_risk_forecast_' + benchmark_name + '.csv'
		
		calculate_small_random_fund_results_and_write_to_file(
			benchmark_name, date_range, results_file_name, number_of_months=60, number_of_stocks=10)