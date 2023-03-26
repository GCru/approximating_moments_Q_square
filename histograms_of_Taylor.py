import numpy, math
from make_histogram_procedure import make_histogram

from drs import drs

from bokeh.plotting import figure, output_file, show
from bokeh.models import Range1d

from bokeh.plotting import figure, output_file, show
from bokeh.models import Label, LinearAxis, Title, PrintfTickFormatter, NumeralTickFormatter

from bokeh.models.tickers import FixedTicker
from bokeh.models import Arrow, NormalHead, OpenHead, VeeHead, LabelSet
from bokeh.models import Range1d, Div
from bokeh.models import ColumnDataSource
from bokeh.layouts import row, column
from bokeh.io import export_png

from bokeh_constants import *



def draw_error_histogram(error_list, my_title='', max_y=0.2):
	frequency_list, bin_edges_list = make_histogram(error_list, bins=50)
	
	left_edges = bin_edges_list[:-1]
	right_edges = bin_edges_list[1:]
	data = {'frequency': frequency_list, 'left_edges': left_edges, 'right_edges': right_edges}
	source = ColumnDataSource(data=data)
	
	# Set up plot
	plot = figure(width=500, height=500)
	plot.toolbar.logo = None
	plot.toolbar_location = None
	
	
	the_title = Title(text=my_title, align='left',
					 text_font=font, text_font_size="18px",
					 text_line_height=1, vertical_align='middle')
	plot.add_layout(the_title, "above")
	
	plot.x_range = Range1d(left_edges[0], right_edges[-1])
	plot.xaxis.axis_label = ' '
	plot.xaxis[0].formatter = NumeralTickFormatter(format="0%")
	plot.min_border_right = 20
	plot.min_border_bottom = 30
	
	plot.xaxis.axis_label = 'Approximation error'
	plot.yaxis.axis_label = 'Relative Frequency'
	plot.yaxis[0].formatter = NumeralTickFormatter(format="0%")
	plot.y_range = Range1d(0.0, max_y)
	
	plot.axis.axis_label_text_font_size = double_graph_axis_label_font_size
	plot.xaxis.axis_label_text_font = font
	plot.yaxis.axis_label_text_font = font
	
	plot.axis.major_label_text_font_size = double_graph_major_label_font_size
	plot.xaxis.major_label_text_font = font
	plot.yaxis.major_label_text_font = font
	
	
	# construct histogram
	r0 = plot.quad(
		bottom=0, top='frequency', left='left_edges', right='right_edges', source=source,
		fill_color='lightgrey', line_color='black')  # legend='Underperforming monkey portfolios')
	
	return plot


def monte_carlo_simulations_lin_comb_chi(eigenvalues, iterations=10000):
	n = len(eigenvalues)
	
	I = numpy.zeros((n, n))
	numpy.fill_diagonal(I, 1)
	
	chi_list = numpy.empty([iterations])
	
	for count in range(iterations):
		z = numpy.random.default_rng().multivariate_normal(numpy.zeros(n), (1 / n) * I)
		
		chi_list[count] = (numpy.dot(eigenvalues, numpy.square(z))) ** 0.5
	
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

def load_data_and_create_plots(n):
	eigenvalue_sum = 1
	
	# create file with monte carlo expected value and variance
	from pathlib import Path
	
	# fname='monte_carlo-'+str(n)+'.npy'
	fname = Path('monte_carlo_list_' + str(n) + '.npy')
	
	mean_taylor_2_errors = []
	mean_taylor_3_errors = []
	
	var_taylor_2_errors = []
	var_taylor_3_errors = []
	results_list = []
	
	with open(fname, 'rb') as fp:
		counter = 0
		while True:
			try:
				item = numpy.load(fp)
				mean_lin_comb_chi_monte_carlo = item[-2]
				var_lin_comb_chi_monte_carlo = item[-1]
				eigenvalues = item[list(range(0, n))]
				
				print(item)
				counter = counter + 1
				print('Eingevalue set number', counter)
				
				print('Monte carlo mean', mean_lin_comb_chi_monte_carlo, 'expected mean', 0.1 ** 0.5 * (1 - 1 / 40))
				print('Monte carlo var', var_lin_comb_chi_monte_carlo, 'expected var', 0.1 * (1 / 20))
				
				mean_sqrt_taylor_2 = calculate_cumulant(1, eigenvalues) ** 0.5
				mean_taylor_2_errors.append(
					(mean_sqrt_taylor_2 - mean_lin_comb_chi_monte_carlo) / mean_lin_comb_chi_monte_carlo)
				print('Mean taylor 2 error: ', mean_taylor_2_errors[-1], '%')
				
				mean_sqrt_taylor_3 = calculate_taylor_3_mean_sqrt(eigenvalues)
				mean_taylor_3_errors.append(
					(mean_sqrt_taylor_3 - mean_lin_comb_chi_monte_carlo) / mean_lin_comb_chi_monte_carlo)
				print('Mean taylor 3 error: ', mean_taylor_3_errors[-1], '%')
				print('Here')
				
				print('Then here')
				var_sqrt_taylor_2 = calculate_taylor_2_var_sqrt(eigenvalues)
				print('Var Taylor 2: ', var_sqrt_taylor_2)
				var_taylor_2_errors.append(
					(var_sqrt_taylor_2 - var_lin_comb_chi_monte_carlo) / var_lin_comb_chi_monte_carlo)
				print('Var Taylor 2 error: ', var_taylor_2_errors[-1], '%')
				
				var_sqrt_taylor_3 = calculate_taylor_3_var_sqrt(eigenvalues)
				# print('Taylor 3 variance: ', var_sqrt_taylor_3)
				var_taylor_3_errors.append(
					(var_sqrt_taylor_3 - var_lin_comb_chi_monte_carlo) / var_lin_comb_chi_monte_carlo)
				print('Var Taylor 3  error: ', var_taylor_3_errors[-1], '%')
			
			
			
			except:
				print("EoF", counter)
				break
	
	return mean_taylor_3_errors, var_taylor_2_errors



if __name__ == '__main__':
	

	mean_taylor_3_errors_for_n_2, var_taylor_2_errors_for_n_2 = load_data_and_create_plots(2)
	
	plot_mean_for_2 = draw_error_histogram(mean_taylor_3_errors_for_n_2, my_title=r"$$n=2$$",max_y=0.16)
	
	plot_var_for_2 = draw_error_histogram(var_taylor_2_errors_for_n_2, my_title=r"$$n=2$$",max_y=0.09)
	
	mean_taylor_3_errors_for_n_10, var_taylor_2_errors_for_n_10 = load_data_and_create_plots(10)
	
	plot_mean_for_10 = draw_error_histogram(mean_taylor_3_errors_for_n_10,
										   my_title=r"$$n=10$$", max_y=0.16)
	
	
	plot_var_for_10 = draw_error_histogram(var_taylor_2_errors_for_n_10,
										  my_title=r"$$n=10$$",max_y=0.09)
	
	the_row=row(plot_mean_for_2, plot_mean_for_10)
		
	export_plot = column(
		Div(text=r"<h2> &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp Approximation errors: Three-term Taylor expansion for  $${\bf E} \left [\sqrt{Q} \right ]$$</h2>", ),
			the_row)
	
	show(export_plot)
			
	export_png(export_plot, filename="error_histograms_three_term_expectation.png")
	exit()
	
	the_row=row(plot_var_for_2, plot_var_for_10)

	export_plot = column(
			Div(text=r"<h2> &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp  Approximation errors: Two-term Taylor expansion for  $$\mathtt{Var} \left [\sqrt{Q} \right ]$$</h2>", ),
				the_row)
		
	show(export_plot)
		
	export_png(export_plot, filename="error_histograms_two_term_var.png")