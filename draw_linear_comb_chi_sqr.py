import numpy, math
from datetime import date
from statistics import stdev, mean
from math import pi, log, sqrt
from scipy.stats import kstest, kurtosis, skew, skewtest, jarque_bera, wasserstein_distance, median_abs_deviation, chi2
from scipy.special import erf
from scipy import optimize

from momentchi2 import hbe

from bokeh.plotting import figure, output_file, show
from bokeh.models import Range1d, PrintfTickFormatter, NumeralTickFormatter, Title
from bokeh.layouts import row
from bokeh.models import ColumnDataSource

"""from fundskill_utilities.fundskill_utilities import change_to_another_month_end, \
	change_to_month_end

from latex_label_utility import LatexLabel

from fundskill_utilities.fundskill_shared_bokeh_utilities import setup_date_ticker_jump

from bokeh_constants import *"""


def draw_sum_of_chi_cdf(w):
	
	grain = 100
	
	plot0 = figure(plot_width=int(500), plot_height=500)
	
	
	
	w = [0.01 / 3, 0.01 / 3, 0.01 / 3, 0.999]
	# w=[0.01,0.01,0.01,0.97]
	# w=[0.2,0.2,0.2,0.2,0.2]
	# w=[0.40,0.2,0.3,0.1]
	print('*****', hbe(coeff=w, x=1))
	
	# w=[0.499,0.001,0.499,0.001]
	print(sum(w))
	x_axis = [i * 4 * sum(w) / grain for i in range(1, grain + 1)]
	
	y_axis = [hbe(coeff=w, x=item) for item in x_axis]
	plot0.line(x_axis, y_axis, line_width=2)
	
	y_axis = [chi2.cdf(item, 1, scale=1 / 1) for item in x_axis]
	plot0.line(x_axis, y_axis, line_width=2, line_color="black")
	#show(plot0)
	#exit()
	
	w = [0.106, 0.106, 0.106, 0.682]
	w = [0.25, 0.25, 0.25, 0.25]
	# w=[0.02,0.02,0.01,0.95]
	print(sum(w))
	x_axis = [i * 4 * sum(w) / grain for i in range(1, grain + 1)]
	
	y_axis = [hbe(coeff=w, x=item) for item in x_axis]
	plot0.line(x_axis, y_axis, line_width=2, line_color="red")
	
	y_axis = [chi2.cdf(item,100,scale=1/100) for item in x_axis]
	plot0.line(x_axis, y_axis, line_width=2, line_color="black")
	
	
	show(plot0)
	
	return


def draw_function():

	grain = 100
	
	plot0 = figure(plot_width=int(500), plot_height=500)
	
	x_min=0
	x_max=10
	
	x_axis = [i*(x_max-x_min)/ grain for i in range(1, grain + 1)]
	
	y_axis = [math.exp(-item) for item in x_axis]
	plot0.line(x_axis, y_axis, line_width=2)
	
	y_axis= [math.exp(-item/10) for item in x_axis]
	plot0.line(x_axis, y_axis, line_width=2,line_color='black')
	
	show(plot0)
	return

if __name__ == '__main__':
	draw_function()
	exit()
	
	draw_sum_of_chi_cdf([0, 1])
	exit()
	
	