import numpy, math
from datetime import date
from statistics import stdev, mean
from math import pi, log, sqrt
from scipy.stats import kstest, kurtosis, skew, skewtest, jarque_bera, wasserstein_distance, median_abs_deviation, chi2
from scipy.special import erf
from scipy import optimize

from momentchi2 import hbe, lpb4,sw,wf
from drs import drs
from imhoff import imhoff

import bokeh_constants
from bokeh.plotting import figure, output_file, show
from bokeh.models import Range1d, PrintfTickFormatter, NumeralTickFormatter, Title
from bokeh.layouts import row
from bokeh.layouts import row,column, Spacer
from bokeh.models import Range1d, Div
from bokeh.models import ColumnDataSource

"""from fundskill_utilities.fundskill_utilities import change_to_another_month_end, \
	change_to_month_end

from latex_label_utility import LatexLabel

from fundskill_utilities.fundskill_shared_bokeh_utilities import setup_date_ticker_jump

from bokeh_constants import *"""


def a_h(w):

	n=len(w)
	sum=numpy.sum(numpy.multiply(numpy.multiply(w,w),w))
	
	sum = sum/ (numpy.dot(w,w)*n)
	
	print('a_h',sum)
	
	return sum

def n_h(w):
	
	d=numpy.sum(numpy.multiply(numpy.multiply(w,w),w))
	sum = numpy.dot(w,w)**3
	sum=sum/ d**2
	
	print('n_h', sum)
	return sum

def beta_h(w):
	
	
	return n_h(w) - (2*n_h(w)-1)**0.5


def draw_sum_of_chi_cdf(w, n):
	
	grain = 900 #( note 1000 does not work)
	#plot0 = figure(height=250, title=r"\[\sin(x)\text{ for }x\text{ between }-2\pi\text{ and }2\pi\]")
	plot0 = figure(width=500, height=500, title = r"$$\chi^2$$" , tools="",
           toolbar_location=None)
	
	if n==2:
		plot0.title = r"$$\text{CDFs of weighted sum of }n=2\text{ chi-squared variables}$$"
	if n==10:
		plot0.title = r"$$\text{CDFs of weighted sum of }n=10\text{ chi-squared variables}$$"
	
	plot0.title.text_font_size="16px" #bokeh_constants.double_graph_sub_title_font_size
	
	plot0.xaxis.axis_label = r"$$Q$$"
	plot0.xaxis.axis_label_text_font_size = bokeh_constants.double_graph_axis_label_font_size
	
	plot0.yaxis.axis_label = r"$$F_{Q}$$"
	plot0.yaxis.axis_label_text_font_size = bokeh_constants.double_graph_axis_label_font_size
	
	plot0.axis.major_label_text_font_size = bokeh_constants.double_graph_major_label_font_size
	
	# Maximum
	####################################################
	
	x_axis = [i * 5 * sum(w) / grain for i in range(1, grain + 1)]
	plot0.x_range = Range1d(0,5)
	
	#y_axis = [lpb4(coeff=w, x=item) for item in x_axis]
	#plot0.line(x_axis, y_axis, line_width=1)
	
	y_axis = [chi2.cdf(item, 1, scale=1 / 1) for item in x_axis]
	plot0.line(x_axis, y_axis, line_width=2, line_color="blue")
	plot0.y_range = Range1d(0,1.01)

	
	# Inbetween
	#####################################################
	#n=2
	if n == 2:

		w=[0,1]
		for idx in range(0,1):
			w[0] = 0.1*(idx+1)
			w[1] = 1-w[0]
			print(w)
			#w = drs(n, 1)
			x_axis = [i * 5 * sum(w) / grain for i in range(1, grain + 1)]
	
			#y_axis = [lpb4(coeff=w, x=item, p=4) for item in x_axis]
		
			y_axis = [(1-imhoff(item, w, [1] * len(w), [0] * len(w))) for item in x_axis]
			plot0.line(x_axis, y_axis, line_width=1, line_color="black", line_dash='dashed')
			
	else:
		
		for idx in range(0, 1):
			w = drs(n, 1)
			#w= [0.70,0.3]
			print(w)
			
			x_axis = [i * 5 * sum(w) / grain for i in range(1, grain + 1)]
			
			# y_axis = [lpb4(coeff=w, x=item, p=4) for item in x_axis]
			
			y_axis = [(1 - imhoff(item, w, [1] * len(w), [0] * len(w))) for item in x_axis]
			plot0.line(x_axis, y_axis, line_width=1, line_color="black", line_dash='dashed')
		
	# Minimum
	##########################################################

	w = [1/n]*n
	
	x_axis = [i * 5 * sum(w) / grain for i in range(1, grain + 1)]
	
	#y_axis = [lpb4(coeff=w, x=item) for item in x_axis]
	#plot0.line(x_axis, y_axis, line_width=2, line_color="red")
	
	y_axis = [chi2.cdf(item,n,scale=1/n) for item in x_axis]
	plot0.line(x_axis, y_axis, line_width=2, line_color="red")
	
	return plot0


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
	
	#show(plot0)
	return


if __name__ == '__main__':
	#draw_function()
	#exit()
	
	plot_left = draw_sum_of_chi_cdf([0, 1],2)
	
	plot_right = draw_sum_of_chi_cdf([0, 1],10)
	
	the_row = row(plot_left, Spacer(width=15), plot_right)
	show(the_row)
	
	
	