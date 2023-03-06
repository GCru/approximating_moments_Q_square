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

from bokeh.plotting import figure, output_file, show
from bokeh.models import Range1d, PrintfTickFormatter, NumeralTickFormatter, Title
from bokeh.layouts import row
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
	
	grain = 100
	
	plot0 = figure(width=500, height=500)
	
	
	
	# Maximum
	####################################################
	
	x_axis = [i * 5 * sum(w) / grain for i in range(1, grain + 1)]
	
	#y_axis = [lpb4(coeff=w, x=item) for item in x_axis]
	#plot0.line(x_axis, y_axis, line_width=1)
	
	y_axis = [chi2.cdf(item, 1, scale=1 / 1) for item in x_axis]
	plot0.line(x_axis, y_axis, line_width=2, line_color="blue")

	
	# Inbetween
	#####################################################
	n=2
	if n == 2:

		w=[0,1]
		for idx in range(0,4):
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
			w= [0.70,0.3]
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
	#draw_function()
	#exit()
	
	draw_sum_of_chi_cdf([0, 1],2)
	exit()
	
	