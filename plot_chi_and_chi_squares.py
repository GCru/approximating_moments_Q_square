import numpy, math
from datetime import date
from statistics import stdev, mean
from math import pi, log, sqrt
from scipy.stats import kstest, kurtosis, skew, skewtest, jarque_bera, wasserstein_distance, median_abs_deviation, chi2,chi
from scipy.special import erf
from scipy import optimize

from momentchi2 import hbe, lpb4, sw, wf
from drs import drs
from imhoff import imhoff

import bokeh_constants
from bokeh.plotting import figure, output_file, show
from bokeh.models import Range1d, PrintfTickFormatter, NumeralTickFormatter, Title
from bokeh.layouts import row
from bokeh.layouts import row, column, Spacer
from bokeh.models import Range1d, Div
from bokeh.models.annotations.labels import Label
from bokeh.io import export_png
from bokeh.models import ColumnDataSource

"""from fundskill_utilities.fundskill_utilities import change_to_another_month_end, \
	change_to_month_end

from latex_label_utility import LatexLabel

from fundskill_utilities.fundskill_shared_bokeh_utilities import setup_date_ticker_jump

from bokeh_constants import *"""


def a_h(w):
	n = len(w)
	sum = numpy.sum(numpy.multiply(numpy.multiply(w, w), w))
	
	sum = sum / (numpy.dot(w, w) * n)
	
	print('a_h', sum)
	
	return sum


def n_h(w):
	d = numpy.sum(numpy.multiply(numpy.multiply(w, w), w))
	sum = numpy.dot(w, w) ** 3
	sum = sum / d ** 2
	
	print('n_h', sum)
	return sum


def beta_h(w):
	return n_h(w) - (2 * n_h(w) - 1) ** 0.5


def draw_chi_square_cdf(n):
	grain = 900  # ( note 1000 does not work)
	# plot0 = figure(height=250, title=r"\[\sin(x)\text{ for }x\text{ between }-2\pi\text{ and }2\pi\]")
	plot0 = figure(width=500, height=500, title=r"$$\chi^2$$", tools="",
				   toolbar_location=None)
	
	
	plot0.title = r"$$n=10$$"
	
	plot0.title.text_font_size = "22px"  # bokeh_constants.double_graph_sub_title_font_size
	
	plot0.xaxis.axis_label = r"$$Q$$"
	plot0.xaxis.axis_label_text_font_size = '22px'  # bokeh_constants.double_graph_axis_label_font_size
	
	plot0.yaxis.axis_label = r"$$F_{Q}$$"
	plot0.yaxis.axis_label_text_font_size = '22px'  # bokeh_constants.double_graph_axis_label_font_size
	
	plot0.axis.major_label_text_font_size = '20px'  # bokeh_constants.double_graph_major_label_font_size
	
	# Maximum
	####################################################
	
	x_axis = [i * 5 / grain for i in range(1, grain + 1)]
	plot0.x_range = Range1d(0, 5)
	
	y1_axis = [chi2.cdf(item, 1, scale=1 / 1) for item in x_axis]
	plot0.line(x_axis, y1_axis, line_width=2, line_color="black", alpha=1)
	plot0.y_range = Range1d(0, 1.002)
	
	label = Label(
		text=r"$$\chi^2(1)$$", x=3.23, y=0.84, text_font_size='22px')
	
	plot0.add_layout(label)
	
	
	# Minimum
	##########################################################
	
	x_axis = [i * 5  / grain for i in range(1, grain + 1)]
	
	y2_axis = [chi2.cdf(item, n, scale=1 / n) for item in x_axis]
	plot0.line(x_axis, y2_axis, line_width=0.8, line_color="black", alpha=0.8)
	
	label = Label(
			text=r"$${\small \frac{1}{10}}\chi^2(10)$$", x=0.9, y=0.25, text_font_size='22px')
	
	for i in range (len(x_axis)):
		if x_axis[i]>2:
			break
			
	
	plot0.varea(x=x_axis[i:], y1=y1_axis[i:], y2=y2_axis[i:])
	
	#plot0.varea(x=x_axis[:20], y1=y1_axis[:20], y2=y2_axis[:20])
	
	plot0.add_layout(label)
	
	return plot0


def draw_chi_square_pdf(n):
	grain = 900  # ( note 1000 does not work)
	# plot0 = figure(height=250, title=r"\[\sin(x)\text{ for }x\text{ between }-2\pi\text{ and }2\pi\]")
	plot0 = figure(width=500, height=500, title=r"$$\chi^2$$", tools="",
				   toolbar_location=None)
	
	plot0.title = r"$$n=10$$"
	
	plot0.title.text_font_size = "22px"  # bokeh_constants.double_graph_sub_title_font_size
	
	plot0.xaxis.axis_label = r"$$Q$$"
	plot0.xaxis.axis_label_text_font_size = '22px'  # bokeh_constants.double_graph_axis_label_font_size
	
	plot0.yaxis.axis_label = r"$$F_{Q}$$"
	plot0.yaxis.axis_label_text_font_size = '22px'  # bokeh_constants.double_graph_axis_label_font_size
	
	plot0.axis.major_label_text_font_size = '20px'  # bokeh_constants.double_graph_major_label_font_size
	
	# Maximum
	####################################################
	
	x_axis = [i * 5 / grain for i in range(1, grain + 1)]
	plot0.x_range = Range1d(0, 5)
	
	y1_axis = [chi2.pdf(item, 1, scale=1 / 1) for item in x_axis]
	plot0.line(x_axis, y1_axis, line_width=2, line_color="black", alpha=1)
	plot0.y_range = Range1d(0, 1.602)
	
	label = Label(
		text=r"$$\chi^2(1)$$", x=3.23, y=0.84, text_font_size='22px')
	
	plot0.add_layout(label)
	
	# Minimum
	##########################################################
	
	x_axis = [i * 5 / grain for i in range(1, grain + 1)]
	
	y2_axis = [chi2.pdf(item, n, scale=1 / n) for item in x_axis]
	plot0.line(x_axis, y2_axis, line_width=0.8, line_color="black", alpha=0.8)
	
	label = Label(
		text=r"$${\small \frac{1}{10}}\chi^2(10)$$", x=0.9, y=0.25, text_font_size='22px')
	
	for i in range(len(x_axis)):
		if x_axis[i] > 2:
			break
	
	plot0.varea(x=x_axis[i:], y1=y1_axis[i:], y2=y2_axis[i:])
	
	# plot0.varea(x=x_axis[:20], y1=y1_axis[:20], y2=y2_axis[:20])
	
	plot0.add_layout(label)
	
	return plot0


def draw_chi_cdf(n):
	grain = 900  # ( note 1000 does not work)
	# plot0 = figure(height=250, title=r"\[\sin(x)\text{ for }x\text{ between }-2\pi\text{ and }2\pi\]")
	plot0 = figure(width=500, height=500, title=r"$$\chi^2$$", tools="",
				   toolbar_location=None)
	
	plot0.title = r"$$n=10$$"
	
	plot0.title.text_font_size = "22px"  # bokeh_constants.double_graph_sub_title_font_size
	
	plot0.xaxis.axis_label = r"$$Q$$"
	plot0.xaxis.axis_label_text_font_size = '22px'  # bokeh_constants.double_graph_axis_label_font_size
	
	plot0.yaxis.axis_label = r"$$F_{Q}$$"
	plot0.yaxis.axis_label_text_font_size = '22px'  # bokeh_constants.double_graph_axis_label_font_size
	
	plot0.axis.major_label_text_font_size = '20px'  # bokeh_constants.double_graph_major_label_font_size
	
	# Maximum
	####################################################
	
	x_axis = [i * 5 / grain for i in range(1, grain + 1)]
	plot0.x_range = Range1d(0, 3)
	
	y1_axis = [chi.cdf(item, 1, scale=1 / 1) for item in x_axis]
	plot0.line(x_axis, y1_axis, line_width=2, line_color="black", alpha=1)
	plot0.y_range = Range1d(0, 1.002)
	
	label = Label(
		text=r"$$\chi^2(1)$$", x=3.23, y=0.84, text_font_size='22px')
	
	plot0.add_layout(label)
	
	# Minimum
	##########################################################
	
	x_axis = [i * 5 / grain for i in range(1, grain + 1)]
	
	y2_axis = [chi.cdf(item, n, scale=1 / n**0.5) for item in x_axis]
	plot0.line(x_axis, y2_axis, line_width=0.8, line_color="black", alpha=0.8)
	
	label = Label(
		text=r"$${\small \frac{1}{10}}\chi^2(10)$$", x=0.9, y=0.25, text_font_size='22px')
	
	for i in range(len(x_axis)):
		if x_axis[i] > 2:
			break
	
	plot0.varea(x=x_axis[i:], y1=y1_axis[i:], y2=y2_axis[i:])
	
	# plot0.varea(x=x_axis[:20], y1=y1_axis[:20], y2=y2_axis[:20])
	
	plot0.add_layout(label)
	
	return plot0


def draw_chi_pdf(n):
	grain = 900  # ( note 1000 does not work)
	# plot0 = figure(height=250, title=r"\[\sin(x)\text{ for }x\text{ between }-2\pi\text{ and }2\pi\]")
	plot0 = figure(width=500, height=500, title=r"$$\chi^2$$", tools="",
				   toolbar_location=None)
	
	plot0.title = r"$$n=10$$"
	
	plot0.title.text_font_size = "22px"  # bokeh_constants.double_graph_sub_title_font_size
	
	plot0.xaxis.axis_label = r"$$Q$$"
	plot0.xaxis.axis_label_text_font_size = '22px'  # bokeh_constants.double_graph_axis_label_font_size
	
	plot0.yaxis.axis_label = r"$$F_{Q}$$"
	plot0.yaxis.axis_label_text_font_size = '22px'  # bokeh_constants.double_graph_axis_label_font_size
	
	plot0.axis.major_label_text_font_size = '20px'  # bokeh_constants.double_graph_major_label_font_size
	
	# Maximum
	####################################################
	
	x_axis = [i * 5 / grain for i in range(1, grain + 1)]
	plot0.x_range = Range1d(0, 3)
	
	y1_axis = [chi.pdf(item, 1, scale=1 / 1) for item in x_axis]
	plot0.line(x_axis, y1_axis, line_width=2, line_color="black", alpha=1)
	plot0.y_range = Range1d(0, 2)
	
	label = Label(
		text=r"$$\chi^2(1)$$", x=3.23, y=0.84, text_font_size='22px')
	
	plot0.add_layout(label)
	
	# Minimum
	##########################################################
	
	x_axis = [i * 5 / grain for i in range(1, grain + 1)]
	
	y2_axis = [chi.pdf(item, n, scale=1 / n ** 0.5) for item in x_axis]
	plot0.line(x_axis, y2_axis, line_width=0.8, line_color="black", alpha=0.8)
	
	label = Label(
		text=r"$${\small \frac{1}{10}}\chi^2(10)$$", x=0.9, y=0.25, text_font_size='22px')
	
	for i in range(len(x_axis)):
		if x_axis[i] > 2:
			break
	
	plot0.varea(x=x_axis[i:], y1=y1_axis[i:], y2=y2_axis[i:])
	
	# plot0.varea(x=x_axis[:20], y1=y1_axis[:20], y2=y2_axis[:20])
	
	plot0.add_layout(label)
	
	return plot0


def draw_function():
	grain = 100
	
	plot0 = figure(plot_width=int(500), plot_height=500)
	
	x_min = 0
	x_max = 10
	
	x_axis = [i * (x_max - x_min) / grain for i in range(1, grain + 1)]
	
	y_axis = [math.exp(-item) for item in x_axis]
	plot0.line(x_axis, y_axis, line_width=2)
	
	y_axis = [math.exp(-item / 10) for item in x_axis]
	plot0.line(x_axis, y_axis, line_width=2, line_color='black')
	
	# show(plot0)
	return


if __name__ == '__main__':
	# draw_function()
	# exit()
	
	plot_left = draw_chi_cdf(10)
	
	plot_right = draw_chi_pdf(10)
	
	
	the_row = row(plot_left, Spacer(width=15), plot_right)
	
	show(the_row)
	exit()
	
	plot_left = draw_chi_square_cdf(10)
	
	plot_right = draw_chi_square_pdf(10)
	
	#show(plot_left)
	
	#exit()
	
	
	the_row = row(plot_left, Spacer(width=15), plot_right)
	
	show(the_row)
	
	exit()
	
	export_plot = column(
		Div(text=r"<h1>&nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp $$\text{CDFs of a weighted sum of chi-squared variables}$$</h1>", ),
		the_row)
	
	show(export_plot)
	
	export_png(export_plot, filename="CDFs_weighted_chi_squares.png")