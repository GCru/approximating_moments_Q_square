import numpy

from bokeh.plotting import figure, show
from bokeh.models import Range1d

x = numpy.linspace(0, 5, 500)
y = x**0.5

x0=1.0  # expansion point)
taylor_1 = x0**0.5
taylor_2 = taylor_1+0.5*x0**(-0.5)*(x-x0)
taylor_3 = taylor_2-(1/8)*x0**(-3/2)*(x-x0)**2
taylor_4 = taylor_3+(3/48)*x0**(-5/2)*(x-x0)**3
taylor_5 = taylor_4-(15/(16*24))*x0**(-7/2)*(x-x0)**4
taylor_6 = taylor_5+ (105/(32*120))*x0**(-9/2)*(x-x0)**5

p = figure(width=800, height=350, title="", tools="",
           toolbar_location=None)

p.x_range=Range1d(0,4)
p.y_range=Range1d(0,3.2)

p.line(x, y, color="navy", alpha=0.4, line_width=4)
#p.line(x, taylor_1, color="navy", alpha=0.4, line_width=2)
p.line(x, taylor_2, color="navy", alpha=0.4, line_width=2)
p.line(x, taylor_3, color="navy", alpha=0.4, line_width=2)
p.line(x, taylor_4, color="navy", alpha=0.4, line_width=2)
#p.line(x, taylor_5, color="navy", alpha=0.4, line_width=2)
#p.line(x, taylor_6, color="navy", alpha=0.4, line_width=2)



p.background_fill_color = "#efefef"
#p.xaxis.fixed_location = 0
#p.yaxis.fixed_location = 0

show(p)