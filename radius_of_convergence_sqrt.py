import numpy

from bokeh.plotting import figure, show
from bokeh.layouts import row,column, Spacer
from bokeh.models import Range1d, Div
from bokeh.models import Legend, LegendItem
from bokeh.models.annotations.labels import Label
from bokeh.io import export_png

import bokeh_constants
from math import pi
x = numpy.linspace(0, 6, 600)
y = x**0.5

#x=2.72
x0=1.0  # expansion point)
taylor_1 = x0**0.5
taylor_2 = taylor_1+0.5*x0**(-0.5)*(x-x0)
taylor_3 = taylor_2-(1/8)*x0**(-3/2)*(x-x0)**2
taylor_4 = taylor_3+(3/48)*x0**(-5/2)*(x-x0)**3
taylor_5 = taylor_4-(15/(16*24))*x0**(-7/2)*(x-x0)**4
taylor_6 = taylor_5+ (105/(32*120))*x0**(-9/2)*(x-x0)**5

# check which Tayulor is best, by settin x = to various values above
#print(abs(taylor_1-x**0.5), abs(taylor_2-x**0.5),abs(taylor_3-x**0.5), abs(taylor_4-x**0.5))

p_left = figure(width=500, height=500, tools="",toolbar_location=None)
p_left.x_range=Range1d(0,5)
p_left.y_range=Range1d(0,3)

p_sqrt=p_left.line(x, y, color="black", alpha=1, line_width=3) # legend_label=r"$$\sqrt{x} $$")
p_left.line(x, taylor_2, color="black", alpha=0.6, line_width=2, line_dash="dashed", legend_label="Two-term Taylor")
p_left.line(x, taylor_3, color="black", alpha=0.6, line_width=2, line_dash="solid", legend_label="Three-term Taylor")
p_left.line(x, taylor_4, color="black", alpha=0.6, line_width=2, line_dash='dotted',legend_label="Four-term Taylor")

label_left = Label(
			text=r"$$ \sqrt{x}$$", x=4.4, y=2.25,text_font_size='22px' )
p_left.add_layout(label_left)


#legend =Legend(items=[LegendItem(label=Label(text=r"$$\sqrt{x}$$"),renderers=[p_sqrt],index=0)])

#p_left.add_layout(legend)
p_left.legend.location="top_left"
p_left.legend.label_text_font_size='22px' #bokeh_constants.double_graph_major_label_font_size

p_left.xaxis.axis_label = r"$$x$$"
p_left.xaxis.axis_label_text_font_size='22px'  #bokeh_constants.double_graph_major_label_font_size
p_left.axis.major_label_text_font_size='20px'  #bokeh_constants.double_graph_major_label_font_size


#show(p_left)
#exit()
p_right = figure(width=500, height=500, title="", tools="",
           toolbar_location=None)

p_right.x_range=Range1d(0,6)
p_right.y_range=Range1d(-5,25)

p_right.line(x, y, color="black", alpha=1, line_width=3)
p_right.line(x, taylor_3, color="black", alpha=0.6, line_width=2, line_dash="solid", legend_label="Three-term Taylor")
p_right.line(x, taylor_4, color="black", alpha=0.6, line_width=2, line_dash='dotted', legend_label= "Four-term Taylor")
p_right.line(x, taylor_5, color="black", alpha=0.6, line_width=2, line_dash='dotdash',legend_label= "Five-term Taylor")
p_right.line(x, taylor_6, color="black", alpha=0.6, line_width=2, line_dash='dashed', legend_label= "Six-term Taylor")

label_right = Label(
			text=r"$$ \sqrt{x}$$", x=5.5, y=3, text_font_size='22px')

p_right.add_layout(label_right)

p_right.legend.location="top_left"
p_right.legend.label_text_font_size='22px' #bokeh_constants.double_graph_major_label_font_size

p_right.xaxis.axis_label = r"$$x$$"
p_right.xaxis.axis_label_text_font_size='22px' #bokeh_constants.double_graph_major_label_font_size
p_right.axis.major_label_text_font_size='20px' #bokeh_constants.double_graph_major_label_font_size

the_row = row(p_left,Spacer(width=15),p_right)
#p_left.background_fill_color = "#efefef"
#p.xaxis.fixed_location = 0
#p.yaxis.fixed_location = 0

export_plot = column(Div(text=r"<h1>&nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp$$\text{Approximating } \sqrt{x} \text{ with truncated Taylor expansions}$$</h1>",), the_row)

show(export_plot)

export_png(export_plot, filename="truncated_taylor_expansions_for_root_x.png")