import matplotlib.pyplot as plt
from math import cos, sin, pi
import matplotlib.lines as mlines

r = 3

class restriction_map:
	def __init__ (rm, ref_len):
		rm.ref_l = ref_len

		rm.fig, rm.axes = plt.subplots(num=None, figsize=(9, 9), dpi=80, facecolor='w', edgecolor='k')

		draw_circle = plt.Circle((5, 5), r, fill=False)

		rm.axes.set_aspect(1)
		rm.axes.add_artist(draw_circle)

		rm.axes.set_xlim([0, 10])
		rm.axes.set_ylim([0, 10])

		rm.axes.get_xaxis().set_visible(False)
		rm.axes.get_yaxis().set_visible(False)

		plt.title('Circle')

	def show_map(rm):
		plt.show()

	def add_enzyme(rm, enzyme, x):
		arc = x / rm.ref_l

		#For one line
		x1 = 5 + (r - .25) * cos(arc * 2 * pi)
		y1 = 5 + (r - .25) * sin(arc * 2 * pi)
		x2 = 5 + (r + .25) * cos(arc * 2 * pi)
		y2 = 5 + (r + .25) * sin(arc * 2 * pi)

		x1, y1 = [x1, x2], [y1, y2]

		plt.plot(x1, y1, linestyle='-', linewidth=2, color='blue')

		x1 = 5 + (r + 1) * cos(arc * 2 * pi)
		y1 = 5 + (r + 1) * sin(arc * 2 * pi)

		plt.text(x1, y1, enzyme + ' @ ' + str(x), horizontalalignment='center', verticalalignment='center')
