from PIL import Image
from matplotlib import pyplot as plot

import matplotlib.animation as animation
import io
import numpy as np



def save_animation(data, save_path):

	# each sublist should be a list of images; so data should be a 4D array [subplot_series, timepoint, img_width, img_height].
	data = np.asarray(data)
	data = data.expand_dims() if len(data.shape)==3 else data

	num_subplots = data.shape[0]

	# animation update functio.
	def update(i, ax_idx, fig):
		plot.subplot(1, num_subplots, ax_idx, figure=fig)
		im.set_array(data[ax_idx][0])
		ax.set_title(f"time={i}")
		return im, 

	fig = plot.figure()
	for ax_idx in range(num_subplots):
		ax = plot.subplot(1, num_subplots, ax_idx+1)
		# set the initial image.
		im = ax.imshow(data[ax_idx][0], animated=True, vmin=0, vmax=10)
		ax.set_title(f"time={1}")
		fig.colorbar(im)
		# make animation.
		animation_fig = animation.FuncAnimation(
			fig, lambda i: update(i, ax_idx, fig), frames=len(data), interval=500, blit=True, repeat_delay=50)
	
	# save figure.
	animation_fig.save(save_path, dpi=100)


def img_bytes2array(fig_bytes):
    buf = io.BytesIO(fig_bytes)
    img = Image.open(buf)
    return np.asarray(img)
