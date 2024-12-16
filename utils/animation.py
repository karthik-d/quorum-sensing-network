from PIL import Image
from matplotlib import pyplot as plot

import matplotlib.animation as animation
import io
import numpy as np



def save_animation(data, save_path):

	# each sublist should be a list of images; so data should be a 4D array [subplot_series, timepoint, img_width, img_height].
	num_subplots = len(data)

	# animation update functio.
	def update(i, axes, imgs):
		for ax_idx in range(num_subplots):
			imgs[ax_idx].set_array(data[ax_idx][i])
			axes[ax_idx].set_title(f"time={i}")
		return imgs[ax_idx], 

	fig, axes = plot.subplots(1, num_subplots)
	axes = [axes] if not isinstance(axes, (list, np.ndarray)) else axes
	imgs = [None for _ in range(num_subplots)]
	for ax_idx in range(num_subplots):
		# set the initial image.
		imgs[ax_idx] = axes[ax_idx].imshow(data[ax_idx][0], animated=True, vmin=0, vmax=10, aspect=1)
		axes[ax_idx].set_title(f"time={1}")
		fig.colorbar(imgs[ax_idx])
	
	# make animation.
	animation_fig = animation.FuncAnimation(
		fig, lambda i: update(i, axes, imgs), frames=len(data[0]), interval=500, blit=True, repeat_delay=50)
	# save figure.
	animation_fig.save(save_path, dpi=100)


def img_bytes2array(fig_bytes):
    buf = io.BytesIO(fig_bytes)
    img = Image.open(buf)
    return np.asarray(img)
