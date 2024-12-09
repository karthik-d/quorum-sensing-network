from PIL import Image
from matplotlib import pyplot as plot

import matplotlib.animation as animation
import io
import numpy as np



def save_animation(imgs_l, save_path):

	# animation update functio.
	def update(i):
		im.set_array(imgs_l[i])
		ax.set_title(f"time={i}")
		return im, 

	fig = plot.figure()
	ax = plot.gca()
	# set the initial image.
	im = ax.imshow(imgs_l[0], animated=True, vmin=0, vmax=10)
	ax.set_title(f"time={1}")
	fig.colorbar(im)

	# make animation.
	animation_fig = animation.FuncAnimation(fig, update, frames=len(imgs_l), interval=500, blit=True, repeat_delay=50)
	animation_fig.save(save_path, dpi=100)


def img_bytes2array(fig_bytes):
    buf = io.BytesIO(fig_bytes)
    img = Image.open(buf)
    return np.asarray(img)
