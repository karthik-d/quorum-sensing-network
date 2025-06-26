from matplotlib import pyplot as plot
from scipy.spatial import KDTree
import numpy as np 
import pandas as pd 
import os


def compute_density():
	return None 


def _old_slide_window(mat, r=10):
	"""
	previous method, sliding a constant window over the entire array
	"""
	windows = np.lib.stride_tricks.sliding_window_view(mat, window_shape=(r, r))
	means_l = []
	stdevs_l = []
	n_cells_l = [] 
	# print(mat.shape, windows.shape)
	for win_row in windows:
		for win in win_row:
			vals = win[np.nonzero(win)].flatten()
			n_cells_l.append(len(vals))
			means_l.append(np.mean(vals))
			stdevs_l.append(np.std(vals))
	return means_l, stdevs_l, n_cells_l


def slide_window(mat, ref_mat, neighbor_thresh=10):
	"""
	go cell-to-cell and find neighbors, instead of sliding the window everywhere
	- neighbor_thresh: how far another cell can be to be considered a neighbor.
	"""
	assert mat.shape==ref_mat.shape, "matrix and ref. matrix should have the same shape."
	
	# stats to collect.
	means_l = []
	stdevs_l = []
	n_cells_l = [] 
	# locate positions of cells.
	cell_posns = np.nonzero(ref_mat)
	# iteration helpers.
	r = neighbor_thresh		# shorthand.
	x_max, y_max = [dim-1 for dim in mat.shape]
	windows_l = []
	for posn_x, posn_y in zip(*cell_posns):
		# get window -- region of interest.
		top, bottom = max(0, posn_x-r), min(x_max, posn_x+r)	# rows.
		left, right = max(0, posn_y-r), min(y_max, posn_y+r)	# cols.
		mat_window = mat[top:bottom+1, left:right+1]
		ref_mat_window = ref_mat[top:bottom+1, left:right+1]
		# find locations of cells, and pick the values at those positions.
		vals = mat_window[np.nonzero(ref_mat_window)].flatten()
		# store stats.
		n_cells_l.append(len(vals))
		means_l.append(np.mean(vals))
		stdevs_l.append(np.std(vals))
		windows_l.append(((left, top), right-left, bottom-top))

	# return stats.
	return (means_l, stdevs_l, n_cells_l)


def slide_window_using_tree(mat, ref_mat, n_neighbors=10):

	assert mat.shape==ref_mat.shape, "matrix and ref. matrix should have the same shape."

	# construct a KD tree of cell positions.
	cell_posns = np.array(tuple(zip(*np.nonzero(ref_mat))))
	tree = KDTree(cell_posns)
	# include the cell itself at index 0.
	dists, indices = tree.query(cell_posns, k=n_neighbors+1) 
	# collect stats for all neighborhoods.
	means_l = []
	stdevs_l = []
	for i in range(len(cell_posns)):
		# exclude self.
		neighbor_signals = mat[cell_posns[indices[i][1:]]] 
		# compute stats.
		mu = np.mean(neighbor_signals)
		sigma = np.std(neighbor_signals)
		# accumulate stats.
		means_l.append(mu)
		stdevs_l.append(sigma)

	return (means_l, stdevs_l, n_neighbors)


# file names of `clouds` file is extrapolated based on these file paths as well.
levels_l = [
	"/home/kd766/quorum-sensing/outputs/04112025064926_size-100x100_select-0.3_seed-0.025/04112025064926_size-100x100_select-0.3_seed-0.025_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/04112025064724_size-100x100_select-0.3_seed-0.0333/04112025064724_size-100x100_select-0.3_seed-0.0333_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/05292025212917_size-100x100_select-1_seed-0.0667/05292025212917_size-100x100_select-1_seed-0.0667_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/05292025180028_size-100x100_select-0.3_seed-0.1/05292025180028_size-100x100_select-0.3_seed-0.1_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/05292025180016_size-100x100_select-0.3_seed-0.125/05292025180016_size-100x100_select-0.3_seed-0.125_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/05292025174657_size-100x100_select-0.4_seed-0.15/05292025174657_size-100x100_select-0.4_seed-0.15_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/05292025180030_size-100x100_select-0.3_seed-0.175/05292025180030_size-100x100_select-0.3_seed-0.175_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/05292025180015_size-100x100_select-0.3_seed-0.2/05292025180015_size-100x100_select-0.3_seed-0.2_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/06262025142707_size-100x100_select-0.3_seed-0.225/06262025142707_size-100x100_select-0.3_seed-0.225_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/06262025142751_size-100x100_select-0.3_seed-0.25/06262025142751_size-100x100_select-0.3_seed-0.25_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/06262025142714_size-100x100_select-0.3_seed-0.275/06262025142714_size-100x100_select-0.3_seed-0.275_levels_final.csv",
]

# set the reqd. cells per window range -- window size will be scaled based on dennsity.
reqd_cells_per_win_range = list(range(4, 19, 2))
for levels_fpath in levels_l:
	
	# infer simulation config.
	dirpath, fname = os.path.split(levels_fpath)
	clouds_fpath = os.path.join(dirpath, '_'.join(fname.split('_')[:-2]) + '_clouds_final.csv')
	seeding_density = round(float(fname.split('_')[-3].split('-')[-1]), 4)
	
	# read files.
	levels = pd.read_csv(levels_fpath, header=None).to_numpy()
	clouds = pd.read_csv(clouds_fpath, header=None).to_numpy()

	# compute window/neighborhood size range.
	neighborhood_range = [x for x in range(4, 12, 1)]
	print(neighborhood_range)
	
	# run sliding window on levels and clouds.
	levels_fig = plot.figure(figsize=(16, 12))
	clouds_fig = plot.figure(figsize=(16, 12))
	for idx, r in enumerate(neighborhood_range):
		
		# levels.
		means_l, stdevs_l, n_cells_l = slide_window_using_tree(levels, ref_mat=levels, n_neighbors=r)
		plot.figure(levels_fig)
		ax = plot.subplot(2, len(neighborhood_range)//2, idx+1)
		ax.scatter(means_l, stdevs_l, s=0.8)
		ax.set_xlabel("mean signal")
		ax.set_ylabel("std dev. signal")
		ax.set_title(f"# cells / window = {round(np.mean(n_cells_l), 2)}")
		
		# clouds.
		means_l, stdevs_l, n_cells_l = slide_window_using_tree(clouds, ref_mat=levels, n_neighbors=r)
		plot.figure(clouds_fig)
		ax = plot.subplot(2, len(neighborhood_range)//2, idx+1)
		ax.scatter(means_l, stdevs_l, s=0.8)
		ax.set_xlabel("mean cloud")
		ax.set_ylabel("std dev. cloud")
		ax.set_title(f"# cells / window = {round(np.mean(n_cells_l), 2)}")
	
	# save the figures.
	plot.figure(levels_fig)
	plot.suptitle(f"signaling intensity at density={seeding_density*100}%")
	plot.savefig(os.path.join(
		'analysis_outputs', f'local-mean-stdev_levels_select-{str(seeding_density).ljust(4, '0')}.png'), dpi=100)
	plot.figure(clouds_fig)
	plot.suptitle(f"cloud intensity at density={seeding_density*100}%")
	plot.savefig(os.path.join(
		'analysis_outputs', f'local-mean-stdev_cloud_select-{str(seeding_density).ljust(4, '0')}.png'), dpi=100)

exit()
n_cells_l = []
# 10% -- 7, 8.
fig = plot.figure(figsize=(8, 8))
levels = pd.read_csv(levels_l[0], header=None).to_numpy()
for idx, r in enumerate(range(2, 8)):
	mean, stdev, n_cells = slide_window(levels, neighbor_thresh=r)
	plot.figure(fig)
	ax = plot.subplot(2, 3, idx+1)
	ax.scatter(mean, stdev, s=0.8)
	ax.set_xlabel("mean signal")
	ax.set_ylabel("std dev. signal")
	ax.set_title(f"avg. cells / window = {round(np.mean(n_cells), 2)}")
	n_cells_l.append(np.mean(n_cells))
plot.tight_layout()
plot.savefig("frac10_sliding.png")

# TODO: reference cell posns using levels for clouds; they cannot be inferred using the cloud array alone!
# plot.figure(figsize=(8, 8))
# clouds = pd.read_csv(clouds_l[0], header=None).to_numpy()
# for idx, r in enumerate(range(2, 8)):
# 	mean, stdev, n_cells = slide_window(clouds, neighbor_thresh=r)
# 	ax = plot.subplot(2, 3, idx+1)
# 	ax.scatter(mean, stdev, s=0.8)
# 	ax.set_xlabel("mean signal")
# 	ax.set_ylabel("std dev. signal")
# 	ax.set_title(f"avg. cells / window = {round(n_cells_l[0], 2)}")
# plot.tight_layout()
# plot.savefig("frac10_CLOUD.png")

exit()

fig = plot.figure(figsize=(8, 8))
levels = pd.read_csv(levels_l[1], header=None).to_numpy()
for idx, r in enumerate(range(12, 18)):
	mean, stdev, n_cells = slide_window(levels, neighbor_thresh=r)
	plot.figure(fig)
	ax = plot.subplot(2, 3, idx+1)
	ax.scatter(mean, stdev, s=0.8)
	ax.set_xlabel("mean signal")
	ax.set_ylabel("std dev. signal")
	ax.set_title(f"avg. cells / window = {round(np.mean(n_cells), 2)}")
	n_cells_l.append(np.mean(n_cells))
plot.tight_layout()
plot.savefig("frac3_sliding.png")

plot.figure(figsize=(8, 8))
levels = pd.read_csv(levels_l[2], header=None).to_numpy()
for idx, r in enumerate(range(14, 20)):
	mean, stdev, n_cells = slide_window(levels, neighbor_thresh=r)
	ax = plot.subplot(2, 3, idx+1)
	ax.scatter(mean, stdev, s=0.8)
	ax.set_xlabel("mean signal")
	ax.set_ylabel("std dev. signal")
	ax.set_title(f"avg. cells / window = {round(np.mean(n_cells), 2)}")
	n_cells_l.append(np.mean(n_cells))
plot.tight_layout()
plot.savefig("frac2_sliding.png")

# for cloud_file in levels_l:
	# levels = pd.read_csv(level_file, header=None).to_numpy()
	# print(levels)


plot.figure(figsize=(8, 8))
clouds = pd.read_csv(clouds_l[0], header=None).to_numpy()
for idx, r in enumerate(range(6, 12)):
	mean, stdev, n_cells = slide_window(clouds, neighbor_thresh=r)
	ax = plot.subplot(2, 3, idx+1)
	ax.scatter(mean, stdev, s=0.8)
	ax.set_xlabel("mean signal")
	ax.set_ylabel("std dev. signal")
	ax.set_title(f"avg. cells / window = {round(n_cells_l[0], 2)}")
plot.tight_layout()
plot.savefig("frac10_CLOUD.png")

plot.figure(figsize=(8, 8))
clouds = pd.read_csv(clouds_l[1], header=None).to_numpy()
for idx, r in enumerate(range(12, 18)):
	mean, stdev, n_cells = slide_window(clouds, neighbor_thresh=r)
	ax = plot.subplot(2, 3, idx+1)
	ax.scatter(mean, stdev, s=0.8)
	ax.set_xlabel("mean signal")
	ax.set_ylabel("std dev. signal")
	ax.set_title(f"avg. cells / window = {round(n_cells_l[1], 2)}")
plot.tight_layout()
plot.savefig("frac3_CLOUD.png")

plot.figure(figsize=(8, 8))
clouds = pd.read_csv(clouds_l[2], header=None).to_numpy()
for idx, r in enumerate(range(14, 20)):
	mean, stdev, n_cells = slide_window(clouds, neighbor_thresh=r)
	ax = plot.subplot(2, 3, idx+1)
	ax.scatter(mean, stdev, s=0.8)
	ax.set_xlabel("mean signal")
	ax.set_ylabel("std dev. signal")
	ax.set_title(f"avg. cells / window = {round(n_cells_l[2], 2)}")
plot.tight_layout()
plot.savefig("frac2_CLOUD.png")