from matplotlib import pyplot as plot
from scipy.spatial import KDTree
import numpy as np 
import pandas as pd 
import os


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
	return covs_l, mean_dists_l


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
	return (covs_l, mean_dists_l)


def slide_window_using_tree(mat, ref_mat, n_neighbors=10):

	assert mat.shape==ref_mat.shape, "matrix and ref. matrix should have the same shape."
	
	# construct a KD tree of cell positions.
	cell_posns = np.array(tuple(zip(*np.nonzero(ref_mat))))
	tree = KDTree(cell_posns)
	# include the cell itself at index 0.
	dists, indices = tree.query(cell_posns, k=n_neighbors+1) 
	# collect stats for all neighborhoods.
	covs_l = []
	mean_dists_l = []
	for i in range(len(cell_posns)):
		# exclude self.
		neighbor_signals = mat[cell_posns[indices[i][1:]]] 
		# compute stats.
		mu = np.mean(neighbor_signals)
		sigma = np.std(neighbor_signals)
		cov = sigma / mu if mu > 0 else 0
		avg_dist = np.mean(dists[i][1:])
		# accumulate stats.
		covs_l.append(cov)
		mean_dists_l.append(avg_dist)

	return (covs_l, mean_dists_l)


# file names of `clouds` file is extrapolated based on these file paths as well.
levels_l = [
	"/home/kd766/quorum-sensing/outputs/06262025215938_size-100x100_select-0.3_seed-0.025/06262025215938_size-100x100_select-0.3_seed-0.025_levels_all.npy",
	"/home/kd766/quorum-sensing/outputs/06262025220039_size-100x100_select-0.3_seed-0.0333/06262025220039_size-100x100_select-0.3_seed-0.0333_levels_all.npy",
	"/home/kd766/quorum-sensing/outputs/06262025220139_size-100x100_select-0.3_seed-0.0667/06262025220139_size-100x100_select-0.3_seed-0.0667_levels_all.npy",
	"/home/kd766/quorum-sensing/outputs/05292025180028_size-100x100_select-0.3_seed-0.1/05292025180028_size-100x100_select-0.3_seed-0.1_levels_all.npy",
	"/home/kd766/quorum-sensing/outputs/05292025180016_size-100x100_select-0.3_seed-0.125/05292025180016_size-100x100_select-0.3_seed-0.125_levels_all.npy",
	"/home/kd766/quorum-sensing/outputs/05292025174657_size-100x100_select-0.4_seed-0.15/05292025174657_size-100x100_select-0.4_seed-0.15_levels_all.npy",
	"/home/kd766/quorum-sensing/outputs/05292025180030_size-100x100_select-0.3_seed-0.175/05292025180030_size-100x100_select-0.3_seed-0.175_levels_all.npy",
	"/home/kd766/quorum-sensing/outputs/05292025180015_size-100x100_select-0.3_seed-0.2/05292025180015_size-100x100_select-0.3_seed-0.2_levels_all.npy",
	# "/home/kd766/quorum-sensing/outputs/06262025142707_size-100x100_select-0.3_seed-0.225/06262025142707_size-100x100_select-0.3_seed-0.225_levels_all.npy",
	# "/home/kd766/quorum-sensing/outputs/06262025142751_size-100x100_select-0.3_seed-0.25/06262025142751_size-100x100_select-0.3_seed-0.25_levels_all.npy",
	# "/home/kd766/quorum-sensing/outputs/06262025142714_size-100x100_select-0.3_seed-0.275/06262025142714_size-100x100_select-0.3_seed-0.275_levels_all.npy",
]

# timepoints at which to analyze each of the above experiments; one set per setup.
timepoints_mat = [
	[x for x in range(1, 49, 5)],  	# 0.025
	[x for x in range(1, 49, 5)],		# 0.033
	[x for x in range(1, 49, 4)],		# 0.0667
	[x for x in range(1, 49, 4)],		# 0.0667
	[x for x in range(1, 49, 4)],		# 0.0667
	[x for x in range(1, 49, 4)],		# 0.0667
	[x for x in range(1, 49, 4)],		# 0.0667
	[x for x in range(1, 49, 4)],		# 0.0667
]

# set the reqd. cells per window range -- window size will be scaled based on density.
for levels_fpath, timepoints_l in zip(levels_l, timepoints_mat):
	
	# infer simulation config.
	dirpath, fname = os.path.split(levels_fpath)
	clouds_fpath = os.path.join(dirpath, '_'.join(fname.split('_')[:-2]) + '_clouds_all.npy')
	seeding_density_str = fname.split('_')[-3].split('-')[-1]
	seeding_density = round(float(seeding_density_str), 4)
	
	# read files.
	levels = np.load(levels_fpath)
	clouds = np.load(clouds_fpath)

	# compute window/neighborhood size range.
	neighborhood_range = [x for x in range(6, 16, 1)]
	print(neighborhood_range)
	
	# run sliding window on levels and clouds.
	levels_fig = plot.figure(figsize=(16, 12))
	clouds_fig = plot.figure(figsize=(16, 12))
	cloud_overlay_idx = 0
	for idx, r in enumerate(neighborhood_range):
		
		# levels.
		plot.figure(levels_fig)
		ax = plot.subplot(2, len(neighborhood_range)//2, idx+1)
		for tpoint in timepoints_l:
			covs_l, mean_dists_l = slide_window_using_tree(
				mat=levels[tpoint, :, :], ref_mat=levels[tpoint, :, :], n_neighbors=r)
			ax.scatter(mean_dists_l, covs_l, s=0.8, label=f"t={tpoint}", alpha=0.5)
		ax.set_xlabel("mean cell–cell distance [plate units]")
		ax.set_ylabel("CoV (signal level)")
		ax.set_title(f"# neighbors = {r}")
		
		# clouds.
		plot.figure(clouds_fig)
		ax = plot.subplot(2, len(neighborhood_range)//2, idx+1)
		for tpoint in timepoints_l:
			covs_l, mean_dists_l = slide_window_using_tree(
				mat=clouds[tpoint, :, :], ref_mat=levels[tpoint, :, :], n_neighbors=r)
			ax.scatter(mean_dists_l, covs_l, s=0.8, label=f"t={tpoint}", alpha=0.5)
		ax.set_xlabel("mean cell–cell distance [plate units]")
		ax.set_ylabel("CoV (cloud intensity)")
		ax.set_title(f"# neighbors = {r}")
	
	# save the figures.
	plot.figure(levels_fig)
	plot.legend()
	plot.suptitle(f"signaling intensity at density={seeding_density*100}%")
	plot.savefig(os.path.join(
		'analysis_outputs/temporal', f'temporal_local-cov-dist_levels_select-{str(seeding_density).ljust(4, '0')}.png'), dpi=100)
	plot.figure(clouds_fig)
	plot.legend()
	plot.suptitle(f"cloud intensity at density={seeding_density*100}%")
	plot.savefig(os.path.join(
		'analysis_outputs/temporal', f'temporal_local-cov-dist_cloud_select-{str(seeding_density).ljust(4, '0')}.png'), dpi=100)