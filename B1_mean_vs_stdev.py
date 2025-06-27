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
	"/home/kd766/quorum-sensing/outputs/06262025214226_size-100x100_select-0.8_seed-0.03/06262025214226_size-100x100_select-0.8_seed-0.03_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/04112025064724_size-100x100_select-0.3_seed-0.0333/04112025064724_size-100x100_select-0.3_seed-0.0333_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/06262025214211_size-100x100_select-0.8_seed-0.05/06262025214211_size-100x100_select-0.8_seed-0.05_levels_final.csv",
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

# select some density to overlay on a single plot.
# strings must match how they occur in the filename.
cloud_overlay_densities_l = ["0.025", "0.03", "0.0333", "0.05", "0.0667", "0.1", "0.125", "0.15", "0.175"]
levels_overlay_densities_l = ["0.025", "0.03", "0.0333", "0.05", "0.0667", "0.1", "0.125", "0.15", 
	"0.175", "0.2", "0.225", "0.25", "0.275"]


# set the reqd. cells per window range -- window size will be scaled based on density.
cloud_overlay_fig = plot.figure(figsize=(16, 12))
levels_overlay_fig = plot.figure(figsize=(16, 12))
for levels_fpath in levels_l:
	
	# infer simulation config.
	dirpath, fname = os.path.split(levels_fpath)
	clouds_fpath = os.path.join(dirpath, '_'.join(fname.split('_')[:-2]) + '_clouds_final.csv')
	seeding_density_str = fname.split('_')[-3].split('-')[-1]
	seeding_density = round(float(seeding_density_str), 4)
	
	# read files.
	levels = pd.read_csv(levels_fpath, header=None).to_numpy()
	clouds = pd.read_csv(clouds_fpath, header=None).to_numpy()

	# compute window/neighborhood size range.
	neighborhood_range = [x for x in range(4, 12, 1)]
	print(neighborhood_range)
	
	# run sliding window on levels and clouds.
	levels_fig = plot.figure(figsize=(16, 12))
	clouds_fig = plot.figure(figsize=(16, 12))
	cloud_overlay_idx = 0
	for idx, r in enumerate(neighborhood_range):
		
		# levels.
		means_l, stdevs_l, n_cells_l = slide_window_using_tree(levels, ref_mat=levels, n_neighbors=r)
		plot.figure(levels_fig)
		ax = plot.subplot(2, len(neighborhood_range)//2, idx+1)
		ax.scatter(means_l, stdevs_l, s=0.8)
		ax.set_xlabel("mean signal")
		ax.set_ylabel("std dev. signal")
		ax.set_title(f"# cells / window = {round(np.mean(n_cells_l), 2)}")

		# [levels] add plot to overlay figure, if in list.
		if seeding_density_str in levels_overlay_densities_l:
			plot.figure(levels_overlay_fig)
			ax = plot.subplot(2, len(neighborhood_range)//2, idx+1)
			ax.scatter(means_l, stdevs_l, s=0.8, label=seeding_density_str,
				alpha=0.5)
			ax.set_xlabel("mean level")
			ax.set_ylabel("std dev. level")
			ax.set_title(f"# cells / window = {round(np.mean(n_cells_l), 2)}")
		
		# clouds.
		means_l, stdevs_l, n_cells_l = slide_window_using_tree(clouds, ref_mat=levels, n_neighbors=r)
		plot.figure(clouds_fig)
		ax = plot.subplot(2, len(neighborhood_range)//2, idx+1)
		ax.scatter(means_l, stdevs_l, s=0.8)
		ax.set_xlabel("mean cloud")
		ax.set_ylabel("std dev. cloud")
		ax.set_title(f"# cells / window = {round(np.mean(n_cells_l), 2)}")

		# [clouds] add plot to overlay figure, if in list.
		if seeding_density_str in cloud_overlay_densities_l:
			plot.figure(cloud_overlay_fig)
			ax = plot.subplot(2, len(neighborhood_range)//2, idx+1)
			ax.scatter(means_l, stdevs_l, s=0.8, label=seeding_density_str,
				alpha=0.5)
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

# save levels overlay figure.
plot.figure(levels_overlay_fig)
plot.legend()
plot.suptitle("signal levels overlayed for select densities")
plot.savefig(os.path.join(
	'analysis_outputs', f'local-mean-stdev_levels_overlay.png'), dpi=100)

# save cloud overlay figure.
plot.figure(cloud_overlay_fig)
plot.legend()
plot.suptitle("cloud intensities overlayed for select densities")
plot.savefig(os.path.join(
	'analysis_outputs', f'local-mean-stdev_clouds_overlay.png'), dpi=100)