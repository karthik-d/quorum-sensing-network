from matplotlib import pyplot as plot
from scipy.spatial import KDTree
import numpy as np 
import pandas as pd 
import os


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
	"/home/kd766/quorum-sensing/outputs/04112025064926_size-100x100_select-0.3_seed-0.025/04112025064926_size-100x100_select-0.3_seed-0.025_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/04112025064724_size-100x100_select-0.3_seed-0.0333/04112025064724_size-100x100_select-0.3_seed-0.0333_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/05292025212917_size-100x100_select-1_seed-0.0667/05292025212917_size-100x100_select-1_seed-0.0667_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/05292025180028_size-100x100_select-0.3_seed-0.1/05292025180028_size-100x100_select-0.3_seed-0.1_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/05292025180016_size-100x100_select-0.3_seed-0.125/05292025180016_size-100x100_select-0.3_seed-0.125_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/05292025174657_size-100x100_select-0.4_seed-0.15/05292025174657_size-100x100_select-0.4_seed-0.15_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/05292025180030_size-100x100_select-0.3_seed-0.175/05292025180030_size-100x100_select-0.3_seed-0.175_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/05292025180015_size-100x100_select-0.3_seed-0.2/05292025180015_size-100x100_select-0.3_seed-0.2_levels_final.csv",
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
	
	# run sliding window on levels and clouds.
	levels_fig = plot.figure(figsize=(16, 12))
	clouds_fig = plot.figure(figsize=(16, 12))
	for idx, r in enumerate(neighborhood_range):
		
		# levels.
		covs_l, mean_dists_l = slide_window_using_tree(levels, ref_mat=levels, n_neighbors=r)
		plot.figure(levels_fig)
		ax = plot.subplot(2, len(neighborhood_range)//2, idx+1)
		ax.scatter(mean_dists_l, covs_l, s=0.8)
		ax.set_xlabel("mean cell–cell distance (plate units)")
		ax.set_ylabel("CoV(signal)")
		ax.set_title(f"# neighbors = {r}")
		
		# clouds.
		covs_l, mean_dists_l = slide_window_using_tree(clouds, ref_mat=levels, n_neighbors=r)
		plot.figure(clouds_fig)
		ax = plot.subplot(2, len(neighborhood_range)//2, idx+1)
		ax.scatter(mean_dists_l, covs_l, s=0.8)
		ax.set_xlabel("mean cell–cell distance (plate units)")
		ax.set_ylabel("CoV(cloud)")
		ax.set_title(f"# neighbors = {r}")
	
	# save the figures.
	plot.figure(levels_fig)
	plot.suptitle(f"signaling intensity at density={seeding_density*100}%")
	plot.savefig(os.path.join(
		'analysis_outputs', f'local-cov-dist_levels_select-{str(seeding_density).ljust(4, '0')}.png'), dpi=100)
	plot.figure(clouds_fig)
	plot.suptitle(f"cloud intensity at density={seeding_density*100}%")
	plot.savefig(os.path.join(
		'analysis_outputs', f'local-cov-dist_cloud_select-{str(seeding_density).ljust(4, '0')}.png'), dpi=100)
	
	# break