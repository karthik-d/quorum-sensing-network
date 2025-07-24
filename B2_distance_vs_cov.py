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
	# find the enclosing areas on the matrix.
	window_minmax_coords = [(np.min(np.squeeze(arr), axis=1), np.max(np.squeeze(arr), axis=1)) 
		for arr in np.split(np.array(cell_posns[indices]), 2, axis=-1)]
	# collect stats for all neighborhoods.
	covs_l = []
	mean_dists_l = []
	for idx, (min_x, max_x, min_y, max_y) in enumerate(
		zip(*(window_minmax_coords[0] + window_minmax_coords[1]))):
		
		mat_window = mat[min_x:max_x+1, min_y:max_y+1]
		ref_mat_window = ref_mat[min_x:max_x+1, min_y:max_y+1]
		# find locations of cells, and pick the values at those positions.
		# vals = mat_window[np.nonzero(ref_mat_window)].flatten()
		vals = mat_window.flatten()				# use all posns; not just where cells are.
		mu = np.mean(vals, where=(vals!=0))		# compute using non-zero values only.
		sigma = np.std(vals, where=(vals!=0))	# compute using non-zero values only.
		cov = sigma / mu if mu > 0 else 0
		avg_dist = np.mean(dists[idx][1:])
		# store stats.
		covs_l.append(cov)
		mean_dists_l.append(avg_dist)
	
	## older method that simply uses the cell posns.
	# for i in range(len(cell_posns)):
	# 	# exclude self.
	# 	neighbor_signals = mat[cell_posns[indices[i][1:]]] 
	# 	# compute stats.
	# 	mu = np.mean(neighbor_signals)
	# 	sigma = np.std(neighbor_signals)
	# 	# accumulate stats.
	# 	means_l.append(mu)
	# 	stdevs_l.append(sigma)

	return (covs_l, mean_dists_l)


# file names of `clouds` file is extrapolated based on these file paths as well.
levels_withneg_l = [
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

levels_noneg_l = [
	"/home/kd766/quorum-sensing/outputs/07232025152722_size-100x100_select-0.3_seed-0.025_noneg/07232025152722_size-100x100_select-0.3_seed-0.025_noneg_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/07232025162644_size-100x100_select-0.3_seed-0.03_noneg/07232025162644_size-100x100_select-0.3_seed-0.03_noneg_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/07232025152745_size-100x100_select-0.3_seed-0.0333_noneg/07232025152745_size-100x100_select-0.3_seed-0.0333_noneg_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/07232025162614_size-100x100_select-0.3_seed-0.05_noneg/07232025162614_size-100x100_select-0.3_seed-0.05_noneg_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/07232025152552_size-100x100_select-0.3_seed-0.0667_noneg/07232025152552_size-100x100_select-0.3_seed-0.0667_noneg_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/07232025152554_size-100x100_select-0.3_seed-0.1_noneg/07232025152554_size-100x100_select-0.3_seed-0.1_noneg_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/07232025152556_size-100x100_select-0.3_seed-0.125_noneg/07232025152556_size-100x100_select-0.3_seed-0.125_noneg_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/07232025152621_size-100x100_select-0.3_seed-0.15_noneg/07232025152621_size-100x100_select-0.3_seed-0.15_noneg_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/07232025152648_size-100x100_select-0.3_seed-0.175_noneg/07232025152648_size-100x100_select-0.3_seed-0.175_noneg_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/07232025152708_size-100x100_select-0.3_seed-0.2_noneg/07232025152708_size-100x100_select-0.3_seed-0.2_noneg_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/07232025164017_size-100x100_select-0.3_seed-0.225_noneg/07232025164017_size-100x100_select-0.3_seed-0.225_noneg_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/07232025164023_size-100x100_select-0.3_seed-0.25_noneg/07232025164023_size-100x100_select-0.3_seed-0.25_noneg_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/07232025164056_size-100x100_select-0.3_seed-0.275_noneg/07232025164056_size-100x100_select-0.3_seed-0.275_noneg_levels_final.csv",
]


# select some density to overlay on a single plot.
# strings must match how they occur in the filename.
cloud_overlay_densities_l = ["0.025", "0.03", "0.0333", "0.05", "0.0667", "0.1", "0.125", "0.15", 
	"0.175", "0.2", "0.225", "0.25", "0.275"]
levels_overlay_densities_l = ["0.025", "0.03", "0.0333", "0.05", "0.0667", "0.1", "0.125", "0.15", 
	"0.175", "0.2", "0.225", "0.25", "0.275"]

cloud_density_vals = [float(x) for x in cloud_overlay_densities_l]
levels_density_vals = [float(x) for x in levels_overlay_densities_l]
print(cloud_density_vals)


for feedback_str in ["", "noneg"]:
	# set the reqd. cells per window range -- window size will be scaled based on dennsity.
	cloud_overlay_fig = plot.figure(figsize=(16, 12))
	levels_overlay_fig = plot.figure(figsize=(16, 12))
	reqd_cells_per_win_range = list(range(4, 19, 2))
	levels_l = levels_noneg_l if feedback_str=="noneg" else levels_withneg_l
	
	cloud_overlay_idx = 0
	levels_overlay_idx = 0
	for levels_fpath in levels_l:
		
		# infer simulation config.
		dirpath, fname = os.path.split(levels_fpath)
		clouds_fpath = os.path.join(dirpath, '_'.join(fname.split('_')[:-2]) + '_clouds_final.csv')
		seeding_density_str = fname.split('_')[3].split('-')[-1]
		seeding_density = round(float(seeding_density_str), 4)
		
		# read files.
		levels = pd.read_csv(levels_fpath, header=None).to_numpy()
		clouds = pd.read_csv(clouds_fpath, header=None).to_numpy()

		# compute window/neighborhood size range.
		neighborhood_range = [x for x in range(4, 37, 4)]
		
		# run sliding window on levels and clouds.
		levels_fig = plot.figure(figsize=(16, 12))
		clouds_fig = plot.figure(figsize=(16, 12))
		levels_overlayed = False
		cloud_overlayed = False 
		for idx, r in enumerate(neighborhood_range):
			
			# levels.
			# covs_l, mean_dists_l = slide_window(levels, ref_mat=levels, neighbor_thresh=r)
			covs_l, mean_dists_l = slide_window_using_tree(levels, ref_mat=levels, n_neighbors=r)
			plot.figure(levels_fig)
			ax = plot.subplot(2, len(neighborhood_range)//2+1, idx+1)
			ax.scatter(mean_dists_l, covs_l, s=0.8)
			ax.set_xlabel("mean cell–cell distance (plate units)")
			ax.set_ylabel("CoV (signal level)")
			ax.set_title(f"# neighbors = {r}")

			# [levels] add plot to overlay figure, if in list.
			if seeding_density_str in levels_overlay_densities_l:
				plot.figure(levels_overlay_fig)
				ax = plot.subplot(2, len(neighborhood_range)//2+1, idx+1)
				levels_overlay = ax.scatter(mean_dists_l, covs_l, s=0.9, cmap='viridis', 
					c=[levels_density_vals[levels_overlay_idx]]*len(covs_l),
					vmin=min(levels_density_vals), vmax=max(levels_density_vals),
					label=seeding_density_str, alpha=0.5)
				ax.set_xlabel("mean cell–cell distance (plate units)")
				ax.set_ylabel("CoV (signal level)")
				ax.set_title(f"# neighbors = {r}")
				levels_overlayed = True
			
			# clouds.
			# covs_l, mean_dists_l = slide_window(clouds, ref_mat=levels, neighbor_thresh=r)
			covs_l, mean_dists_l = slide_window_using_tree(clouds, ref_mat=levels, n_neighbors=r)
			plot.figure(clouds_fig)
			ax = plot.subplot(2, len(neighborhood_range)//2+1, idx+1)
			ax.scatter(mean_dists_l, covs_l, s=0.8)
			ax.set_xlabel("mean cell–cell distance (plate units)")
			ax.set_ylabel("CoV (cloud intensity)")
			ax.set_title(f"# neighbors = {r}")

			# [clouds] add plot to overlay figure, if in list.
			if seeding_density_str in cloud_overlay_densities_l:
				plot.figure(cloud_overlay_fig)
				ax = plot.subplot(2, len(neighborhood_range)//2+1, idx+1)
				cloud_overlay = ax.scatter(mean_dists_l, covs_l, s=0.9, cmap='viridis', 
					c=[cloud_density_vals[cloud_overlay_idx]]*len(covs_l), 
					vmin=min(cloud_density_vals), vmax=max(cloud_density_vals),
					label=seeding_density_str, alpha=0.5)
				ax.set_xlabel("mean cell–cell distance (plate units)")
				ax.set_ylabel("CoV (cloud intensity)")
				ax.set_title(f"# neighbors = {r}")
				cloud_overlayed = True
		
		# save the figures.
		plot.figure(levels_fig)
		plot.suptitle(f"signaling intensity at density={seeding_density*100}%")
		plot.savefig(os.path.join(
			f'analysis_outputs/kdtree/{feedback_str}', 
			f'kdtree_{feedback_str}_local-cov-dist_levels_select-{str(seeding_density).ljust(4, '0')}.png'), dpi=100)
		plot.close()

		plot.figure(clouds_fig)
		plot.suptitle(f"cloud intensity at density={seeding_density*100}%")
		plot.savefig(os.path.join(
			f'analysis_outputs/kdtree/{feedback_str}', 
			f'kdtree_{feedback_str}_local-cov-dist_cloud_select-{str(seeding_density).ljust(4, '0')}.png'), dpi=100)
		plot.close()

		# update scatter color index.
		levels_overlay_idx += 1 if levels_overlayed else 0
		cloud_overlay_idx += 1 if cloud_overlayed else 0
		

	# save levels overlay figure.
	plot.figure(levels_overlay_fig)
	# plot.legend()
	plot.colorbar(levels_overlay)
	plot.suptitle("signal levels overlayed for select densities")
	plot.savefig(os.path.join(
		f'analysis_outputs/kdtree/{feedback_str}', 
		f'kdtree_{feedback_str}_local-cov-dist_levels_overlay.png'), dpi=100)
	plot.close()

	# save cloud overlay figure.
	plot.figure(cloud_overlay_fig)
	# plot.legend()
	plot.colorbar(cloud_overlay)
	plot.suptitle("cloud intensities overlayed for select densities")
	plot.savefig(os.path.join(
		f'analysis_outputs/kdtree/{feedback_str}', 
		f'kdtree_{feedback_str}_local-cov-dist_clouds_overlay.png'), dpi=100)
	plot.close()