from matplotlib import pyplot as plot
import numpy as np 
import pandas as pd 


def compute_density():
	return None 


def slide_window(mat, r=10):
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


clouds_l = [
	"/home/kd766/quorum-sensing/outputs/sliding-window/04112025064544_size-100x100_select-1_seed-0.1_clouds_final.csv",
	"/home/kd766/quorum-sensing/outputs/sliding-window/04112025064724_size-100x100_select-0.3_seed-0.0333_clouds_final.csv",
	"/home/kd766/quorum-sensing/outputs/sliding-window/04112025064926_size-100x100_select-0.3_seed-0.025_clouds_final.csv",
]

levels_l = [
	"/home/kd766/quorum-sensing/outputs/sliding-window/04112025064544_size-100x100_select-1_seed-0.1_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/sliding-window/04112025064724_size-100x100_select-0.3_seed-0.0333_levels_final.csv",
	"/home/kd766/quorum-sensing/outputs/sliding-window/04112025064926_size-100x100_select-0.3_seed-0.025_levels_final.csv",
]

# for level_file in levels_l:
# 	levels = pd.read_csv(level_file, header=None).to_numpy()
# 	print(np.sum(levels!=0)/1e4)
	# print(slide_window(levels))

n_cells_l = []
# 10% -- 7, 8.
plot.figure(figsize=(8, 8))
levels = pd.read_csv(levels_l[0], header=None).to_numpy()
for idx, r in enumerate(range(6, 12)):
	mean, stdev, n_cells = slide_window(levels, r=r)
	ax = plot.subplot(2, 3, idx+1)
	ax.scatter(mean, stdev, s=0.8)
	ax.set_xlabel("mean signal")
	ax.set_ylabel("std dev. signal")
	ax.set_title(f"avg. cells per window = {round(np.mean(n_cells), 2)}")
	n_cells_l.append(np.mean(n_cells))
plot.tight_layout()
plot.savefig("frac10_sliding.png")

plot.figure(figsize=(8, 8))
levels = pd.read_csv(levels_l[1], header=None).to_numpy()
for idx, r in enumerate(range(12, 18)):
	mean, stdev, n_cells = slide_window(levels, r=r)
	ax = plot.subplot(2, 3, idx+1)
	ax.scatter(mean, stdev, s=0.8)
	ax.set_xlabel("mean signal")
	ax.set_ylabel("std dev. signal")
	ax.set_title(f"avg. cells per window = {round(np.mean(n_cells), 2)}")
	n_cells_l.append(np.mean(n_cells))
plot.tight_layout()
plot.savefig("frac3_sliding.png")

plot.figure(figsize=(8, 8))
levels = pd.read_csv(levels_l[2], header=None).to_numpy()
for idx, r in enumerate(range(14, 20)):
	mean, stdev, n_cells = slide_window(levels, r=r)
	ax = plot.subplot(2, 3, idx+1)
	ax.scatter(mean, stdev, s=0.8)
	ax.set_xlabel("mean signal")
	ax.set_ylabel("std dev. signal")
	ax.set_title(f"avg. cells per window = {round(np.mean(n_cells), 2)}")
	n_cells_l.append(np.mean(n_cells))
plot.tight_layout()
plot.savefig("frac2_sliding.png")

# for cloud_file in levels_l:
	# levels = pd.read_csv(level_file, header=None).to_numpy()
	# print(levels)


plot.figure(figsize=(8, 8))
levels = pd.read_csv(clouds_l[0], header=None).to_numpy()
for idx, r in enumerate(range(6, 12)):
	mean, stdev, n_cells = slide_window(levels, r=r)
	ax = plot.subplot(2, 3, idx+1)
	ax.scatter(mean, stdev, s=0.8)
	ax.set_xlabel("mean signal")
	ax.set_ylabel("std dev. signal")
	ax.set_title(f"avg. cells per window = {round(n_cells_l[0], 2)}")
plot.tight_layout()
plot.savefig("frac10_CLOUD.png")

plot.figure(figsize=(8, 8))
levels = pd.read_csv(clouds_l[1], header=None).to_numpy()
for idx, r in enumerate(range(12, 18)):
	mean, stdev, n_cells = slide_window(levels, r=r)
	ax = plot.subplot(2, 3, idx+1)
	ax.scatter(mean, stdev, s=0.8)
	ax.set_xlabel("mean signal")
	ax.set_ylabel("std dev. signal")
	ax.set_title(f"avg. cells per window = {round(n_cells_l[1], 2)}")
plot.tight_layout()
plot.savefig("frac3_CLOUD.png")

plot.figure(figsize=(8, 8))
levels = pd.read_csv(clouds_l[2], header=None).to_numpy()
for idx, r in enumerate(range(14, 20)):
	mean, stdev, n_cells = slide_window(levels, r=r)
	ax = plot.subplot(2, 3, idx+1)
	ax.scatter(mean, stdev, s=0.8)
	ax.set_xlabel("mean signal")
	ax.set_ylabel("std dev. signal")
	ax.set_title(f"avg. cells per window = {round(n_cells_l[2], 2)}")
plot.tight_layout()
plot.savefig("frac2_CLOUD.png")