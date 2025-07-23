import glob 
import os
import pandas as pd 
import numpy as np
from matplotlib import pyplot as plot


levels_l = np.array(glob.glob(os.path.join(
	"/home/kd766/quorum-sensing/outputs/subexp",
	"*/*levels_final.csv*"
)))
clouds_l = np.array(glob.glob(os.path.join(
	"/home/kd766/quorum-sensing/outputs/subexp",
	"*/*clouds_final.csv*"
)))
print(f"Total # trials: {len(levels_l)}")


density_str_l = []
# infer all densities and pre-order them into a list.
for idx, (level_file) in enumerate(levels_l):

	fname = os.path.basename(level_file)
	parts = fname.split('_')
	density = float(parts[3].split('-')[-1])
	density_str = str(round(density, 4))
	density_str_l.append(density_str)

files_sorted_order_idx = np.argsort(density_str_l)


density_d = dict()
level_means_d = dict()
cloud_means_d = dict()
levels_fig = plot.figure(1, figsize=(16, 16))
clouds_fig = plot.figure(2, figsize=(16, 16))
for idx, (cloud_file, level_file) in enumerate(zip(
	clouds_l[files_sorted_order_idx], levels_l[files_sorted_order_idx])):

	fname = os.path.basename(level_file)
	parts = fname.split('_')
	density = float(parts[3].split('-')[-1])
	density_str = str(round(density, 4))

	cloud = pd.read_csv(cloud_file, header=None).to_numpy()
	levels = pd.read_csv(level_file, header=None).to_numpy()
	n_cells = np.sum(levels!=0)

	cloud_means_d[density_str] = cloud_means_d.get(density_str, []) + [np.sum(cloud)/n_cells]
	level_means_d[density_str] = level_means_d.get(density_str, []) + [np.sum(levels)/n_cells]
	density_d[density_str] = density

	plot.figure(levels_fig)
	ax = plot.subplot(9, 10, idx+1)
	ax.hist(levels.flatten(), bins=8)
	ax.set_title(f"density = {round(density, 4)}")
	ax.set_yscale('log')

	plot.figure(clouds_fig)
	ax = plot.subplot(9, 10, idx+1)
	ax.hist(cloud.flatten(), bins=12)
	ax.set_title(f"density = {round(density, 4)}")
	ax.set_yscale('log')

plot.figure(levels_fig)
plot.tight_layout()
plot.savefig(os.path.join(
	'analysis_outputs/kdtree', f'kdtree_global-mean_levels_hist.png'), dpi=100)

plot.figure(clouds_fig)
plot.tight_layout()
plot.savefig(os.path.join(
	'analysis_outputs/kdtree', f'kdtree_global-mean_clouds_hist.png'), dpi=100)


# plot levels trend.
plot.figure(figsize=(10, 10))
for k in level_means_d.keys():
	trials_mean = np.mean(level_means_d[k])
	trials_std = np.std(level_means_d[k])
	# plot centrality measures.
	plot.errorbar(density_d[k], trials_mean, yerr=trials_std, ecolor='red', color='blue')
	# plot all pts.
	for signal_mean in level_means_d[k]:
		plot.scatter(density_d[k], signal_mean, s=1.4, c='green', marker='x')

plot.xlabel("cell seeding density")
plot.ylabel("global mean levels")
plot.savefig(os.path.join(
	'analysis_outputs/kdtree', f'kdtree_global-density-mean_levels.png'), dpi=100)


# plot clouds trend.
plot.figure(figsize=(10, 10))
for k in cloud_means_d.keys():
	trials_mean = np.mean(cloud_means_d[k])
	trials_std = np.std(cloud_means_d[k])
	# plot centrality measures.
	plot.errorbar(density_d[k], trials_mean, yerr=trials_std, ecolor='red', color='blue')
	# plot all pts.
	for signal_mean in cloud_means_d[k]:
		plot.scatter(density_d[k], signal_mean, s=1.4, c='green', marker='x')

plot.xlabel("cell seeding density")
plot.ylabel("global mean cloud")
plot.savefig(os.path.join(
	'analysis_outputs/kdtree', f'kdtree_global-density-mean_clouds.png'), dpi=100)
