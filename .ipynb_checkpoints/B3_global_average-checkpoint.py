import glob 
import os
import math
import pandas as pd 
import numpy as np
from matplotlib import pyplot as plot


# empty string is for the case with negative feedback.
for feedback_str in ['', 'noneg']:

	if feedback_str == '':
		filter_func = lambda x: 'noneg' not in x 
	else:
		filter_func = lambda x: 'noneg' in x

	levels_l = np.array(list(filter(filter_func, glob.glob(os.path.join(
		"/home/kd766/quorum-sensing/outputs/subexp",
		f"*/*levels_final.csv")))))
	clouds_l = np.array(list(filter(filter_func, glob.glob(os.path.join(
		"/home/kd766/quorum-sensing/outputs/subexp",
		f"*/*clouds_final.csv")))))

	n_trials = len(levels_l)
	print(f"Total # trials: {n_trials}")
	n_hists = 2*24
	fig_grid_size = (math.ceil(n_hists**0.5), math.ceil(n_hists**0.5))


	density_str_l = []
	# infer all densities and pre-order them into a list.
	for idx, (level_file) in enumerate(levels_l):

		fname = os.path.basename(level_file)
		parts = fname.split('_')
		density = float(parts[3].split('-')[-1])
		density_str = str(round(density, 4))
		density_str_l.append(density_str)

	files_sorted_order_idx = np.argsort(density_str_l)


	hist_counter_d = dict()
	hist_ctr = 0
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

		hist_counter_d[density_str] = hist_counter_d.get(density_str, 0) + 1
		# plot only first two histograms.
		if hist_counter_d[density_str] <= 2:
			hist_ctr += 1
			# plot levels.
			plot.figure(levels_fig)
			ax = plot.subplot(*fig_grid_size, hist_ctr)
			ax.hist(levels.flatten(), bins=8, density=True)
			ax.set_title(f"d = {round(density, 4)}")
			ax.set_yscale('log')
			ax.set_ylim(1e-5, 1)
			# plot clouds.
			plot.figure(clouds_fig)
			ax = plot.subplot(*fig_grid_size, hist_ctr)
			ax.hist(cloud.flatten(), bins=12, density=True)
			ax.set_title(f"d = {round(density, 4)}")
			ax.set_yscale('log')
			ax.set_ylim(1e-5, 1)

	plot.figure(levels_fig)
	plot.tight_layout()
	plot.savefig(os.path.join(
		f'analysis_outputs/kdtree/hist/{feedback_str}', f'kdtree_{feedback_str}_global-mean_levels_hist.png'), dpi=100)
	plot.clf()

	plot.figure(clouds_fig)
	plot.tight_layout()
	plot.savefig(os.path.join(
		f'analysis_outputs/kdtree/hist/{feedback_str}', f'kdtree_{feedback_str}_global-mean_clouds_hist.png'), dpi=100)
	plot.clf()


	# plot levels trend.
	plot.figure(figsize=(10, 10))
	for density_str in level_means_d.keys():
		trials_mean = np.mean(level_means_d[density_str])
		trials_std = np.std(level_means_d[density_str])
		# plot centrality measures.
		plot.errorbar(density_d[density_str], trials_mean, yerr=trials_std, ecolor='red', color='blue')
		# plot all pts.
		mean_signal_mean = np.mean(level_means_d[density_str])
		plot.scatter(density_d[density_str], mean_signal_mean, s=2.5, c='blue', marker='*')
		for signal_mean in level_means_d[density_str]:
			plot.scatter(density_d[density_str], signal_mean, s=2, c='green', marker='x')

	plot.xlabel("cell seeding density")
	plot.ylabel("global mean levels")
	plot.savefig(os.path.join(
		f'analysis_outputs/kdtree/{feedback_str}', f'kdtree_{feedback_str}_global-density-mean_levels.png'), dpi=100)
	plot.clf()


	# plot clouds trend.
	plot.figure(figsize=(10, 10))
	for density_str in cloud_means_d.keys():
		trials_mean = np.mean(cloud_means_d[density_str])
		trials_std = np.std(cloud_means_d[density_str])
		# plot centrality measures.
		plot.errorbar(density_d[density_str], trials_mean, yerr=trials_std, ecolor='red', color='blue')
		# plot all pts.
		mean_signal_mean = np.mean(cloud_means_d[density_str])
		plot.scatter(density_d[density_str], mean_signal_mean, s=2.5, c='blue', marker='*')
		for signal_mean in cloud_means_d[density_str]:
			plot.scatter(density_d[density_str], signal_mean, s=2, c='green', marker='x')

	plot.xlabel("cell seeding density")
	plot.ylabel("global mean cloud")
	plot.savefig(os.path.join(
		f'analysis_outputs/kdtree/{feedback_str}', f'kdtree_{feedback_str}_global-density-mean_clouds.png'), dpi=100)
	plot.clf()
