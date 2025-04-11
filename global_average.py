import glob 
import os
import pandas as pd 
import numpy as np
from matplotlib import pyplot as plot


levels_l = glob.glob(os.path.join(
	"/home/kd766/quorum-sensing/outputs/subexp",
	"*levels_final.csv*"
))

clouds_l = glob.glob(os.path.join(
	"/home/kd766/quorum-sensing/outputs/subexp",
	"*clouds_final.csv*"
))

print(len(levels_l))
n_cells_l = []
density_l = []
means_l = []
plot.figure(1, figsize=(10, 10))
for idx, level_file in enumerate(levels_l):

	fname = os.path.basename(level_file)
	parts = fname.split('_')
	density = float(parts[3].split('-')[-1])
	density_l.append(density)

	levels = pd.read_csv(level_file, header=None).to_numpy()
	n_cells = np.sum(levels!=0)
	n_cells_l.append(n_cells)

	means_l.append(np.sum(levels)/n_cells)

	ax = plot.subplot(5, 6, idx+1)
	ax.hist(levels.flatten(), bins=8)
	ax.set_title(f"density = {round(density, 4)}")
	ax.set_yscale('log')

plot.tight_layout()
plot.savefig("check_hist.png")


pts_order = np.argsort(density_l)
density_l = np.array(density_l)
means_l = np.array(means_l)
plot.figure(2)
plot.scatter(density_l[pts_order], means_l[pts_order])
plot.xlabel("cell seeding density")
plot.ylabel("global mean signal output")
plot.savefig("check.png")



# --
n_cells_l = []
density_l = []
means_l = []
plot.figure(3, figsize=(10, 10))
for idx, level_file in enumerate(clouds_l):

	fname = os.path.basename(level_file)
	parts = fname.split('_')
	density = float(parts[3].split('-')[-1])
	density_l.append(density)

	levels = pd.read_csv(level_file, header=None).to_numpy()
	n_cells = np.sum(levels!=0)
	n_cells_l.append(n_cells)

	means_l.append(np.sum(levels)/n_cells)

	ax = plot.subplot(5, 6, idx+1)
	ax.hist(levels.flatten(), bins=8)
	ax.set_title(f"density = {round(density, 4)}")
	# ax.set_yscale('log')

plot.tight_layout()
plot.savefig("check_hist2.png")


pts_order = np.argsort(density_l)
density_l = np.array(density_l)
means_l = np.array(means_l)
plot.figure(4)
plot.scatter(density_l[pts_order], means_l[pts_order])
plot.xlabel("cell seeding density")
plot.ylabel("global mean signal output")
plot.savefig("check2.png")