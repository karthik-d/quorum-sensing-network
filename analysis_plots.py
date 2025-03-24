import os
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plot 


simulation_id = "03242025054229_size-50x50_select-1_seed-0.0667"
nodetable_df = pd.read_csv(os.path.join("./outputs", simulation_id, f"{simulation_id}_nodetable.csv"),
	index_col=0)

# plotting params.
bin_size = 8

plot.figure(figsize=(10, 4))
# plot 1: histogram of out-degrees.
ax = plot.subplot(1, 3, 1)
plot.hist(nodetable_df["outdegree"], bins=bin_size)
plot.xlabel("outgoing nodes")
plot.ylabel("frequency")

# plot 2: histogram of in-degrees.
ax = plot.subplot(1, 3, 2)
plot.hist(nodetable_df["indegree"], bins=bin_size)
plot.xlabel("incoming nodes")
plot.ylabel("frequency")

# plot 3: (log-log) of out-degree.
ax = plot.subplot(1, 3, 3)
counts, bins = np.histogram(nodetable_df["outdegree"], bins=bin_size)
bin_centers = (bins[:-1] + bins[1:]) / 2
plot.scatter(bin_centers, counts)
plot.xlabel("outgoing nodes")
plot.ylabel("rel. frequency")
ax.set_xscale("log")
ax.set_yscale("log")

plot.savefig(os.path.join("./outputs", simulation_id, f"{simulation_id}_plots.png"), dpi=100)

