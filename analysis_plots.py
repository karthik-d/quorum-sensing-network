import os
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plot 

from utils import utils

"""
# for 50x50.
# simulation_id = "02152025042853_size-50x50_select-1_seed-0.0333"  	# 3.33% - base.
# simulation_id = "03242025054236_size-50x50_select-1_seed-0.0333"  	# 3.33% - base.
# simulation_id = "02152025033708_size-50x50_select-1_seed-0.0667"	 	# 6.67% - base.
# simulation_id = "03242025054229_size-50x50_select-1_seed-0.0667" 		# 6.67% - base.
"""

## for 100x100.
# simulation_id = "03242025061244_size-100x100_select-1_seed-0.025"			# 2.50% - base.
# simulation_id = "03242025062513_size-100x100_select-0.4_seed-0.025"		# 2.50 - 0.4f.
# simulation_id = "03242025063109_size-100x100_select-0.3_seed-0.025"		# 2.50 - 0.3f.
# simulation_id = "03242025062438_size-100x100_select-0.2_seed-0.025"		# 2.50 - 0.2f.

# simulation_id = "03242025060053_size-100x100_select-1_seed-0.0333"		# 3.33% - base.
# simulation_id = "03242025061908_size-100x100_select-0.4_seed-0.0333"		# 3.33 - 0.4f.
# simulation_id = "03242025064846_size-100x100_select-0.3_seed-0.0333"		# 3.33 - 0.3f.
# simulation_id = "03242025061850_size-100x100_select-0.2_seed-0.0333"		# 3.33 - 0.2f.

# simulation_id = "03242025060107_size-100x100_select-1_seed-0.0667"		# 6.67% - base.	
# simulation_id = "03242025062007_size-100x100_select-0.4_seed-0.0667"		# 6.67 - 0.4f.	
# simulation_id = "03242025061947_size-100x100_select-0.3_seed-0.0667"		# 6.67 - 0.3f.	
# simulation_id = "03242025061955_size-100x100_select-0.2_seed-0.0667"		# 6.67 - 0.2f.

"""
## for 150x150.
# simulation_id = "03242025065647_size-150x150_select-1_seed-0.025"			# 2.50 - base.
# simulation_id = "03242025075315_size-150x150_select-0.4_seed-0.025"		# 2.50 - 0.4f.
# simulation_id = "03242025075319_size-150x150_select-0.3_seed-0.025"		# 2.50 - 0.3f.
# simulation_id = "03242025075344_size-150x150_select-0.2_seed-0.025"		# 2.50 - 0.2f.

# simulation_id = "03242025065738_size-150x150_select-1_seed-0.0333"		# 3.33 - base.
# simulation_id = ""		# 3.33 - 0.4f.

# simulation_id = "03242025065720_size-150x150_select-1_seed-0.0667"		# 6.67 - base.
"""

# read data file.
nodetable_df = pd.read_csv(os.path.join("./outputs", simulation_id, f"{simulation_id}_nodetable.csv"),
	index_col=0)


# input plotting params.
bin_size = 1


plot.figure(figsize=(20, 10))
# plot 1: histogram of out-degrees.
ax = plot.subplot(1, 3, 1)
data = nodetable_df["outdegree"]
plot.hist(nodetable_df["outdegree"], bins=utils.get_nbins_hist(data, bin_size))
plot.xlabel("outgoing nodes")
plot.ylabel("frequency")

# plot 2: histogram of in-degrees.
ax = plot.subplot(1, 3, 2)
data = nodetable_df["indegree"]
plot.hist(data, bins=utils.get_nbins_hist(data, bin_size))
plot.xlabel("incoming nodes")
plot.ylabel("frequency")

# plot 3: (log-log) of out-degree.
ax = plot.subplot(1, 3, 3)
data = nodetable_df["outdegree"]
counts, bins = np.histogram(data, bins=utils.get_nbins_hist(data, bin_size))
print(len(bins), max(data), min(data))
bin_centers = (bins[:-1] + bins[1:]) / 2
plot.scatter(bin_centers, counts)
plot.xlabel("outgoing nodes")
plot.ylabel("rel. frequency")
ax.set_xscale("log")
ax.set_yscale("log")

plot.savefig(os.path.join("./outputs", simulation_id, f"{simulation_id}_plots.png"), dpi=100)

