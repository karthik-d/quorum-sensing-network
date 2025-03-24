import json 
import os
import numpy as np


def get_config_of_simulation(seeding_src, **kwargs):
	with open(os.path.join("./outputs", seeding_src, "config.json")) as f:
		config = json.loads(json.load(f))
	return config


def get_nbins_hist(data, bin_size):
	return np.arange(np.min(data), np.max(data)+bin_size, bin_size)