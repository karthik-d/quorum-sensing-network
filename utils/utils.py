import json 
import os


def get_config_of_simulation(seeding_src, **kwargs):
	with open(os.path.join("./outputs", seeding_src, "config.json")) as f:
		config = json.loads(json.load(f))
	return config