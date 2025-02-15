import numpy as np
import os
import pathlib
import math
import json

from matplotlib import pyplot as plot
from scipy import signal
from datetime import datetime


from utils import animation
from utils import network
from utils import seeding
from utils import utils

rand_gen = np.random.default_rng()


class QSNetwork:

	def __init__(self, config):
		""" parameters:
		- cell_area_dim: determines the dimension of the square cell matrix.
		- cell_density: determines ...
		- cell_posn_encoding: a binary encoding of cell positions.
		- negative_feedback: use negative feedback?
		"""

		kwargs = config.copy()
		self.size = tuple(kwargs.get('cell_area_dim', (5, 5)))
		self.cell_density = kwargs.get('cell_density', 0.2) 

		# adjust to include the last value.
		domain_range_raw = kwargs.get('domain_range', (1, 10))
		domain_range = (domain_range_raw[0], domain_range_raw[-1]+1) 
		self.domain = np.array(range(*domain_range))

		# next levels for each corresponding current level (these levels are drawn from `domain`).
		self.has_negative_feedback = kwargs.get('negative_feedback')
		if self.has_negative_feedback:
			cell_response_mapping = kwargs.get('cell_response_mapping', [1, 3, 5, 7, 6, 5, 4, 3, 2, 1])
		else:
			cell_response_mapping = kwargs.get('cell_response_mapping', [1, 2, 3, 4, 5, 6, 7, 7, 7, 7]), 
		conc_response_mapping=[x+2 for x in cell_response_mapping]
		
		# initialize the network.
		self.init_net(
			kwargs.get('cell_posn_encoding', '00000.00100.01110.01100.00000'), 
			cell_response_mapping, conc_response_mapping)
	

	def init_net(self, cells, cell_response_mapping, conc_response_mapping):
	
		# indicator matrix denoting positions of "active" cells.
		self.cells = np.array([list(map(int, list(row_posns))) for row_posns in cells.split('.')])
		assert all(x==y for x, y in zip(self.cells.shape, self.size)), f"cell positions must match area dim = {self.size}."

		## NOTE: `levels` is always indexed through the values of domains. 
		## so it has an additional unused element at idx=0.
		## and it is initialized to have len(domain)+1 size along its first dimension.
		
		# levels of signaling at each posn in area, per level.
		self.levels = np.zeros((self.domain[-1]-self.domain[0]+2, ) + self.size)
		# set the first 2D matrix to cell positions -- initial level is 1 at all the cells.
		self.levels[self.domain[0]] = self.cells.copy()

		# define how current level defines the next level for "cell".
		# equivalent: aiCell.
		self.cell_response_map = dict(zip(
			self.domain, cell_response_mapping
		))

		# define how current level defines the next level for "concentrations".
		# equivalent: aiConcentration.
		self.conc_response_map = dict(zip(
			self.domain, conc_response_mapping
		))


	def get_neighborhood_kernel(self, curr_signal):
		"""
		INPUT: takes current signal; returns neighborhood effect corresponding to the next signal level.
		OUTPUT: matrix encoding the effect of a cell's level distributes over area.
		"""

		next_signal = self.conc_response_map.get(curr_signal)
		return np.array([
			[max([0, next_signal - abs(x) - abs(y)]) for x in range(-self.size[0], self.size[0]+1)]
			for y in range(-self.size[1], self.size[1]+1)
		])


# TODO: 
# - Add a logging data structure for the simulator to track updates made to the network.

class QSNetworkSimulator:

	def __init__(self, qs_net, config):
		""" parameters:
		- qs_net: the network to simulate; a QSNetwork instance.
		- obs_duration: total simulation length.
		- signaling_interval: frequency of update steps.
		- signaling_frac: number of cells to select as signaling cells in each update step (default: all cells signal).
		- verbose: boolean to indicate if simulation comments should be displayed.
		"""
		
		kwargs = config.copy()
		self.config = config.copy()
		self.net = qs_net 
		self.obs_duration = kwargs.get("obs_duration", 10)
		self.signaling_interval = kwargs.get("signaling_interval", 1)
		self.level_update_thresh = kwargs.get("level_update_thresh", 3)
		self.signaling_frac = kwargs.get("signaling_frac", 1)
		self.fixed_signalers = kwargs.get('fixed_signalers', False)
		self.seeding_frac = kwargs.get("cell_seeding_frac", 0.1)
		self.verbose = kwargs.get("verbose", False)

		# initialize simulator by seeding network and running one step.
		self.init_simulation()


	def _convolve_with_neighborhood(self, curr_signal):
		return signal.convolve2d(
			self.net.levels[curr_signal], self.net.get_neighborhood_kernel(curr_signal), 
			mode='same', boundary='fill', fillvalue=0)
	

	@classmethod 
	def find_hub_cells(cls, edge_matrix, hub_cell_thresh=5):
		num_connections_l = np.sum(edge_matrix, axis=1) # sum each row.
		hub_cells = [int(num_connections_l[i] >= hub_cell_thresh) for i in range(edge_matrix.shape[0])]
		return hub_cells


	def get_network_graph(self, simulation_log, time_step=None):
		"""
		Returns: 
		- Edge matrix of constructed graph - 2D, undirected.
		- List of index positions of cells (nodes) on plate.
		Note: Nodes occur in the same relative order in the edge matrix and index positions list.
		"""
		
		# use the last time step by default.
		time_step = max(simulation_log.keys()) if time_step is None else time_step
		graph = simulation_log[time_step]

		n_cells = self.net.cells.sum()
		cellposn_idx_l = np.where(self.net.cells==1)
		edge_matrix = np.zeros((n_cells, n_cells, ))

		assert len(cellposn_idx_l[0]) == n_cells, "[ERROR] number of cells must match size of node positions."

		for i in range(n_cells):
			# the threshold is basically how far this level has "influence" 
			# (directly translates, since level dissipates directly with Manhattan distance).
			max_dist = graph[cellposn_idx_l[0][i], cellposn_idx_l[1][i]]
			for j in range(n_cells):
				if i!=j: 
					dist = abs(cellposn_idx_l[0][i] - cellposn_idx_l[0][j]) + (
						abs(cellposn_idx_l[1][i] - cellposn_idx_l[1][j])
					)	# manhattan distance.
					edge_matrix[i, j] = 1 if (dist<=max_dist) else 0
		
		return edge_matrix, list(zip(*cellposn_idx_l))


	def init_simulation(self):
		self.graph = list(range(self.obs_duration))
		self.production_list = list(range(self.obs_duration))   # not used yet.

		# generate the initial signal cloud by running one convolution with neighborhood.
		signal_cloud = self._convolve_with_neighborhood(curr_signal=1)

		# run a single step of simulation.
		self.net.levels = self.update_signal_levels(signal_cloud, cells_to_update=self.net.cells.copy())


	def update_signal_levels(self, signal_cloud, cells_to_update, levels=None, thresh=3):
		"""
		INPUTS: signal_cloud, cells_to_update, levels (use from net if None), thresh (for level updation).
		OUTPUTS: updated levels.
		> updates `levels` of the network based on current `signal_cloud`.
		> does NOT mutate any instance members.
		"""

		# override `levels` and `thresh` to use if passed.
		levels = self.net.levels if levels is None else levels
		thresh = self.net.level_update_thresh if thresh is None else thresh

		# copy current levels and modify the new object.
		updated_levels = levels.copy()

		# update-routine for cells whose current level is 1.
		# select only cells that should have signals updated.
		cellposn_idx_l = list(zip(*[ 
			(i, j)
			for i, j in zip(*np.where(levels[self.net.domain[0]]==1))
			if cells_to_update[i, j]
		]))
		perceived_signals_l = signal_cloud[cellposn_idx_l[0], cellposn_idx_l[1]] if cellposn_idx_l else []
		for i, perceived_signal in enumerate(perceived_signals_l):
			# if the signal at this cell OVER (min+thresh) --> increase cell's level by 1.
			if perceived_signal > self.net.domain[0] + thresh:
				updated_levels[self.net.domain[0] + 1, cellposn_idx_l[0][i], cellposn_idx_l[1][i]] = 1
				updated_levels[self.net.domain[0], cellposn_idx_l[0][i], cellposn_idx_l[1][i]] = 0

		# update-routine for cells whose current level is 10.
		cellposn_idx_l = list(zip(*[ 
			(i, j)
			for i, j in zip(*np.where(levels[self.net.domain[-1]]==1))
			if cells_to_update[i, j]
		]))
		perceived_signals_l = signal_cloud[cellposn_idx_l[0], cellposn_idx_l[1]] if cellposn_idx_l else []
		for i, perceived_signal in enumerate(perceived_signals_l):
			# if the signal at this cell UNDER (max+thresh) --> decrease cell's level by 1.
			if perceived_signal < self.net.domain[-1] + thresh:
				updated_levels[self.net.domain[-1] - 1, cellposn_idx_l[0][i], cellposn_idx_l[1][i]] = 1
				updated_levels[self.net.domain[-1], cellposn_idx_l[0][i], cellposn_idx_l[1][i]] = 0

		# update-routine for cells at all other levels.
		for my_val in self.net.domain[1:-1]:
			
			cellposn_idx_l = list(zip(*[ 
				(i, j)
				for i, j in zip(*np.where(levels[my_val]==1))
				if cells_to_update[i, j]
			]))
			perceived_signals_l = signal_cloud[cellposn_idx_l[0], cellposn_idx_l[1]] if cellposn_idx_l else []
			for i, perceived_signal in enumerate(perceived_signals_l):
				
				# if the signal at this cell UNDER (curr+thresh) --> decrease cell's level by 1.
				if perceived_signal < my_val + thresh:
					updated_levels[my_val - 1, cellposn_idx_l[0][i], cellposn_idx_l[1][i]] = 1
					updated_levels[my_val, cellposn_idx_l[0][i], cellposn_idx_l[1][i]] = 0

				# if the signal at this cell is OVER (curr+thresh) --> increase cell's level by 1.
				elif perceived_signal > my_val + thresh:
					updated_levels[my_val + 1, cellposn_idx_l[0][i], cellposn_idx_l[1][i]] = 1
					updated_levels[my_val, cellposn_idx_l[0][i], cellposn_idx_l[1][i]] = 0

		# return the updated levels matrix.
		return updated_levels.copy()


	def run_qs_simulation(self, obs_duration=None, signaling_interval=None, save_outputs=False):
		"""
		driver to run simulation with all set parameters.
		"""

		# -- arg processing --
		# override simulator defaults if args passed.
		obs_duration = self.obs_duration if obs_duration is None else obs_duration
		signaling_interval = self.signaling_interval if signaling_interval is None else signaling_interval
		# -----

		# -- run simulation --
		_ = print(f"simulating {obs_duration} time steps ...") if self.verbose else None
		log = dict(cloud=dict(), levels=dict())

		# save initial.
		init_levels = np.sum([
			self.net.levels[i]*(self.net.cell_response_map.get(i)+1) for i in self.net.domain
		], axis=0) 
		log.get("levels")[0] = init_levels
		log.get("cloud")[0] = self._convolve_with_neighborhood(curr_signal=1)

		# if fixed_signalers is True, then pre-select signaling subsets.
		# shuffle the list and select cells with overlap.
		if self.fixed_signalers:
			signaling_set_idx = 0
			signaling_sets = []
			n_cells = self.net.cells.sum()
			n_cells_per_set = int(n_cells*self.signaling_frac)
			# overlap between subsets to have common signalers.
			n_overlaps = int(0.5*n_cells_per_set)   # overlapping half the set; change if required.
			all_cell_posns = rand_gen.permutation(list(zip(*np.nonzero(self.net.cells))))
			for i in range(0, n_cells-n_overlaps, n_cells_per_set-n_overlaps):
				responding_cells = np.zeros(self.net.cells.shape)
				responding_cells[*list(zip(*all_cell_posns[i:i+n_cells_per_set]))] = 1
				signaling_sets.append(responding_cells)

		for time in range(1, obs_duration+1):

			# decide whether or not to signal.
			if time % signaling_interval != 0:
				continue

			# -- signaling routine --

			# alternate within fixed_signalers OR select the signaling cells at random.
			if self.fixed_signalers:
				responding_cells = signaling_sets[signaling_set_idx]
				signaling_set_idx = (signaling_set_idx+1)%len(signaling_sets)
			else:
				n_cells = self.net.cells.sum()
				responding_cells = np.zeros(self.net.cells.shape)
				responding_cells[*list(zip(*rand_gen.choice(
					list(zip(*np.nonzero(self.net.cells))), int(self.signaling_frac*n_cells), replace=False
				)))] = 1
			
			# generate signal cloud.
			signal_cloud = np.sum([
				self._convolve_with_neighborhood(curr_signal=i) for i in self.net.domain
			], axis=0) 

			# `static_levels` defines levels of all cells that won't update in this time step.
			static_levels = self.net.levels - np.tile(
				responding_cells, reps=((self.net.levels.shape[0],) + (1,)*responding_cells.ndim)) 
			static_levels[static_levels==-1] = 0
			# `dynamic_levels` defines levels of all cells that will update in this time step.
			dynamic_levels = self.net.levels - static_levels

			# update cell levels.
			self.net.levels = self.update_signal_levels(
				signal_cloud, cells_to_update=responding_cells, levels=self.net.levels
			)
			
			# log observation.
			# NOTE: (cell_response_map[x] + 1) is essentially `aiDistance`.
			log.get("levels")[time] = np.sum([
				self.net.levels[i]*(self.net.cell_response_map.get(i)+1) for i in self.net.domain
			], axis=0) 
			log.get("cloud")[time] = signal_cloud.copy()


		# optionally, save time evolution outputs.
		if save_outputs:
			_ = print("saving plots...") if self.verbose else None
			self.save_outputs(log, obs_duration, self.config, os.path.join(
				"./outputs", f"""{
					datetime.now().strftime("%m%d%Y%H%M%S")}_size-{
						'x'.join(map(str, self.net.size))}_select-{round(self.signaling_frac, 2)}_seed-{
							round(self.seeding_frac, 4)}{'_noneg' if not self.net.has_negative_feedback else ''}"""))
		else:
			_ = print("done running. not saving.") if self.verbose else None


		return log
	

	def save_outputs(self, log, obs_duration, config, savedir):
		pathlib.Path(savedir).mkdir(
			exist_ok=False, parents=True
		)

		# plot the logged matrix.
		levels_l = []
		clouds_l = []
		delta_clouds_l = []
		graphs_l = []
		subplot_dim = math.ceil((obs_duration+1)**0.5)
		for time in range(obs_duration+1):	# to include plot of time=0.
			levels_l.append(log["levels"][time])
			clouds_l.append(log["cloud"][time])

			# compute and save how the cloud changed in this time step.
			cloud_delta = (log["cloud"][time] - log["cloud"][time-1]) if (
				time>0) else np.zeros(log["cloud"][time].shape)
			delta_clouds_l.append(cloud_delta)

			# plot each time step of `levels`.
			plot.figure(1, figsize=(30, 30))
			plot.subplot(subplot_dim, subplot_dim, time+1)
			plot.title(f"levels; time={time}")
			plot.imshow(log["levels"][time], vmin=0, vmax=10)
			plot.colorbar()

			# plot each time step of `cloud`.
			plot.figure(2, figsize=(30, 30))
			plot.subplot(subplot_dim, subplot_dim, time+1)
			plot.title(f"cloud; time={time}")
			plot.imshow(log["cloud"][time], vmin=0, vmax=10)
			plot.colorbar()

			# plot each delta of `cloud`.
			plot.figure(3, figsize=(30, 30))
			plot.subplot(subplot_dim, subplot_dim, time+1)
			plot.title(f"cloud delta; time={time}")
			plot.imshow(cloud_delta, vmin=0, vmax=10)
			plot.colorbar()

			# make graphs at each time step.
			edge_matrix, node_posns = self.get_network_graph(log["levels"], time_step=time)
			hub_cells = self.find_hub_cells(edge_matrix, hub_cell_thresh=5)
			_ = print(f"no. of hubs at t={time}: {sum(hub_cells)}") if self.verbose else None

			# draw and save as np objects.
			agr = network.get_graph(edge_matrix, hub_nodes = np.where(np.array(hub_cells)==1)[0])
			agr = graph = network.get_graph(
				edge_matrix, 
				hub_nodes = np.where(np.array(hub_cells)==1)[0],   # returns tuple of list of indices; unpack.
				node_posns = node_posns
			)
			agr_bytes = network.plot_graph(agr, savepath=None)
			graphs_l.append(animation.img_bytes2array(agr_bytes))

	
		# save progression.
		plot.figure(1)
		plot.savefig(os.path.join(savedir, f"levels_time-evolution.png"), dpi=100)

		plot.figure(2)
		plot.savefig(os.path.join(savedir, f"cloud_time-evolution.png"), dpi=100)

		plot.figure(3)
		plot.savefig(os.path.join(savedir, f"cloud-deltas_time-evolution.png"), dpi=100)

		# make animations from saved plots.
		animation.save_animation(
			data = [
				levels_l,
				clouds_l,
				delta_clouds_l,
			],
			save_path = os.path.join(savedir, f"levels_cloud_time-evolution.gif")
		)

		# save final graph.
		plot.figure(3, figsize=(30, 30))
		plot.clf()
		# set the initial image.
		plot.imshow(graphs_l[-1], vmin=0, vmax=10)
		plot.savefig(os.path.join(savedir, f"graph_final.png"), dpi=100)

		# save graphs as an animation.
		animation.save_animation(
			data = [graphs_l],
			save_path = os.path.join(savedir, f"graph_time-evolution.gif")
		)

		# save configs.
		with open(os.path.join(savedir, "config.json"), 'w') as f:
			json.dump(json.dumps(config, indent='\t'), f)


# -----
# Template run parameters.
# -----

## 1. basic 5x5 array.
# simulator = QSNetworkSimulator(
# 	qs_net = QSNetwork(
# 		cells="00100.00000.01101.00010.01000"
# 	),
# 	obs_duration = 10,
# 	signaling_frac = 1,
# 	verbose = True
# )

## 2. constant seeding.
simulation_config = dict(
	# network params.
	cell_seeding_frac = 0.0333,
	cell_area_dim = (50, 50),
	negative_feedback = True,

	# set simulation id to load seeding from; None for random.
	seeding_src = None,
	# seeding_src = "02152025033708_size-50x50_select-1_seed-0.0667", 	# 6.67%
	# seeding_src = "02152025042853_size-50x50_select-1_seed-0.0333",  	# 3.33%

	# params for graded seeding; set to `None` if using uniform seeding.
	seeding_transition_frac = None,
	n_seeding_transitions = None,

	# simulator params.
	obs_duration = 99,		# set as (perfect_sq - 1) for good formatting.
	signaling_frac = 0.5,

	# when True, cells are divided (based on signaling_frac) into pre-defined sets; 
	# during updation, a set is chosen cyclically to respond.
	fixed_signalers = False,

	# other params 
	verbose = True
)

## 3. graded seeding.
# simulation_config = dict(
# 	# network params.
# 	cell_seeding_frac = 0.02,
# 	cell_area_dim = (96, 50),
# 	negative_feedback = True,

# 	# params for graded seeding; set to `None` if using uniform seeding.
# 	seeding_transition_frac = 0.005,
# 	n_seeding_transitions = 7,

# 	# simulator params.
# 	obs_duration = 24,		# set as (perfect_sq - 1) for good formatting.
# 	signaling_frac = 1,

# 	# other params 
# 	verbose = True
# )


# -----
# cell seeding options.
# -----

if simulation_config.get('seeding_src') is not None:
	## 1. Load seeding and related params from parameter file (using simulation id).
	config = utils.get_config_of_simulation(**simulation_config)
	# update seeding params based on source.
	simulation_config.update({k: config[k] for k in [
		'cell_posn_encoding', 'cell_seeding_frac', 'cell_area_dim',
		'seeding_transition_frac', 'n_seeding_transitions'
	]})
else:
	## 2. Randomly seed cells.
	seeding_func = seeding.graded_density_array if (
		simulation_config.get('n_seeding_transitions', False)) else seeding.uniform_density_array
	simulation_config.update(dict(
		cell_posn_encoding = seeding_func(**simulation_config)
	))


# -----
# run simulation.
# -----
simulator = QSNetworkSimulator(
	qs_net = QSNetwork(simulation_config),
	config = simulation_config)
log = simulator.run_qs_simulation(
	save_outputs = True
)


# edge_matrix, node_posns = simulator.get_network_graph(log["levels"])
# hub_cells = QSNetworkSimulator.find_hub_cells(edge_matrix, hub_cell_thresh=5)
# print(f"No. of Hubs: {sum(hub_cells)}")

# graph = network.get_graph(
# 	edge_matrix, 
# 	hub_nodes = np.where(np.array(hub_cells)==1)[0],   # returns tuple of list of indices; unpack.
# 	node_posns = node_posns
# )
# network.plot_graph(graph, savepath = 'check.png')



# =======
### previous graphing logic.
# make graphs at each time step.
# edge_matrix, node_posns = self.get_network_graph(log["levels"], time_step=time)
# hub_cells = self.find_hub_cells(edge_matrix, hub_cell_thresh=5)
# _ = print(f"no. of hubs at t={time}: {sum(hub_cells)}") if self.verbose else None

# # draw and save as np objects.
# agr = network.get_graph(edge_matrix, hub_nodes = np.where(np.array(hub_cells)==1)[0])
# agr_bytes = network.plot_graph(agr, savepath=None) 
# graphs_l.append(animation.img_bytes2array(agr_bytes))