import numpy as np
from matplotlib import pyplot as plot
from scipy import signal


rand_gen = np.random.default_rng()


class QSNetwork:

	def __init__(
		self, size=5, cell_density=0.2, domain_range=(1, 10),
		cells='00000.00100.01110.01100.00000',
		cell_response_mapping=[1, 3, 5, 7, 6, 5, 4, 3, 2, 1],
		conc_response_mapping=[3, 5, 7, 9, 8, 7, 6, 5, 4, 3]
	):
		""" parameters:
		- size: determines the dimension of the square cell matrix.
		- cell_density: determines ...
		- cells: a binary encoding of cell positions.
		- response_mapping: next levels for each corresponding current level (these levels are drawn from `domain`).
		"""

		self.size = size 
		self.cell_density = cell_density

		domain_range = (domain_range[0], domain_range[-1]+1) # adjust to include the last value.
		self.domain = np.array(range(*domain_range))

		# initialize the network.
		self.init_net(cells, cell_response_mapping, conc_response_mapping)
	

	def init_net(self, cells, cell_response_mapping, conc_response_mapping):
	
		# indicator matrix denoting positions of "active" cells.
		self.cells = np.array([list(map(int, list(row_posns))) for row_posns in cells.split('.')])
		assert self.cells.shape == (self.size, self.size), f"cell positions must match area dim = {self.size}."

		## NOTE: `levels` is always indexed through the values of domains. 
		## so it has an additional unused element at idx=0.
		## and it is initialized to have len(domain)+1 size along its first dimension.
		
		# levels of signaling at each posn in area, per level.
		self.levels = np.zeros((self.domain[-1]-self.domain[0]+2, self.size, self.size))
		# set the first 2D matrix to cell positions -- initial level is 1 at the cells?
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
			[max([0, next_signal - abs(x) - abs(y)]) for x in range(-self.size, self.size+1)]
			for y in range(-self.size, self.size+1)
		])


# TODO: 
# - Add a logging data structure for the simulator to track updates made to the network.

class QSNetworkSimulator:

	def __init__(
		self, qs_net, obs_duration=10, signaling_interval=1, 
		level_update_thresh=3, signaling_frac=1,
		verbose=False
	):
		""" parameters:
		- qs_net: the network to simulate; a QSNetwork instance.
		- obs_duration: total simulation length.
		- signaling_interval: frequency of update steps.
		- signaling_frac: number of cells to select as signaling cells in each update step (default: all cells signal).
		- verbose: boolean to indicate if simulation comments should be displayed.
		"""
		
		self.net = qs_net 
		self.obs_duration = obs_duration
		self.signaling_interval = signaling_interval
		self.level_update_thresh = level_update_thresh
		self.signaling_frac = signaling_frac
		self.verbose = verbose

		# initialize simulator by seeding network and running one step.
		self.init_simulation()


	def _convolve_with_neighborhood(self, curr_signal):
		return signal.convolve2d(
			self.net.cells, self.net.get_neighborhood_kernel(curr_signal), 
			mode='same', boundary='fill', fillvalue=0)


	def init_simulation(self):
		self.graph = list(range(self.obs_duration))
		self.production_list = list(range(self.obs_duration))

		# generate the initial signal cloud by running one convolution with neighborhood.
		signal_cloud = self._convolve_with_neighborhood(curr_signal=1)

		# run a single step of simulation.
		self.net.levels = self.step_simulation(signal_cloud)


	def step_simulation(self, signal_cloud, levels=None, thresh=3):
		"""
		INPUTS: signal_cloud, levels (use from net if None), thresh (for level updation).
		OUTPUTS: updated levels.
		> updates `levels` of the network based on current `signal_cloud`.
		> does NOT mutate any instance members.
		"""

		# override levels to use if passed.
		levels = self.net.levels if levels is None else levels
		thresh = self.net.level_update_thresh if thresh is None else thresh

		# copy current levels and modify the new object.
		updated_levels = levels.copy()

		# update-routine for first domain value.
		cellposn_idx_l = np.where(levels[self.net.domain[0]]==1)
		conc_l = signal_cloud[cellposn_idx_l[0], cellposn_idx_l[1]]
		for i, conc in enumerate(conc_l):
			if conc > self.net.domain[0] + thresh:
				updated_levels[self.net.domain[0] + 1, cellposn_idx_l[0][i], cellposn_idx_l[1][i]] = 1
				updated_levels[self.net.domain[0], cellposn_idx_l[0][i], cellposn_idx_l[1][i]] = 0

		# update-routine for last domain value.
		cellposn_idx_l = np.where(levels[self.net.domain[-1]]==1)
		conc_l = signal_cloud[cellposn_idx_l[0], cellposn_idx_l[1]]
		for i, conc in enumerate(conc_l):
			if conc < self.net.domain[-1] + thresh:
				updated_levels[self.net.domain[-1] - 1, cellposn_idx_l[0][i], cellposn_idx_l[1][i]] = 1
				updated_levels[self.net.domain[-1], cellposn_idx_l[0][i], cellposn_idx_l[1][i]] = 0

		# update-routine for rest of the domain values.
		for domain_val in self.net.domain[1:-1]:
			
			cellposn_idx_l = np.where(levels[domain_val]==1)
			conc_l = signal_cloud[cellposn_idx_l[0], cellposn_idx_l[1]]
			for i, conc in enumerate(conc_l):
				
				# conc. is within 13.
				if conc < domain_val + 13:
					updated_levels[domain_val - 1, cellposn_idx_l[0][i], cellposn_idx_l[1][i]] = 1
					updated_levels[domain_val, cellposn_idx_l[0][i], cellposn_idx_l[1][i]] = 0

				# conc. is 13 or more over.
				elif conc > domain_val + 13:
					updated_levels[domain_val + 1, cellposn_idx_l[0][i], cellposn_idx_l[1][i]] = 1
					updated_levels[domain_val, cellposn_idx_l[0][i], cellposn_idx_l[1][i]] = 0

		# return the updated levels matrix.
		return updated_levels.copy()


	def run_qs_simulation(self, obs_duration=None, signaling_interval=None):

		# override simulator defaults if args passed.
		obs_duration = self.obs_duration if obs_duration is None else obs_duration
		signaling_interval = self.signaling_interval if signaling_interval is None else signaling_interval

		_ = print(f"simulating {obs_duration} time steps ...") if self.verbose else None
		# run simulation.
		log = dict()
		plot.clf()
		plot.figure(figsize=(12, 12))

		# plot initial.
		init_network = np.sum([
			self.net.levels[i]*(self.net.cell_response_map.get(i)+1) for i in self.net.domain
		], axis=0) 
		plot.subplot(5, 5, 1)
		plot.title(f"time=0")
		plot.imshow(init_network, vmin=0, vmax=10)
		plot.colorbar()
		
		for time in range(1, obs_duration+1):

			# decide whether or not to signal.
			if time % signaling_interval != 0:
				continue

			# == signaling routine ==

			# select the signaling cells.
			# NOTE: this is where the probability-based selection of signaling cells would go.
			num_cells = self.net.cells.sum()
			# signaling_cells = self.net.cells.copy()
			signaling_cells = np.zeros(self.net.cells.shape)
			signaling_cells[*list(zip(*rand_gen.choice(
				list(zip(*np.nonzero(self.net.cells))), int(self.signaling_frac*num_cells), replace=False
			)))] = 1
			
			# generate signal cloud.
			signal_cloud = np.sum([
				self._convolve_with_neighborhood(curr_signal=i) for i in self.net.domain
			], axis=0) 

			# TODO: vectorize the for-loop using `apply_along_axis` and reshaping the inner 2D matrices in the function.
			# TODO: can also simply replicate `signaling_cells` are perform a matrix operation.
			# set any level that has become -1 to 0. (why though??).
			# reduce the level by 1 for signaling cells -- call it `static_levels`. (why??).
			static_levels = np.zeros(self.net.levels.shape)
			for domain in self.net.domain:
				static_levels[domain] = self.net.levels[domain] - signaling_cells
				static_levels[domain][np.where(static_levels[domain]==-1)] = 0

			# `dynamic_levels` simply recovers `levels`, but has +1 wherever levels were 0 at a signaling cell position.
			# (rationale unclear??)
			dynamic_levels = self.net.levels - static_levels
			# use signaling_cloud and static_levels to run a simulation step.
			dynamic_levels = self.step_simulation(signal_cloud, levels=dynamic_levels)
			self.net.levels = dynamic_levels + static_levels
			
			# log observation.
			# NOTE: (cell_response_map + 1) is essentially `aiDistance`.
			log[time] = np.sum([
				self.net.levels[i]*(self.net.cell_response_map.get(i)+1) for i in self.net.domain
			], axis=0) 

			# plot the logged matrix.
			plot.subplot(5, 5, time+1)
			plot.title(f"time={time}")
			plot.imshow(log[time], vmin=0, vmax=10)
			plot.colorbar()
		
		_ = print("saving plots...") if self.verbose else None
		# save progression.
		plot.tight_layout()
		plot.savefig(f"levels_duration-{obs_duration}_select-{self.signaling_frac}.png", dpi=100)


simulator = QSNetworkSimulator(
	qs_net = QSNetwork(
		cells="00100.00000.01101.00010.01000"
	),
	obs_duration = 24,
	signaling_frac = 1,
	verbose = True
)

simulator.run_qs_simulation()
