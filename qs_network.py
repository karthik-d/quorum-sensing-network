import numpy as np
from scipy import signal


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
	
		# indicator matrix denoting cell positions.
		self.cells = np.array([list(map(int, list(row_posns))) for row_posns in cells.split('.')])
		assert self.cells.shape == (self.size, self.size), f"cell positions must match area dim = {self.size}."

		# levles of signaling at each area, per level.
		self.levels = np.zeros((self.domain[-1]-self.domain[0]+1, self.size, self.size))
		# set the first 2D matrix to cell positions -- initial level is 1 at the cells?
		self.levels[0] = self.cells.copy()

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
		returns a matrix encoding the effect of a cell's level distributes over area.
		takes current signal; returns neighborhood effect corresponding to the next signal level.
		"""

		next_signal = self.conc_response_map.get(curr_signal)
		return np.array([
			[max([0, next_signal - abs(x) - abs(y)]) for x in range(-self.size, self.size+1)]
			for y in range(-self.size, self.size+1)
		])


class QSNetworkSimulator:

	def __init__(self, qs_net, obs_duration=10, signaling_interval=1):
		self.net = qs_net 
		self.obs_duration = obs_duration
		self.signaling_interval = signaling_interval

		self.init_simulator()


	def init_simulator(self):
		self.graph = list(range(self.obs_duration))
		self.production_list = list(range(self.obs_duration))

		self.signal_cloud = signal.convolve2d(
			self.net.cells, self.net.get_neighborhood_kernel(curr_signal=1), 
			mode='same', boundary='fill', fillvalue=0)

		print(self.net.cells)
		print(self.net.get_neighborhood_kernel(curr_signal=1))
		print(self.signal_cloud)


	def run_main_qs_cycle(self):

		for clock in range(self.obs_duration):

			# decide whether or not to signal.
			if clock % self.signaling_interval == 0:
				pass 

			# routine, regardless of signaling.
			self.graph[clock] = None


simulator = QSNetworkSimulator(
	qs_net = QSNetwork()
)
