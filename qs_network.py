import numpy as np


class QSNetwork:

	def __init__(
		self, area_dim=5, cell_density=0.2, domain_range=(1, 10),
		cells='00000.00100.01110.01100.00000'):
		""" Parameters:
		- area_dim: determines the dimension of the square cell matrix.
		- cell_density: determines ...
		"""

		self.area_dim = area_dim 
		self.cell_density = cell_density
		self.domain = np.array(range(*domain_range))

		# initialize the network.
		self.init_net(cells)
	

	def init_net(self, cells):
	
		# indicator matrix denoting cell positions.
		self.cells = np.array(map(int, [list(row_posns) for row_posns in cells.split('.')]))
		assert self.cells.shape == (self.area_dim, self.area_dim), f"cell positions must match area dim = {self.area_dim}."

		# levles of signaling at each area, per level.
		self.levels = np.zeros((self.domain[-1]-self.domain[0]+1, self.area_dim, self.area_dim))
		# set the first 2D matrix to cell positions -- initial level is 1 at the cells?
		self.levels[0] = self.cells.copy()



class QSNetworkSimulator:

	def __init__(self, qs_net, obs_duration=10, signaling_interval=1):
		self.net = qs_net 
		self.obs_duration = obs_duration
		self.signaling_interval = signaling_interval


	def init_simulator(self):
		self.graph = list(range(self.obs_duration))
		self.production_list = list(range(self.obs_duration))


	def run_main_qs_cycle(self):

		for clock in range(self.obs_duration):

			# decide whether or not to signal.
			if clock % self.signaling_interval == 0:
				pass 

			# routine, regardless of signaling.
			self.graph[clock] = None


	@classmethod
	def get_neighborhood_effect(cls, signal_level, area=None):
		"""
		Returns a matrix encoding the effect of a cell's level distributes over area.
		"""

		return np.array([
			[max([0, signal_level - abs(x) - abs(y)]) for x in range(-area, area+1)]
			for y in range(-area, area+1)
		])