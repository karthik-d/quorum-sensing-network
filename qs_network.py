class QSNetwork:

	def __init__(self, area=5, cell_density=0.2):
		""" Parameters:
		- area: determines the dimension of the square cell matrix.
		- cell_density: determines ...
		"""

		self.area = area 
		self.cell_density = cell_density
		
	
	def init_net(self):
		pass 



class QSNetworkSimulator:

	def __init__(self, qs_net, obs_duration=10, signaling_interval=1):
		self.net = qs_net 
		self.obs_duration = obs_duration
		self.signaling_interval = signaling_interval


	def init_simulator(self):
		self.graph = [x for x in range(self.obs_duration)]


	def run_main_qs_cycle(self):

		for clock in range(self.obs_duration):

			# decide whether or not to signal.
			if clock % self.signaling_interval == 0:
				pass 

			# routine, regardless of signaling.
			self.graph[clock] = 