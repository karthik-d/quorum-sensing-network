import numpy as np 


def encode_as_string(cell_array):
	return ".".join(["".join(map(str, row)) for row in cell_array])


def uniform_density_array(cell_area_dim, cell_seeding_frac, **kwargs):
	# ignore extra kwargs.
	cells = np.random.choice([0, 1], size=cell_area_dim, p=[1-cell_seeding_frac, cell_seeding_frac])
	print("actual seeding fraction:", cells.sum()/np.prod(cell_area_dim))
	return encode_as_string(cells)
	

def graded_density_array(cell_area_dim, cell_seeding_frac, **kwargs):
	"""
	- seeding_transition_frac: by how much should seeding fraction be stepped.
	- n_seeding_transitions: number of transitions to include; half below and half above `cell_seeding_frac`.
	"""
	transition_frac = kwargs.get("seeding_transition_frac", 0.05)
	n_transitions = kwargs.get("n_seeding_transitions", 4)
	
	sub_area_dim = (cell_area_dim[0]//(n_transitions+1), cell_area_dim[1])
	lower_frac = cell_seeding_frac - (n_transitions*transition_frac)//2
	sub_cells_l = [
		uniform_density_array(sub_area_dim, frac)
		for frac in [lower_frac] + [
			(lower_frac + i*transition_frac) for i in range(1, n_transitions+1)]
	]
	return encode_as_string(np.vstack(sub_cells_l))
	