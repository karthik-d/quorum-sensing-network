from matplotlib import pyplot as plot
from networkx.drawing.nx_agraph import graphviz_layout, to_agraph

import networkx as nx
import numpy as np
import pandas as pd


def get_cytoscape_tables(adjacency_matrix, node_posns=None, return_adjacency_df=False):
	"""
	OUTPUT: (edgetable, nodetable, [adjacency_df]).
	`node_posns` are included as `x` and `y` columns in the output nodetable.
	edge weights are just set to 1 for now.
	"""

	node_posns = np.array(node_posns)
	nodetable = pd.DataFrame(dict(
		x=node_posns[:, 0],
		y=node_posns[:, 1]
	), index=[f"C{i+1}" for i in range(adjacency_matrix.shape[0])])

	adjacency_df = pd.DataFrame(adjacency_matrix.astype(int), 
		index=nodetable.index,
		columns=nodetable.index)

	sources_l = []
	targets_l = []
	weights_l = []
	xs_l, ys_l = [], []
	for col in adjacency_df.columns:
		target_nodes = adjacency_df.index[adjacency_df[col]==1]
		targets_l.extend(list(target_nodes))
		sources_l.extend([col]*len(target_nodes))
		weights_l.extend([1]*len(target_nodes))

	edgetable = pd.DataFrame(dict(
		source=sources_l,
		target=targets_l,
		weight=weights_l))

	ret = (edgetable, nodetable)
	ret = ret + ((adjacency_df,) if return_adjacency_df else tuple())
	return ret
	


def get_graph(adjacency_matrix, hub_nodes=[], node_posns=None):
	"""
	use node posns to position nodes, if supplied; else, use auto-layouting.
	"""

	rows, cols = np.where(adjacency_matrix == 1)
	edges = zip(rows.tolist(), cols.tolist())
	gr = nx.DiGraph()
	gr.add_edges_from(edges)
	
	# make AGraph; set attributes for hub nodes.
	agr = to_agraph(gr)
	agr.graph_attr['size'] = 100
	agr.graph_attr['ratio'] = 1
	for i, node in enumerate(agr.iternodes()):
		node.attr['label'] = f"C{i+1}" 
		node.attr['style'] = 'filled'
		node.attr['width'], node.attr['height'] = 1.5, 1.5
		node.attr['fontsize'] = 80
		if int(str(node)) in hub_nodes:
			node.attr['color'] = 'red'
			node.attr['fillcolor'] = 'red'
			node.attr['shape'] = 'hexagon'
		else:
			node.attr['color'] = 'green'
			node.attr['fillcolor'] = 'green'
			node.attr['shape'] = 'ellipse'
		# set node positions if supplied.
		if node_posns is not None:
			node.attr['pos'] = f"{node_posns[i][0]},{node_posns[i][1]}!"

	for i, edge in enumerate(agr.iteredges()):
		source, target = edge
		edge.attr['dir'] = "both"
		edge.attr['arrowsize'] = 5
		edge.attr['arrowtail'] = "tee"
		if int(str(source)) in hub_nodes:
			edge.attr['color'] = 'red'
			edge.attr['penwidth'] = 2.5
		else:
			edge.attr['color'] = 'black'
			edge.attr['penwidth'] = 2

	_ = agr.layout(prog='sfdp') if node_posns is None else agr.layout(prog='neato')
	return agr


def plot_graph(agr, savepath=None):
	# savepath=None returns as Bytes object.
	agr.graph_attr['dpi'] = '100'
	agr_bytes = agr.draw(savepath, format="png")
	return agr_bytes

