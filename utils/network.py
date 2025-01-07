from matplotlib import pyplot as plot
from networkx.drawing.nx_agraph import graphviz_layout, to_agraph

import networkx as nx
import numpy as np


def get_graph(adjacency_matrix, hub_nodes=[], node_posns=None):
	"""
	use node posns to position nodes, if supplied; else, use auto-layouting.
	"""

	print(node_posns)
	rows, cols = np.where(adjacency_matrix == 1)
	edges = zip(rows.tolist(), cols.tolist())
	gr = nx.Graph()
	gr.add_edges_from(edges)
	
	# make AGraph; set attributes for hub nodes.
	agr = to_agraph(gr)
	for i, node in enumerate(agr.iternodes()):
		node.attr['label'] = f"C{str(node)}"
		node.attr['style'] = 'filled'
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

	_ = agr.layout(prog='sfdp') if node_posns is None else agr.layout(prog='neato')
	return agr


def plot_graph(agr, savepath=None):
	# savepath=None returns as Bytes object.
	agr.graph_attr['dpi'] = '100'
	agr_bytes = agr.draw(savepath, format="png")
	return agr_bytes

