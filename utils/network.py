from matplotlib import pyplot as plot
from networkx.drawing.nx_agraph import graphviz_layout, to_agraph

import networkx as nx
import numpy as np


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

