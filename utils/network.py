from matplotlib import pyplot as plot
from networkx.drawing.nx_agraph import graphviz_layout, to_agraph

import networkx as nx
import numpy as np


def get_graph_object(adjacency_matrix):
    rows, cols = np.where(adjacency_matrix == 1)
    edges = zip(rows.tolist(), cols.tolist())
    gr = nx.Graph()
    gr.add_edges_from(edges)
    return gr


def plot_graph_with_labels(adjacency_matrix, savepath, hub_nodes=[], labels=None):
	# TODO: fix labeling key-error issue.

	gr = get_graph_object(adjacency_matrix)
	fig = plot.figure()
	# nx.draw(gr, node_size=50, labels=labels, with_labels=True)
	# nx.draw(gr, node_size=10, with_labels=False)
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

	agr.layout(prog='sfdp')
	agr.draw(savepath)
	return fig

