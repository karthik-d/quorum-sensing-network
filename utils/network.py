from matplotlib import pyplot as plot

import networkx as nx
import numpy as np


def get_graph_object(adjacency_matrix):
    rows, cols = np.where(adjacency_matrix == 1)
    edges = zip(rows.tolist(), cols.tolist())
    gr = nx.Graph()
    gr.add_edges_from(edges)
    return gr


def plot_graph_with_labels(adjacency_matrix, labels):
    gr = get_graph_object(adjacency_matrix)
    fig = plot.figure()
    nx.draw(gr, node_size=500, labels=labels, with_labels=True)
    return fig

