import networkx as nx
import numpy as np


def get_tips(G):
    """
    Returns the indices of the leaf nodes of graph G.
    :param G : (networkx.DiGraph) graph
    :return: numpy.array(dtype=int) with indices of leaf nodes in G.
    """
    E = G.edge
    return np.array([e for e in E if E[e] == {}])


def get_branchpoints(G):
    """
    Returns the indices of the branch nodes of graph G.
    :param G: (networkx.DiGraph) graph
    :return: numpy.array(dtype=int) with indices of branch nodes in graph.
    """
    bp_indx = np.where(np.array(np.sum(nx.adjacency_matrix(G), axis=1)).flatten() > 1)
    return np.array(G.nodes())[bp_indx]


def get_nodes(G, type_ix=None, data=False):

    if type_ix is None:
        nodes = G.nodes(data=data)
    else:
        nodes = [k for k in G.node if G.node[k]['type'] == type_ix]
    return nodes


def get_root(G):
    try:
        root = np.min(get_nodes(G, type_ix=1))
    except (ValueError, KeyError):
        print('No node is attributed to being the soma. Returning the smallest node id.')
        root = np.min(G.nodes())

    return root
