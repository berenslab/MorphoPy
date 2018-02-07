
import copy
import random
import numpy as np
import networkx as nx
from .graph_utils import get_tips, get_nodes


def get_persistence_barcode(G, dist='radDist', axon=True, basal_dendrites=True, apical_dendrites=True):

    if dist == 'radDist':
        f = _radial_dist_to_soma
    else:
        raise NotImplementedError

    nodes = get_nodes(G, type_ix=1)
    if axon:
        nodes += get_nodes(G,type_ix=2)

    if basal_dendrites:
        nodes += get_nodes(G, type_ix=3)

    if apical_dendrites:
        nodes += get_nodes(G, type_ix=4)

    graph = nx.subgraph(G, nodes)
    return _get_persistence_barcode(graph, f)


def _get_persistence_barcode(G, f):
    """
    Creates the persistence barcode for the graph G. The algorithm is taken from
    _Quantifying topological invariants of neuronal morphologies_ from Lida Kanari et al
    (https://arxiv.org/abs/1603.08432).


    :param G: networkx.DiGraph
    :param f: A real-valued function defined over the set of nodes in G.
    :return: numpy.array(Nx2) persistence bar code, holding the birth time and the death time for each leaf and
    branch point in G.
    """
    # Initialization
    L = get_tips(G)
    R = 1

    D = []  # holds persistence barcode
    v = dict()  # holds 'aging' function of visited nodes defined by f

    # active nodes
    A = list(copy.copy(L))

    # set the initial value for leaf nodes
    for l in L:
        v[l] = f(G, l)

    while R not in A:
        for l in A:
            p = G.predecessors(l)[0]
            C = G.successors(p)

            # if all children are active
            if all(c in A for c in C):
                # choose randomly from the oldest children
                c_m = _get_oldest_children(v, C)

                A.append(p)

                for c_i in C:
                    A.remove(c_i)
                    if c_i != c_m:
                        D.append([v[c_i], f(G, p)])
                v[p] = v[c_m]
    D.append([v[R], f(G, R)])
    return np.array(D)


def _radial_dist_to_soma(G,n):
    root_pos = G.node[1]['pos']
    return np.sqrt(np.dot(G.node[n]['pos'] - root_pos, G.node[n]['pos'] - root_pos))


def _get_oldest_children(v, C):
    """

    :param v: dict
    :param C: list. Holds the children ids whose values are evaluated.
    :return: id of oldest child according to dictionary v. If there is no clear maximum the id is chosen randomly across
    the set of oldest.
    """
    age = np.array([v[c] for c in C])
    indices = np.where(age == age[np.argmax(age)])[0]

    return C[random.choice(indices)]
