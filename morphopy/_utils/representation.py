
import copy
import random
import numpy as np
import pandas as pd

from .graph_utils import get_tips, get_root


def get_persistence_barcode(G, dist='radDist'):

    if dist == 'radDist':
        f = _radial_dist_to_soma
    else:
        raise NotImplementedError
    return _get_persistence_barcode(G, f)


def _get_persistence_barcode(G, f):
    """
    Creates the persistence barcode for the graph G. The algorithm is taken from
    _Quantifying topological invariants of neuronal morphologies_ from Lida Kanari et al
    (https://arxiv.org/abs/1603.08432).


    :param G: networkx.DiGraph
    :param f: A real-valued function defined over the set of nodes in G.
    :return: pandas.DataFrame with entries node_id | type | birth | death . Where birth and death are defined in
    distance from soma according to the distance function f.
    """
    # Initialization
    L = get_tips(G)
    R = get_root(G)

    D = dict(node_id=[], type=[], birth=[], death=[])   # holds persistence barcode
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
                age = np.array([v[c] for c in C])
                indices = np.where(age == age[np.argmax(age)])[0]
                c_m = C[random.choice(indices)]

                A.append(p)

                for c_i in C:
                    A.remove(c_i)
                    if c_i != c_m:
                        D['node_id'].append(c_i)
                        D['type'].append(G.node[c_i]['type'])
                        D['birth'].append(v[c_i])
                        D['death'].append(f(G, p))
                v[p] = v[c_m]
    D['node_id'].append(R)
    D['type'].append(1)
    D['birth'].append(v[R])
    D['death'].append(f(G, R))
    return pd.DataFrame(D)


def _radial_dist_to_soma(G, n):
    root_pos = G.node[1]['pos']
    return np.sqrt(np.dot(G.node[n]['pos'] - root_pos, G.node[n]['pos'] - root_pos))
