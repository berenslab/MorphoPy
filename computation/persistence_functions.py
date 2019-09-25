import networkx as nx
import numpy as np

# Array for selecting the functions in generating persistence
functions = ["None",  "radial_distance", "height", "path_length", "branch_order"]


# A selection of filter functions
def radial_distance(G, u, v):
    n = G.node[u]['pos']
    r = G.node[v]['pos']

    return np.sqrt(np.dot(n - r, n - r))


def height(G, u, v):
    n = G.node[u]['pos']
    r = G.node[v]['pos']
    return np.abs((n - r))[2]


def path_length(G, u, v):
    return nx.shortest_path_length(G, v, u, weight='path_length')


def branch_order(G, u, v):

    if u == v:
        bo = 0
    else:
        path = nx.shortest_path(G, v, u)
        bo = np.sum(np.array(list(nx.degree(G, path).values())) > 2)
    return bo
