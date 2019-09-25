import copy
import sys
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from os import makedirs
from os.path import exists
from scipy.io import savemat
from scipy.interpolate import interp1d
from scipy import stats
from shapely.geometry import MultiLineString, LineString, Point
from itertools import combinations
from sklearn.decomposition import PCA, FastICA
from neurontree.utils import angle_between, get_rotation_matrix, rotationMatrixToEulerAngles

sys.setrecursionlimit(100000)
matplotlib.rcParams.update({'font.size': 14})

class NeuronTree:

    # resample nodes along the edges of G in distance d, return array of 3D positions
    @staticmethod
    def resample_nodes(G, d):
        P = []
        for (u, v, edata) in G.edges(data=True):
            n1 = G.node[u]
            n2 = G.node[v]
            e = edata['euclidean_dist']
            m = n2['pos'] - n1['pos']
            m /= np.linalg.norm(m)

            g = lambda x: x * m + n1['pos']

            for a in range(1, int(e / d)):
                P.append(g(a * d))

        P += list(nx.get_node_attributes(G, 'pos').values())
        return np.array(P)

    # creates a networkX Tree out of a swc file.
    # scaling denotes the conversion factor needed to convert the units given in swc file to microns
    # soma_rad denotes the radius of the soma given in microns
    def __init__(self, swc=None, scaling=1., node_data=[], edge_data=[], graph=None, post_process=True, nxversion=1):
        """
        Creates a NeuronTree object that contains a networkx.DiGraph with node attributes 'pos' [x,y,z],
        'type' [1:soma,2:axon,3: dendrite], 'radius' and edge attributes 'euclidean_dist' and 'path_length'.
        It can be created from an swc file (as pandas.DataFrame or as numpy,ndarray) from a list of nodes and
        edges or from a given networkx.DiGraph.
        :param swc: pandas.DataFrame or numpy.ndarray containing the data of an swc file.
        :param scaling: Scaling of the coordinates in swc file to microns. Sometimes the coordinates are given in voxels
        or other units.
        :param soma_rad: radius of the soma in microns.
        :param node_data: list of nodes with node id and dictionary of node attributes
        :param edge_data: list of edges with dicionary of edge attributes
        :param graph: networkx.DiGraph. It is assumed that the graph contains the required labels.
        :param post_process: boolean. Determines whether the Neuron is post-processed or not. If set to True the soma
        points are merged into one single point, the branch types are unified and short branches that originate from the
        soma and are labelled as axon will get removed.
        """
        # set version of networkX
        self._nxversion = nxversion
        # initialize tree DIRECTED
        if graph:
            G = graph
        else:
            G = nx.DiGraph()

            if swc is not None:
                node_keys = ['pos', 'type', 'radius']

                if type(swc) == pd.DataFrame:

                    # sort out node data
                    n = swc['n'].values # get node ids
                    pos = np.array([swc['x'].values, swc['y'].values, swc['z'].values]).T / scaling
                    radius = swc['radius'].values / scaling
                    t = swc['type'].values
                    pid = swc['parent'].values

                elif type(swc) == np.ndarray:
                    n = swc['n']
                    pos = np.array([swc['x'], swc['y'], swc['z']]).T / scaling
                    t = swc['type']
                    radius = swc['radius'] / scaling
                    pid = swc['parent']

                else:
                    raise ValueError('type of swc representation unknown')

                # create node data
                t[pid == -1] = 1
                # create a list of nodes of the form [(node_id, {'pos': [x,y,z], 'type': t, 'radius' : r}),...]
                node_data = list(zip(n,
                                     [dict(zip(node_keys, [pos[ix], t[ix], radius[ix]])) for ix in range(pos.shape[0])]))

                # create edge data
                n_ = n.tolist()
                parent_idx = [n_.index(pid[ix]) for ix in range(1,len(pid))]
                # calculate euclidean distance between end points of edge
                ec = np.sqrt(np.sum((pos[parent_idx] - pos[1:]) ** 2, axis=1))
                edge_keys = ['euclidean_dist', 'path_length']

                # create a list of edges of the form [(e1,e2, {'euclidean_dist': ec, 'path_length': pl}), ..]
                edge_data = list(zip(pid[1:], n[1:],
                                     [dict(zip(edge_keys, [ec[ix], ec[ix]])) for ix in range(ec.shape[0])]))

            G.add_nodes_from(node_data)
            G.add_edges_from(edge_data)

        self._G = G

        if G.nodes():

            self._remove_redundant_nodes()
            self._make_tree()   # needed to access the functions predecessor and successor
            if 'type' in self.get_node_attributes() and post_process:
                self._merge_roots_by_type()
                self._unify_type()
                self._clean_axon()

    def _merge_roots_by_type(self):
        """
        Starts at the root node in the soma and successively merges all nodes labeled as soma (type=1) into one node.
        This results in a neuron with one single soma node. The manipulation is done inplace.

        """

        R = self._G
        nodes_to_merge = True
        root_ix = self.get_root()
        root = R.node[root_ix]

        while nodes_to_merge:
            nodes_to_merge = False
            for succ in R.successors(root_ix):
                s = R.node[int(succ)]

                if s['type'] == root['type']:
                    nodes_to_merge = True
                    for e in R.successors(succ):
                        n2 = R.node[int(e)]
                        d = np.sqrt(np.sum((root['pos'] - n2['pos']) ** 2))
                        R.add_edge(root_ix, e, euclidean_dist=d, path_length=d)
                    R.remove_node(succ)
        self._G = R

    def _remove_redundant_nodes(self):
        """
        Remove redundant nodes from the NeuronTree. A node is considered redundant if the edge between two nodes has a
        euclidean distance = 0, so the node has the same 3D position as its predecessor.
        """
        # get the nodes whose edge between them has distance = 0.
        nodeindices = np.array(list(nx.get_edge_attributes(self._G, 'euclidean_dist').values())) == 0
        edgelist = list(np.array(list(nx.get_edge_attributes(self._G, 'euclidean_dist').keys()))[nodeindices])
        while edgelist:
            predecessor, redundantNode = edgelist.pop()

            n1 = self._G.node[predecessor]
            successors = self._G.successors(redundantNode)

            # connect edges across redundant nodes
            for succ in successors:
                n2 = self._G.node[succ]
                d = np.sqrt(np.sum((n1['pos'] - n2['pos']) ** 2))
                self._G.add_edge(predecessor, succ, euclidean_dist=d, path_length=d)

            # remove redundant node from graph
            self._G.remove_node(redundantNode)
            nodeindices = np.array(list(nx.get_edge_attributes(self._G, 'euclidean_dist').values())) == 0
            edgelist = list(np.array(list(nx.get_edge_attributes(self._G, 'euclidean_dist').keys()))[nodeindices])

    def _merge_roots_by_distance(self, dist):
        """
        Merge nodes into one single soma node that are 'dist' microns far away from the soma. Merging is done inplace.
        :param dist: double, distance in microns.
        """

        R = self._G
        nodes_to_remove = True
        root_ix = self.get_root()
        root = R.node[root_ix]

        # merge everything that is dist microns away from root node pos
        # dist for BC cells is < 50 voxel = 2.5 microns
        while (nodes_to_remove):
            nodes_to_remove = False
            for succ in R.successors(root_ix):
                edge = R.get_edge_data(root_ix, succ)
                if edge['euclidean_dist'] <= dist:
                    nodes_to_remove = True
                    for e in R.successors(succ):
                        n2 = R.node[int(e)]
                        d = np.sqrt(np.sum((root['pos'] - n2['pos']) ** 2))
                        R.add_edge(root_ix, e, euclidean_dist=d, path_length=d)
                    R.remove_node(succ)

        self._G = R

    def _merge_edges_on_path_by_displacement(self, start=1, disp=5):
        """
        Reduces the number of nodes along all neurites based on their displacement from the line connecting the outer
        branch points. It deletes all nodes B on path A-->B-->C that are maximally 'disp' microns displaced from the
        edge A-->C. The reduction is done inplace, to preserve the original structure please copy the NeuronTree before.
        :param start: node id, default = 1. Denotes the starting point within the NeuronTree.
        :param disp: double, displacement threshold in microns. It determines the perpendicular distance within which
        all nodes get deleted.
        """

        Tree = self._G
        current_node = start
        successors = Tree.successors(start)

        while (successors):

            pred = Tree.predecessors(current_node)
            succ = successors.pop(0)
            deg = Tree.out_degree()[current_node]
            if deg == 1:
                if pred and pred == Tree.predecessors(Tree.predecessors(succ)[0]):
                    p = Tree.node[pred[0]]
                    s = Tree.node[succ]
                    c = Tree.node[current_node]
                    d = np.linalg.norm(np.cross(c['pos'] - p['pos'], p['pos'] - s['pos'])) / np.linalg.norm(
                        c['pos'] - p['pos'])
                    if d < disp:
                        for s in Tree.successors(current_node):
                            path = nx.shortest_path_length(Tree, pred[0], s, weight='path_length')
                            d = np.sqrt(np.sum((p['pos'] - s['pos']) ** 2))
                            Tree.add_edge(pred[0], s, euclidean_dist=d, path_length=path)

                        Tree.remove_node(current_node)
            S = Tree.successors(succ)
            S[len(S):] = successors
            successors = S

            current_node = succ
        self._G = Tree

    def _merge_edges_on_path_by_edge_length(self, start=1, e=0.01):
        """
         It removes nodes B on path A-->B-->C that do not change the length of a direct edge A-->C by the amount of
         epsilon 'e' (in microns). The reduction is done inplace, to preserve the original structure please
         copy the NeuronTree before.
         changed for use with networkx v2 (works also in old version: edge -> adj)
        :param start: node id, default = 1. Denotes the starting point within the NeuronTree.
        :param e: double, tolerated epsilon (given in microns) of path length that we do not care to remove.
        """
        Tree = self._G
        current_node = start
        successors = list(Tree.successors(start))

        while (successors):

            pred = list(Tree.predecessors(current_node))
            succ = successors.pop(0)
            deg = Tree.out_degree()[current_node]

            if deg == 1:
                if pred and pred == Tree.predecessors(list(Tree.predecessors(succ))[0]):
                    p = Tree.node[pred[0]]
                    s = Tree.node[succ]

                    d = np.sqrt(np.sum((s['pos'] - p['pos']) ** 2))

                    d1 = Tree.adj[pred[0]][current_node]['euclidean_dist']
                    d2 = Tree.adj[current_node][succ]['euclidean_dist']
                    dist = d1 + d2

                    if dist - d < e:
                        for s in Tree.successors(current_node):
                            path = nx.shortest_path_length(Tree, pred[0], s, weight='path_length')
                            Tree.add_edge(pred[0], s, euclidean_dist=d, path_length=path)

                        Tree.remove_node(current_node)
            S = list(Tree.successors(succ))
            S[len(S):] = successors
            successors = S

            current_node = succ
        self._G = Tree

    def _get_branch_type(self, B):
        """
        get the type of a branch based on majority vote.
        :param B: subgraph, a branch within the NeuronTree
        :return: int, type id (1: soma, 2: axon, 3: basal dendrite, 4: apical dendrite). Type of the branch 'B'.
        """
        Y = nx.get_node_attributes(self._G, 'type')
        bt = []
        bt += [Y[k] for k in B]
        return np.round(np.mean(bt))

    def _unify_type(self):
        """
        Fixing routine: Set all branches originating from the soma to their majority vote type.
        """
        for s in self._G.successors(self.get_root()):
            S = nx.dfs_tree(self._G, s).nodes()
            t = self._get_branch_type(S)
            for k in S:
                self._G.node[k]['type'] = t

    def _clean_axon(self):
        """
        Fixing routine: Only keep one axon originating from the soma. If there is multiple present only the longest
        neurite is kept ( in terms of number of edges) the other ones get deleted. The manipulation happens inplace.
        """
        # clean up axon: only longest axon is kept
        axon_edges = self.get_axon_edges()

        if axon_edges:
            axon_edges = np.array(axon_edges)
            edges = axon_edges[(axon_edges[:, 0] == 1)]
            l = []

            for n in edges:
                l.append(len(nx.dfs_tree(self._G, n[1]).nodes()))

            m = max(l)
            ax = edges[l.index(m)][1]
            to_remove = edges[(edges[:, 1] != ax)]

            for e in to_remove:
                self._G.remove_nodes_from(nx.dfs_tree(self._G, e[1]).nodes())

    def make_pos_relative(self):
        """
        Deprecated! Makes all node positions relative to the soma. The soma will then have position (0,0,0).
        """
        root = self.get_root()
        if 'pos' in self.get_node_attributes():
            root_pos = self._G.node[root]['pos']

            for v in self._G.node:
                v_pos = self._G.node[v]['pos']
                self._G.node[v]['pos'] = v_pos - root_pos
        else:
            raise Warning('There is no position data assigned to the nodes.')

    def _make_tree(self):
        """
        Forces the networkx.DiGraph to be a tree. This routine was needed because the networkx.dfs_tree() function
        looses the original node and edge attributes. Manipulation done inplace.
        changed for use with networkx v2 (works also in old version: parameters with names)
        """

        G = self._G
        r = self.get_root()
        T = nx.dfs_tree(G, r)

        for n_attr in self.get_node_attributes():
            attr = nx.get_node_attributes(G, name=n_attr)
            nx.set_node_attributes(T, name=n_attr, values=attr)

        for e_attr in self.get_edge_attributes():
            attr = nx.get_edge_attributes(G, name=e_attr)
            nx.set_edge_attributes(T, name=e_attr, values=attr)

        self._G = T

    def get_graph(self):
        return self._G

    def get_mst(self):
        """
        Returns the minimal spanning tree of the Neuron.
        changed for use with networkx v2 (works also in old version: edge -> adj)
        :return:
            NeuronTree: mst. The minimal spanning tree representation of the original neuron.
        """
        # get the included nodes, which are soma, branch points and tips
        other_points = np.unique(np.append(self.get_branchpoints(), self.get_root()))
        tips = self.get_tips()

        # get the node data
        node_data = self.get_graph().node
        node_data_new = [(node, node_data[node]) for node in np.append(other_points, tips)]

        # get parent of each node and create edge_data
        nodes = set(tips)
        edge_data = self.get_graph().adj
        edge_data_new = []
        while nodes:
            current_node = nodes.pop()
            # if node is  not soma
            if current_node != self.get_root():
                cn = copy.copy(current_node)
                pred = list(nx.DiGraph.predecessors(self.get_graph(), current_node))[0]
                path_length = edge_data[pred][cn]['path_length']
                while pred not in other_points:
                    cn = pred
                    pred = list(nx.DiGraph.predecessors(self.get_graph(), pred))[0]
                    path_length += edge_data[pred][cn]['path_length']

                ec = np.sqrt(np.sum((node_data[current_node]['pos'] - node_data[pred]['pos']) ** 2))
                # add edge to edge_data_new
                edge_data_new.append((pred, current_node, dict(euclidean_dist=ec, path_length=path_length)))

                nodes.add(pred)  # adds the predecessor only once since nodes is a set
        return NeuronTree(node_data=node_data_new, edge_data=edge_data_new, nxversion=self._nxversion)

    def smooth_neurites(self, dim=1, window_size=21):

        from scipy.signal import savgol_filter

        G = copy.copy(self.get_graph())

        positions = nx.get_node_attributes(G, 'pos')

        smoothed = dict(zip(G.nodes(), [False] * len(G.nodes())))
        smoothed[1] = True

        dist_ = nx.single_source_dijkstra_path_length(G, source=1, weight='path_length')
        tips = self.get_tips()

        # get tips sorted by path length in descending order
        to_visit = tips[np.argsort([dist_[t] for t in tips])[::-1]]

        # smoothing
        for c in to_visit:
            predecessors = [c]
            for p in predecessors:
                # add all unsmoothed predecessors
                add = [p_ for p_ in G.predecessors(p) if not smoothed[p_]]
                predecessors += add

            path = [positions[p] for p in predecessors]
            if len(path) > window_size:
                data = np.array(path)
                datahat = savgol_filter(data[:, dim], window_size, 3, mode='nearest').reshape(-1,1)

                if dim == 0:
                    datahat = np.hstack([datahat.reshape(-1, 1), data[:, [1,2]]])
                elif dim == 1:
                    datahat = np.hstack([data[:, 0].reshape(-1, 1),
                                         datahat.reshape(-1, 1),
                                         data[:, 2].reshape(-1, 1)])
                elif dim == 2:
                    datahat = np.hstack([data[:, [0,1]],datahat.reshape(-1, 1)])
                    
                positions.update(dict(zip(predecessors, datahat)))

            smoothed.update(dict(zip(predecessors, [True] * len(predecessors))))

        # Create a new tree
        e_attr = []
        for e in G.edges():
            d = np.sqrt(np.sum((positions[e[0]] - positions[e[1]]) ** 2))
            e_attr.append(d)

        nx.set_node_attributes(G, 'pos', positions)
        nx.set_edge_attributes(G, 'euclidean_dist', dict(zip(G.edges(), e_attr)))
        nx.set_edge_attributes(G, 'path_length', dict(zip(G.edges(), e_attr)))

        S = NeuronTree(node_data=G.nodes(data=True), edge_data=G.edges(data=True), post_process=False, nxversion=self._nxversion)
        return S

    def get_node_attributes(self):
        """ returns the list of attributes assigned to each node.
            If no attributes are assigned it returns an empty list.
            changed for use with networkx v2 (works in all versions)
            :return:
                list : attribute names assigned to nodes
        """
        attr = []
        if self._G.nodes():
            if self._nxversion == 2:
                # changed for version 2.2 of networkX
                node_id = list(self._G.nodes().keys())[0]
            else:
                # this works only with version 1 of networkX
                node_id = self._G.nodes()[0]
            attr = list(self._G.node[node_id].keys())
        return attr

    def get_edge_attributes(self):
        """
            Returns the list of attributes assigned to each edge.
            If no attributes are assigned it returns an empty list.
            changed for use with networkx v2 (works in all versions)
        :return:
            list : attribute names assigned to edges
        """
        attr = []
        if self._G.edges():
            if self._nxversion == 2:
                # changed for version 2.2 of networkX
                e = list(self._G.edges().keys())[0]
            else:
                # this works only with version 1 of networkX
                e = self._G.edges()[0]
            attr = list(self._G.adj[e[0]][e[1]].keys())
        return attr

    def get_root(self):
        """
        Returns the root of the Neuron which is typically the soma with node id = 1.
        :return: int, node id of the soma.
        """
        try:
            root = np.min(self.nodes(type_ix=1))
        except (ValueError, KeyError):
            print('No node is attributed as being the soma. Returning the smallest node id.')
            root = np.min(self.nodes())

        return root

    def reduce(self, method='mst', e=0.01):
        """
        Reduces the number of nodes in the given tree by pruning nodes on paths like A-->B-->C. B is pruned and the
        resulting edge A-->C is inserted under circumstances defined by the keyword 'method'.
        :param method: (default 'mst') Defines the method by which the number of nodes is reduced. Possible methods are
                'mst' -- deletes all nodes in between branch points. Results in the minimal spanning tree of the neuron
                representation.
                'dist' -- deletes all nodes B on path A-->B-->C that do not change the length of edge A-->C by the
                amount of epsilon e
                (in microns).
                'disp' -- deletes all nodes B on path A-->B-->C that are maximally e microns displaced from the edge
                A-->C.
        :param e: float (default 0.01) error margin for methods 'dist' and 'disp'. e is interpreted in microns.
        :return: None. The tree is pruned in place. If the original tree is desired to be conserved then copy the tree
        data beforehand via the copy/deepcopy constructor.
        """

        if type(self._G) == nx.classes.graph.Graph:
            self._make_tree()

        if method == 'mst':
            mst = self.get_mst()
            self._G = mst.get_graph()
        elif method == 'dist':
            self._merge_edges_on_path_by_edge_length(e=e)
        elif method == 'disp':
            self._merge_edges_on_path_by_displacement(disp=e)
        else:
            raise NotImplementedError('Method {0} is not implemented!'.format(method))

    def truncate_nodes(self, perc=.1, no_trunc_nodes=None):
        """
        Truncates the number of nodes by the given percentage. The nodes are pruned equally from the tips inward until
        the approximated percentage has been deleted or at least the given number of nodes is truncated.
        :param perc: float. Percent of nodes that are to be truncated between 0 and 1. Default=0.1
        :param no_trunc_nodes: int . Optional. Number of nodes to be truncated
        :return: The truncated tree.
        """
        T = copy.copy(self)
        nodes = self.nodes()
        total_no_nodes = len(nodes)

        if not no_trunc_nodes:
            no_trunc_nodes = int(total_no_nodes * perc)

        tips = T.get_tips()

        for t in tips:
            nodes.remove(t)

        while total_no_nodes - no_trunc_nodes < len(nodes):
            T = NeuronTree(graph=T.get_graph().subgraph(nodes), nxversion=self._nxversion)
            nodes = T.nodes()

            tips = T.get_tips()

            for t in tips:
                nodes.remove(t)
        
        return T

    def nodes(self, type_ix=None, data=False):
        """
        Wrapper function to the networkx.nodes() function. It returns a list of all nodes.
        :param type_ix: int, optional, default = None. Determines the type of nodes to be returned.
        Options are None (= all nodes), 2 (= axonal nodes) and 3 (= dendritic nodes).
        :param data: boolean, default=False. If set to True the node attribute data is returned as well.
        :return: list of nodes
        """

        if type_ix is None:
            nodes = self._G.nodes(data=data)
        else:
            if type(type_ix) == list:
                nodes = [k for k in self._G.node if self._G.node[k]['type'] in type_ix]
            else:
                nodes = [k for k in self._G.node if self._G.node[k]['type'] == type_ix]
        return nodes

    def edges(self, start=None, type_ix=None, data=False):
        """
        Wrapper function to the networkx.edges() function. It returns a list of edges.
        :param start: int, node id, determines the starting edge within the neuron
        :param type_ix: int, optional, default = None. Determines the type of edges to be returned.
        Options are None (= all edges), 2 (= axonal edges) and 3 (= dendritic edges).
        :param data: boolean, default=False. If set to True the edge attribute data is returned as well.
        :return: list of edges
        """
        if type_ix is None:
            edges = self._G.edges(start, data=data)
        else:
            nodes = self.nodes(type_ix=type_ix)
            edges = [x for x in self._G.edges(start, data=data) if (x[0] in nodes or x[1] in nodes)]
        return edges

    def get_dendrite_nodes(self, data=False):
        """
        Returns all dendritic nodes.
        :param data: boolean, default=False. If set to True the node attribute data is returned as well.
        :return: list of dendritic nodes
        """
        return np.array(self.nodes(type_ix=[3, 4] , data=data))

    def get_axon_nodes(self, data=False):
        """
        Returns all axonal nodes.
        :param data: boolean, default=False. If set to True the node attribute data is returned as well.
        :return: list of axonal nodes
        """
        return np.array(self.nodes(type_ix=2, data=data))

    def get_axon_edges(self, start=None, data=False):
        axon_nodes = self.get_axon_nodes()
        return [x for x in self._G.edges(start,data=data) if (x[0] in axon_nodes or x[1] in axon_nodes)]

    def get_dendrite_edges(self, start=None, data=False):
        dendrite_nodes = self.get_dendrite_nodes()
        return [x for x in self._G.edges(start,data=data) if (x[0] in dendrite_nodes or x[1] in dendrite_nodes)]

    def get_tips(self):
        """
        Returns a list of tips, so the ending nodes within the neuron.
        changed for use with networkx v2 (works also in old version: edge -> adj)
        :return: list of node ids of tips.
        """
        E = self._G.adj
        return np.array([e for e in E if E[e] == {}])

    def get_branchpoints(self):
        """
        Returns a list of branch point ids of the neuron.
        :return: list of node ids of all branch points (nodes that have more than one successor).
        """
        # is this correct???:
        bp_indx = np.where(np.array(np.sum(nx.adjacency_matrix(self._G).toarray(), axis=1)).flatten() > 1)
        #bp_indx = np.where(np.array(np.sum(nx.adjacency_matrix(self._G), axis=1)).flatten() > 1)
        return np.array(self.nodes())[bp_indx]

    def get_dendritic_tree(self):
        """
        Returns the dendrites as a new NeuronTree.
        :return: NeuronTree
        """
        nodes = list(self.get_dendrite_nodes())
        nodes.insert(0, self.get_root())
        subgraph = nx.subgraph(self._G, nodes)
        return NeuronTree(graph=subgraph, nxversion=self._nxversion)

    def get_axonal_tree(self):
        """
        Returns the axon as a new NeuronTree
        :return: NeuronTree
        """
        nodes = list(self.get_axon_nodes())
        nodes.insert(0, self.get_root())
        subgraph = nx.subgraph(self._G, nodes)
        return NeuronTree(graph=subgraph, nxversion=self._nxversion)

    def get_adjacency_matrix(self, weight=None):
        """
        Returns the adjacency matrix of the Tree saved in self._G. weight can be None, 'euclidean_dist' or 'path_length'
        :param weight: edge attribute that is considered for the adjacency matrix A. Default = None, then A only
        contains the structural connectivity between nodes. If set to 'euclidean_dist' or 'path_length' the adjacency
        matrix is weighted accordingly.
        :return: sparse array
        """
        return nx.adjancency_matrix(self._G, weight=weight)

    def get_extend(self):
        """
        Returns the maximal extend in x, y and z direction.
        :return: 1x3 numpy.array
        """
        P = np.array(list(nx.get_node_attributes(self._G, 'pos').values()))
        return np.max(P, axis=0) - np.min(P, axis=0)

    def get_root_angle_dist(self, angle_type='axis_angle', **kwargs):
        """
               Returns the histogram over the root angle distribution over a tree. Root angle denotes the orientation of
               each edge with respect to the root (soma).
               :param self: NeuronTree object
               :param bins: int
                    Number of bins used in the histogram. Default is 10.
               :param angle_type: either 'axis_angle' or 'euler'
                    Defines the type of angle that is calculated. Euler angles are defined as angles around the canonical
                    euler axes (x, y and z). Axis angles are defined with respect to the rotation axis between two
                    vectors.
               :returns:
                    hist: histogram over angles
                    edges: edges of histogram. For further information see numpy.histogramdd()
           """

        angles = []
        if angle_type == 'axis_angle':
            dim = 1
            func = lambda u,v: angle_between(u, v)
        elif angle_type == 'euler':
            dim =3
            func = lambda u,v: rotationMatrixToEulerAngles(get_rotation_matrix(u, v))
        else:
            raise NotImplementedError('Angle type %s is not implemented' % angle_type)

        for n1, n2 in self._G.edges():
            u = self._G.node[n2]['pos'] - self._G.node[n1]['pos']
            v = self._G.node[n1]['pos'] - self._G.node[self.get_root()]['pos']
            angles.append(func(u, v))
        angles = np.array(angles)
        hist = np.histogramdd(angles, range=[[0, np.pi]]*dim, **kwargs)
        return hist

    def get_branch_order(self):
        """
        Returns the dictionary of the branch order of each node. The branch order denotes the number of branch points
        that are crossed when tracing the path back to the soma.
        :return:
            d: dict
            Dictionary of the form {u: branch_order} for each node.
        """
        return self._get_branch_order(1,0)

    def _get_branch_order(self, start, bo):
        """
        Returns the dictionary assigning the right branch order to each node.
        :param start: starting node
        :param bo: starting branch order
        :return:
            d: dict
            Dictionary of the form {u: branch_order} for each node reachable from starting node.
        """
        d = {}
        edges = self.edges(start)
        d[start] = bo
        if len(edges) > 1:
            for e in edges:
                d.update(self._get_branch_order(e[1], bo + 1))
        elif len(edges) == 1:
            d.update(self._get_branch_order(edges[0][1], bo))
        return d

    def _get_distance(self, dist='path_from_soma', weight='euclidean_dist', as_dict=False):
        """
        Returns the distance
        :param dist: String, defines the distance measure to be used (default is 'path_from_soma'), Options are
        'path_from_soma' and 'branch_order'.
        :param weight: String. Defines the edge attribute that is considered for the distance. Options are
        ['euclidean_dist' = default, 'path_lenght']
        :param as_dict: boolean, default = False. Determines whether the distance are returned as a dictionary of the
        form {'node_id': ,distance} or as an numpy.array.
        :return: Dictionary or numpy.array of the defined distance measure from each node to the soma.
        """
        if dist == 'path_from_soma':

            if as_dict:
                dist_ = nx.single_source_dijkstra_path_length(self.get_graph(), source=self.get_root(),
                                                              weight=weight)
            else:
                dist_ = np.array(list(nx.single_source_dijkstra_path_length(self.get_graph(), source=self.get_root(),
                                                                            weight=weight).values()))
        elif dist == 'branch_order':
            dist_ = self._get_branch_order(1, 0)
            if not as_dict:
                dist_ = np.array(list(dist_.values()))

        else:
            raise NotImplementedError

        return dist_

    def get_segment_length(self, dist='path_length'):
        """
        Returns the dictionary of the segment length in microns of each branch where dist denotes the distance measure.
         Possible options are 'path_length' an 'euclidean_dist'. The keys of the dictionary denote the
        tuples of the starting and end node of each segment.
        :param dist: String, options ['path_length', 'euclidean_dist']
        :return:
            d: dict
            Dictionary of the form {(n_s, n_e): segment length[u] in either 'path length' or 'euclidean distance' }
        """

        T = self.get_mst()

        segment_length = nx.get_edge_attributes(T.get_graph(), dist)
        return segment_length

    def get_kde_distribution(self, key, dist=None):
        """
        Returns a Gaussian kernel density estimate of the distribution of the measure defined by _key_.
        changed for use with networkx v2 (works also in old version: edge -> adj)
        :param key: String, measure to calculate the kernel density estimate from. Options are: 'thickness',
        'path_angle', 'branch_angle', 'branch_order' and 'distance'.
        :param dist: Optional. Allows to express the measure as a function of distance. If set to None, the kde is
        one-dimensional, otherwise it is has two dimensions. Options for available distance measures,
        see _get_distance().
        :return: A one- or two-dimensional Gaussian kernel density estimate of measure _key_.
        """

        if key == 'thickness':
            thickness = np.array(list(nx.get_node_attributes(self.get_graph(), 'radius').values()))

            if dist:
                data = np.array(list(zip(self._get_distance(dist), thickness)))
                return stats.gaussian_kde(data)

            else:
                return stats.gaussian_kde(thickness)
        elif key == 'path_angle':

            successors = nx.dfs_successors(self.get_graph(), 1)
            path_angle = []
            nodes = []
            for n1, n2 in self.edges():
                u = self.get_graph().node[n2]['pos'] - self.get_graph().node[n1]['pos']
                try:
                    for succ in successors[n2]:
                        v = self.get_graph().node[succ]['pos'] - self.get_graph().node[n2]['pos']
                        path_angle.append(angle_between(u, v) * 180 / np.pi)  # convert angles into degree
                        nodes.append(n2)
                except KeyError:
                    continue

            if dist:
                distance = self._get_distance(dist, as_dict=True)
                distance = [distance[n] for n in nodes]

                data = np.array(list(zip(distance, path_angle)))
                return stats.gaussian_kde(data)
            else:
                return stats.gaussian_kde(path_angle)

        elif key == 'branch_angle':
            branch_angles = []
            branchpoints = self.get_branchpoints()
            for bp in branchpoints:
                successors = list(self.get_graph().adj[bp].keys())
                branches = []
                for succ in successors:  # create a vector for each branching edge
                    v = self.get_graph().node[succ]['pos'] - self.get_graph().node[bp]['pos']
                    branches.append(v)
                for u, v in combinations(branches, 2):
                    branch_angles.append(angle_between(u, v) * 180 / np.pi)

            if dist:
                distance = self._get_distance(dist, as_dict=True)
                distance = [distance[bp] for bp in branchpoints]
                data = np.array(list(zip(distance, branch_angles)))
                return stats.gaussian_kde(data)
            else:
                return stats.gaussian_kde(branch_angles)

        elif key == 'branch_order':
            branch_order = list(self._get_branch_order(1, 0).values())

            if dist:
                data = np.array(list(zip(self._get_distance(dist), branch_order)))
                return stats.gaussian_kde(data)
            else:
                return stats.gaussian_kde(branch_order)

        elif key == 'distance':
            distance = np.array(list(nx.single_source_dijkstra_path_length(self.get_graph(), source=self.get_root(),
                                                                           weight='euclidean_dist').values()))

            if dist:
                data = np.array(list(zip(self._get_distance(dist), distance)))
                return stats.gaussian_kde(data)
            else:
                return stats.gaussian_kde(distance)
        else:
            raise NotImplementedError

    def get_distance_dist(self, **kwargs):
        """
            Returns histogram of the distance distribution from soma
        :param kwargs: optional parameters passed to histogram calculation (see numpy.histogramdd)
        :return:
            hist    1D histogram
            edges   bin edges of the histogram in hist. Definition as in numpy.histogrammdd
        """

        dist = np.array(list(nx.single_source_dijkstra_path_length(self.get_graph(), source=self.get_root(),
                                                                   weight='euclidean_dist').values()))

        return np.histogram(dist, **kwargs)

    def get_branch_order_dist(self, **kwargs):
        """
            Returns histogram of the branch order distribution
        :param kwargs: optional parameters passed to histogram calculation (see numpy.histogramdd)
        :return:
            hist    1D histogram
            edges   bin edges of the histogram in hist. Definition as in numpy.histogrammdd
        """
        data = list(self._get_branch_order(1, 0).values())

        # TODO make it assignable a dist_measure
        return np.histogram(data, **kwargs)

    def get_thickness_dist(self, dist_measure=None, **kwargs):
        """
            Returns the distribution of the neurons thickness in microns in form of a histogram.
            The distribution can be calculated against a distance measure, namely 'path_from_soma' or 'branch_order'.
        :param dist_measure:
        :param kwargs:
        :return:
            hist    1D or 2D histogram
            edges   bin edges used
        """

        # get the thickness of each node
        thickness_dict = nx.get_node_attributes(self.get_graph(), 'radius')

        # delete the soma since it usually skews the distribution
        thickness_dict.pop(self.get_root())

        thickness = np.array(list(thickness_dict.values()))

        if dist_measure:

            dist = self._get_distance(dist_measure)
            data = np.array(list(zip(dist, thickness)))
            return np.histogramdd(data, **kwargs)
        else:

            return np.histogram(thickness, **kwargs)

    def get_path_angles(self):
        """
        Returns a dictionary of path angles between two edges. Angles are reported in degree
        :return: dict of path angles btw two edges.
                d[u][v][w] returns the angle between edge (u,v) and (v,w)
        """
        # get the depth first search successors from the soma (id=1).
        successors = nx.dfs_successors(self.get_graph(), 1)
        path_angle = {}

        for u, v in self.edges():
            e1 = self.get_graph().node[v]['pos'] - self.get_graph().node[u]['pos']
            path_angle[u] = {}
            try:
                path_angle[u][v] = {}
                for w in successors[v]:
                    e2 = self.get_graph().node[w]['pos'] - self.get_graph().node[v]['pos']
                    path_angle[u][v][w] = angle_between(e1, e2) * 180 / np.pi
            except KeyError:
                continue

        return path_angle

    def get_path_angle_dist(self, dist_measure=None, **kwargs):
        """
            Returns the distribution of the path angle, so the angles that are made between two consecutive segments.
            The distribution can be calculated against a distance measure, namely 'path_from_soma' or 'branch_order'.
            The path angles are returned in degree.
        :param dist_measure: string. default = None
            Defines the distance measure against whom the thickness is calculated. Possible choices are 'path_from_soma'
            or 'branch_order'.
        :param kwargs: additional arguments to be passed to histogram calculation
        :return:
            hist    1- or 2D histogram
            edges   bin edges of the histogram in hist. Definition as in numpy.histogrammdd
        """

        successors = nx.dfs_successors(self.get_graph(), 1)
        path_angle = []
        nodes = []
        for n1, n2 in self.edges():
            u = self.get_graph().node[n2]['pos'] - self.get_graph().node[n1]['pos']
            try:
                for succ in successors[n2]:
                    v = self.get_graph().node[succ]['pos'] - self.get_graph().node[n2]['pos']
                    path_angle.append(angle_between(u, v) * 180 / np.pi)  # convert angles into degree
                    nodes.append(n2)
            except KeyError:
                continue

        if dist_measure:
            distances = self._get_distance(dist_measure, as_dict=True)
            dist = [distances[n] for n in nodes]
            data = np.array(list(zip(dist, path_angle)))
            return np.histogramdd(data, **kwargs)
        else:
            return np.histogram(path_angle, **kwargs)

    def get_branch_angles(self):
        """
        Returns the list of branch angles.
        changed for use with networkx v2 (works also in old version: edge -> adj)
        :return:
            branch_angles   list of branch angles.
        """

        from itertools import combinations

        branch_angles = []
        branchpoints = self.get_branchpoints()
        for bp in branchpoints:
            successors = list(self.get_graph().adj[bp].keys())
            branches = []
            for succ in successors:  # create a vector for each branching edge
                v = self.get_graph().node[succ]['pos'] - self.get_graph().node[bp]['pos']
                branches.append(v)
            for u, v in combinations(branches, 2):
                branch_angles.append(angle_between(u, v) * 180 / np.pi)
        return branch_angles

    def get_branch_angle_dist(self, dist_measure=None, **kwargs):
        """
            Returns the distribution of the branch angles, so the angles that are made between two branching segments.
            The distribution is calculated against a distance measure, namely 'path_from_soma' or 'branch_order'.
            The branch angles are returned in degree.
        :param dist_measure: string. default = None
            Defines the distance measure against whom the thickness is calculated. Possible choices are 'path_from_soma'
            or 'branch_order'.
        :param kwargs: dditional arguments to be passed to histogram calculation
        :return:
            hist    2D histogram
            edges   bin edges of the histogram in hist. Definition as in numpy.histogrammdd
        """

        branch_angles = self.get_branch_angles()
        branchpoints = self.get_branchpoints()
        if dist_measure:
            distances = self._get_distance(dist_measure, as_dict=True)
            dist = [distances[n] for n in branchpoints]
            data = np.array(list(zip(dist, branch_angles)))
            return np.histogramdd(data, **kwargs)
        else:
            return np.histogramdd(branch_angles, **kwargs)

    def get_segment_length_dist(self, segment_dist='path_length', branch_order=False, **kwargs):
        """
        Returns the histogram over segment lengths within the neuron. A segment is the part of a neurite between two
        branch points.
        :param segment_dist:    String, default='path_length',
                                option=['path_length', 'euclidean_dist'] determines if the segment length is measured as
                                euclidean distance between the end points or as the segment's path length.
        :param branch_order:   bool, default=False, Defines whether segment path length is calculated as a function of
                                branch order.
        :param kwargs: parameters that can be passed to numpy.histogram
        :return:
            hist    histogram over segment lengths
            edges   bins used for the histogram _hist_
        """

        segment_length = self.get_segment_length(segment_dist)
        if branch_order:
            bo = self.get_branch_order()
            data = np.array([[bo[item[0][1]], item[1]] for item in segment_length.items()])
            return np.histogramdd(data, **kwargs)
        else:
            return np.histogram(list(segment_length.values()), **kwargs)

    def get_volume(self):
        """
        Returns the volume for each segment in the tree.
        :return: dictionary of the form d[(u,v)] = volume(u,v)
        """

        d = {}
        for e in self.edges(data=True):

            h = e[2]['euclidean_dist']
            r = self.get_graph().node[e[0]]['radius']
            R = self.get_graph().node[e[1]]['radius']

            d[(e[0], e[1])] = (1/3)*np.pi*h*(r*r + r*R + R*R)
        return d

    def get_surface(self):
        """
        Returns the surface for each segment in the tree treating each edge as a pipe (without closing lids!).
        :return: dictionary of the form d[(u,v)] = surface(u,v)
        """
        d = {}
        for e in self.edges(data=True):
            h = e[2]['euclidean_dist']
            r = self.get_graph().node[e[0]]['radius']
            R = self.get_graph().node[e[1]]['radius']

            d[(e[0], e[1])] = np.pi*(r + R)*np.sqrt((R-r)**2 + h**2)
        return d

    def get_sholl_intersection_profile(self, proj='xy', steps=36, centroid='centroid'):
        """
        Calculates the Sholl intersection profile of the neurons projection determined by the parameter _proj_. The
        Sholl intersection profile counts the intersection of the neurites with concentric circles with increasing
        radii. The origin of the concentric circles can be chosen to either be the centroid of the
        projected neuron's convex hull or to be the soma.
        :param proj: 2D projection of all neurites. Options are 'xy', 'xz' and 'yz'
        :param steps: number of concentric circles centered around _centroid_. Their radii are determined as the
        respective fraction of the distance from the centroid to the farthest point of the convex hull.
        :param centroid: Determines the origin of the concentric circles. Options are 'centroid' or 'soma'.
        :return:
            intersections  list, len(steps), count of intersections with each circle.
            intervals list, len(steps +1), radii of the concentric circles.
        """

        G = self.get_graph()
        coordinates = []

        if proj == 'xy':
            indx = [0, 1]
        elif proj == 'xz':
            indx = [0, 2]
        elif proj == 'yz':
            indx = [1, 2]
        else:
            raise ValueError("Projection %s not implemented" % proj)

        # get the coordinates of the points
        for e in G.edges():
            p1 = np.round(G.node[e[0]]['pos'], 2)
            p2 = np.round(G.node[e[1]]['pos'], 2)
            coordinates.append((p1[indx], p2[indx]))

        # remove illegal points
        coords = [c for c in coordinates if (c[0][0] != c[1][0] or c[0][1] != c[1][1])]

        lines = MultiLineString(coords).buffer(0.0001)
        bounds = np.array(lines.bounds).reshape(2, 2).T
        if centroid == 'centroid':
            center = np.array(lines.convex_hull.centroid.coords[0])
            p_circle = Point(center)
        elif centroid == 'soma':
            center = G.node[self.get_root()]['pos'][indx]
            p_circle = Point(center)
        else:
            raise ValueError("Centroid %s is not defined" % centroid)

        # get the maximal absolute coordinates in x and y
        idx = np.argmax(np.abs(bounds), axis=1)
        r_max = np.linalg.norm(center - bounds[idx, [0, 1]])

        intersections = []
        intervals = [0]
        for k in range(1, steps + 1):

            r = (r_max / steps) * k
            c = p_circle.buffer(r).boundary

            i = c.intersection(lines)
            if type(i) in [Point, LineString]:
                intersections.append(1)
            else:
                intersections.append(len(i))
            intervals.append(r)

        return intersections, intervals

    def get_psad(self):
        """
        Returns the proportional sum of absolute deviations (PSAD) for each branching node in the neuron. The PSAD is
        a measure of topological tree asymmetry and is defined for a branching node p as
        PSAD_p = m/(2*(m-1)*(n-m) * sum_m |r_i - n/m|
        where m is the out_degree of node p, n is the number of leaves of the subtree starting at p and r_i is the
        number of leaves of the i_th subtree of p. For a more detailed definition see:
        Verwer, Ronald WH, and Jaap van Pelt.
        "Descriptive and comparative analysis of geometrical properties of neuronal tree structures."
        Journal of neuroscience methods 18.1-2 (1986): 179-206.

        :return:
            w:      dict boolean, dict of weights [0,1] indicating which subtrees have more than 3 leaves.
            psad:   dict, dictionary of the form {bp: PSAD} for each  in neuron.

        """

        def get_no_of_leaves(G, start, out_degree):
            """
            Helper-function: Returns the number of leaves that are attached to start point in graph G.
            downstream.
            :param G: nx.Digraph
            :param start: int   point to get the number of terminals from.
            :param out_degree:  dict . Dictionary of out degrees for each node in G.
            :return: int
            """
            t = 0
            if out_degree[start] == 0:
                t = 1
            else:
                for s in G.successors(start):
                    t += get_no_of_leaves(G, s, out_degree)
            return t

        G = self.get_graph()

        branchpoints = self.get_branchpoints()
        nodes = self.nodes()
        out_degree = nx.DiGraph.out_degree(G)

        degree = dict(zip(nodes, [0] * len(nodes)))
        for n in nodes:
            degree[n] = get_no_of_leaves(G, n, out_degree)

        psad = dict()
        for bp in branchpoints:
            m = out_degree[bp]
            n = degree[bp]

            PSAD = 0
            if m < n:  # only calculate for asymmetric trees
                for s in G.successors(bp):
                    PSAD += np.abs(degree[s] - n / m)
                PSAD = PSAD * m / (2 * (m - 1) * (n - m))
            psad[bp] = PSAD

        w = dict(zip(psad.keys(), [degree[n] > 3 for n in psad.keys()]))
        return w, psad

    def _resample_tree_data(self, dist=1):
        """
        Re-sample new nodes along the tree in equidistant distance dist (given in microns) and calculates the respective
        node and edge attributes to generate a new NeuronTree.
        All original branching points are kept!
        :param dist: distance (in microns) at which each neurite is resampled.
        :return:
            POS     list of positions
            PID     list of parent ids
            TYPE    list of types for each resampled node
            RADIUS  list of radii for each resampled node
            DIST
        """
        P = []
        PID = []
        TYPE = []
        RADIUS = []
        DIST = []

        pid = -1

        for (u, v, edata) in self._G.edges(data=True):
            n1 = self._G.node[u]
            n2 = self._G.node[v]
            e = edata['euclidean_dist']
            v = n2['pos'] - n1['pos']
            v /= np.linalg.norm(v)

            # line between two nodes
            g = lambda x: x * v + n1['pos']

            # radius function

            if n1['radius'] == n2['radius']:
                r = lambda x: n2['radius']
            else:
                x = [0, e]
                y = [n1['radius'], n2['radius']]
                r = interp1d(x, y)

            n = list(n1['pos'])

            if n in P:
                pid = P.index(n)
            else:
                P.append(n)
                PID.append(pid)
                pid += 1
                TYPE.append(n1['type'])
                DIST.append(0)
                RADIUS.append(n1['radius'])

            for a in range(1, int(e / dist)):
                P.append(list(g(a * dist)))
                PID.append(pid)
                pid = len(P) - 1
                TYPE.append(n2['type'])
                DIST.append(dist)
                RADIUS.append(np.array([r(a * dist)])[0])

            P.append(list(n2['pos']))
            PID.append(pid)
            pid += 1
            TYPE.append(n2['type'])
            DIST.append(e % dist)
            RADIUS.append(n2['radius'])

        P = np.array(P)

        return P, PID, TYPE, RADIUS, DIST

    def resample_tree(self, dist=1):
        """
         Re-sample new nodes along the tree in equidistant distance dist (given in microns) and return a new
         NeuronTree with the resampled data. Original branch points are kept.
        :param dist: distance (in microns) at which to resample
        :return: NeuronTree
        """

        (pos, pid, t, r, d) = self._resample_tree_data(dist)
        n_attr = [dict(pos=pos[i], type=t[i], radius=r[i]) for i in range(len(d))]

        nodes = list(zip(range(1, len(pid) + 1), n_attr))

        e_attr = [dict(euclidean_dist=d[i], path_length=d[i]) for i in range(len(d))]
        edges = list(zip(np.array(pid) + 1, range(1, len(pid) + 1), e_attr))

        T = NeuronTree(node_data=nodes, edge_data=edges[1:], nxversion=self._nxversion)

        return T

    def get_histogramdd(self, decomposition='ica', dim=3, proj_axes=None, whiten=True,
                        nbins=100, r=None, sampling_dist=0.01):
        p = NeuronTree.resample_nodes(self._G, sampling_dist)

        if decomposition:
            if decomposition == 'pca':
                # find principal axes of the 3D point cloud using PCA
                pca = PCA(n_components=dim, whiten=whiten)
                results = pca.fit_transform(p)
            elif decomposition == 'ica':
                ica = FastICA(n_components=dim, whiten=whiten)
                results = ica.fit_transform(p)
            else:
                raise ValueError('decomposition {0} is not implemented.'.format(decomposition))
        else:
            if proj_axes:
                results = p[:, proj_axes]
            else:
                results = p
        if r:
            return np.histogramdd(results, bins=(nbins,) * dim, range=r, normed=True)
        else:
            return np.histogramdd(results, bins=(nbins,) * dim, normed=True)


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
    ########################################## DRAWING FUNCTIONS #############################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

    def get_node_colors(self):
        """
        Returns a list of colors for each node as it appears in the graph. Used to give the 2D graph representation
        node colors.
        :return: list of color char tags for each node in graph. 'g' for axon, 'y' for dendrite and 'grey' for soma.
        """
        axon_nodes = self.get_axon_nodes()
        dendrite_nodes = self.get_dendrite_nodes()

        colors = []
        for node in self._G.nodes():
            if node in axon_nodes:
                colors.append('g')
            elif node in dendrite_nodes:
                colors.append('y')
            else:
                colors.append('grey')
        return colors

    def draw_3D(self, fig=None, ix=111, reverse=True, r_axis='z', axon_color='grey', dendrite_color='darkgrey'):
        """
        Draws a stick figure neuron in 3D.
        :param fig: figure in which to draw. If None a new figure is created.
        :param ix: index of the subplot within the figure 'fig'. Default: 111
        :param reverse: Default=True. Determines whether the axis specified in 'r_axis' is inverted.
        :param r_axis: Default = 'z'. Defines the axis that is inverted if 'reverse' is set to True. Possible values are
        'x', 'y' and 'z'.
        :param axon_color: Color, default='grey'. Defines the color of the axon.
        :param dendrite_color: Color, default='darkgrey'. Defines the color of the dendrites.
        """

        nodes = [k for k in nx.get_node_attributes(self._G, 'pos').values()]
        nodes = np.array(nodes)

        t = [axon_color if k == 2 else dendrite_color for k in nx.get_node_attributes(self._G, 'type').values()]

        # plot G
        if not fig:
            fig = plt.figure()
        ax = fig.add_subplot(ix, projection='3d')
        ax.scatter(nodes[:, 0], nodes[:, 1], nodes[:, 2], c=t, marker='.')
        plt.hold
        root_pos = self._G.node[self.get_root()]['pos']
        ax.scatter(root_pos[0], root_pos[1], root_pos[2], c='k', marker='^')

        colors = ['k', axon_color, dendrite_color]

        for k, e in enumerate(self._G.edges()):
            n1 = self._G.node[e[0]]
            n2 = self._G.node[e[1]]
            v = np.array([n1['pos'], n2['pos']])

            ax.plot3D(v[:, 0], v[:, 1], v[:, 2], c=colors[int(n2['type']) - 1])

        ax = plt.gca()
        if reverse:
            if 'z' in r_axis:
                ax.set_zlim(ax.get_zlim()[::-1])
            if 'x' in r_axis:
                ax.set_xlim(ax.get_xlim()[::-1])
            if 'y' in r_axis:
                ax.set_ylim(ax.get_ylim()[::-1])

        ax.set_xlabel('X [microns]')
        ax.set_ylabel('Y [microns]')
        ax.set_zlabel('Z [microns]')

    def draw_3D_volumetric(self, fig=None, ix=111, axon_color='grey', dendrite_color='darkgrey'):
        """
        Draws a volumetric neuron in 3D.
        changed for use with networkx v2 (works also in old version: edge -> adj)
        :param fig: Figure in which to draw. If None a new figure is created.
        :param ix: Index of the subplot within the figure 'fig'. Default: 111
        :param axon_color: Color, default='grey'. Defines the color of the axon.
        :param dendrite_color: Color, default='darkgrey'. Defines the color of the dendrites.
        """

        P = nx.get_node_attributes(self._G, 'pos')          # get 3D position for each node
        Rad = nx.get_node_attributes(self._G, 'radius')     # get radius for each node
        Type = nx.get_node_attributes(self._G, 'type')      # get type for each node

        # number of surfaces for each volume segment
        num = 5

        u = np.linspace(0, 2 * np.pi, num=num)
        v = np.linspace(0, np.pi, num=num)

        # parametric surface of a sphere
        fx = lambda r: r * np.outer(np.cos(u), np.sin(v))
        fy = lambda r: r * np.outer(np.sin(u), np.sin(v))
        fz = lambda r: r * np.outer(np.ones(np.size(u)), np.cos(v))

        # Create a figure
        if not fig:
            fig = plt.figure()

        # Add three dimensional axis
        ax = fig.add_subplot(ix, projection='3d')

        unit_z = np.array([0, 0, 1])
        
        # plot the nodes as spheres
        for i in self.nodes():
            pos = P[i]
            r = Rad[i]
            if Type[i] == 2:
                c = axon_color
            elif Type[i] == 1:  # plot the soma black
                c = 'k'
            else:
                c = dendrite_color
            ax.plot_surface(fx(r) + pos[0], fy(r) + pos[1], fz(r) + pos[2], color=c)

        # plot segments as cone frustums
        for e in self._G.edges():
            if e[0] == 1:
                a = Rad[e[1]]
                b = Rad[e[1]]
            else:
                a = Rad[e[0]]
                b = Rad[e[1]]

            if Type[e[1]] == 2:
                c = axon_color
            else:
                c = dendrite_color

            h = self._G.adj[e[0]][e[1]]['euclidean_dist'] # length of the edge
            # translation
            T = P[e[0]]
            # rotation
            R = get_rotation_matrix(unit_z, P[e[1]] - P[e[0]])  # rotate edge

            k = np.linspace(0, h, num)
            t = np.linspace(0, 2 * np.pi, num)

            # parametric surface of a cone frustrum
            cx = lambda k, t: np.outer((a * (h - k) + b * k) / h, np.cos(t))
            cy = lambda k, t: np.outer((a * (h - k) + b * k) / h, np.sin(t))

            F = np.array([cx(k, t), cy(k, t), np.meshgrid(k, k)[1]]).T

            R_ = np.reshape(np.tile(R.T, (1, num)).T, [num, 3, 3])
            E = np.einsum('lij, lkj->lki', R_, F)
            ax.plot_surface(E[:, :, 0] + T[0], E[:, :, 1] + T[1], E[:, :, 2] + T[2], color=c)

        ax.set_xlabel('X [microns]')
        ax.set_ylabel('Y [microns]')
        ax.set_zlabel('Z [microns]')

    def draw_tree(self, edge_labels=False, **kwds):
        """
        Draw neuron as a planar tree.
        :param edge_labels: Boolean, default=False. Determines if edgelabels are drawn as well.
        :param kwds: arguments that can be passed to the nx.draw_networkx function.
        """

        pos = nx.drawing.nx_agraph.graphviz_layout(self._G, prog='dot')

        colors = self.get_node_colors()

        nx.draw_networkx(self._G, pos, node_color=colors, **kwds)
        if edge_labels:
            # draw graph with weights on edges
            edge_labels = {(n1, n2): self._G[n1][n2]['path_length'] for (n1, n2) in self._G.edges()}
            nx.draw_networkx_edge_labels(self._G, pos, edge_labels=edge_labels, **kwds)

    def draw_2D(self, fig=None, projection='xz', axon_color='grey', dendrite_color='darkgrey',
                apical_dendrite_color='darkyellow', x_offset=0, y_offset=0, **kwargs):
        """
        Plots a 2D projection of the stick figure neuron.
        :param fig: Figure in which to draw. If None a new figure is created.
        :param projection: Default='xz'. Identifier of the ewo dimensional projection plane ['xz', 'xy', 'yz']
        :param axon_color: Color, default='grey'. Defines the color of the axon.
        :param dendrite_color: Color, default='darkgrey'. Defines the color of the dendrites.
        :param x_offset: double, default = 0, define offset for the first projection coordinate.
        :param y_offset: double, default = 0, define offset for the first projection coordinate.
        :param kwargs: further plotting parameters that can be passed to the plt.plot function.
        """
        if not fig:
            fig = plt.figure()
        if projection == 'xy':
            indices = [0,1]
        elif projection == 'xz':
            indices = [0,2]
        elif projection == 'yz':
            indices = [1,2]
        else:
            raise ValueError('projection %s is not defined.'% projection)

        G = self.get_graph()
        V = np.zeros((len(G.edges()), 2, 3))
        colors = ['k', axon_color, dendrite_color, apical_dendrite_color]

        cs = []
        for k, e in enumerate(G.edges()):
            n1 = G.node[e[0]]
            n2 = G.node[e[1]]
            V[k, 0] = n1['pos']
            V[k, 1] = n2['pos']
            cs.append(colors[int(n2['type']) - 1])

        cs = np.array(cs)
        x = indices[0]
        y = indices[1]

        plt_idx = np.array([cs_i == axon_color for cs_i in cs])
        if len(plt_idx.shape) > 1:
            plt_idx = plt_idx[:, 0]

        if plt_idx.any():
            _ = fig.gca().plot(V[plt_idx, :, x].T + x_offset, V[plt_idx, :, y].T + y_offset, c=axon_color, **kwargs)

        plt_idx = np.array([cs_i == dendrite_color for cs_i in cs])
        if len(plt_idx.shape) > 1:
            plt_idx = plt_idx[:, 0]

        if plt_idx.any():
            _ = fig.gca().plot(V[plt_idx, :, x].T + x_offset, V[plt_idx, :, y].T + y_offset, c=dendrite_color, **kwargs)

        plt_idx = np.array([cs_i == apical_dendrite_color for cs_i in cs])
        if len(plt_idx.shape) > 1:
            plt_idx = plt_idx[:, 0]

        if plt_idx.any():
            _ = fig.gca().plot(V[plt_idx, :, x].T + x_offset, V[plt_idx, :, y].T + y_offset, c=apical_dendrite_color,
                               **kwargs)


    ############# SAVING FUNCTIONS #####################


    def to_swc(self):
        """
        Write NeuronTree into swc file compatible pandas.DataFrame.
        :return: pandas.DataFrame with columns 'n', 'type', 'x', 'y', 'z', 'radius' and 'parent'
        """

        # create dataframe with graph data
        G = self._G
        ids = [int(k) for k in G.nodes()]
        pos = np.round(np.array(list(nx.get_node_attributes(G, 'pos').values())), 2)
        r = np.array(list(nx.get_node_attributes(G, 'radius').values()))
        t = np.array(list(nx.get_node_attributes(G, 'type').values())).astype(int)
        pids = [int(list(l.keys())[0]) for l in list(G.pred.values()) if list(l.keys())]
        pids.insert(0, -1)
        # write graph into swc file
        d = {'n': ids, 'type': t, 'x': pos[:, 0], 'y': pos[:, 1], 'z': pos[:, 2], 'radius': r, 'parent': pids}
        df = pd.DataFrame(data=d, columns=['n', 'type', 'x','y', 'z' ,'radius' , 'parent'])

        return df

    def write_to_swc(self, file_name,
                     path='/gpfs01/berens/data/data/anatomy/BC_morphologies/swc_tree/'):
        """
        Write NeuronTree to swc file.
        :param file_name: String. File name without '.swc' tag.
        :param path: String. Path to file

        """

        if not exists(path):
            makedirs(path)

        df = self.to_swc()
        df.to_csv(path + file_name + '.swc', sep=' ', encoding='utf-8', header=False, index=False)

    def write_to_mat(self, file_name,
                     path='/gpfs01/berens/data/data/anatomy/BC_morphologies/csv_luxburg/'):
        """
        Write NeuronTree to mat file.
        :param file_name: String. Filename without .mat tag.
        :param path:  String. Path to file
        """

        if not exists(path):
            makedirs(path)

        data = {}
        G = self._G
        P = list(nx.get_node_attributes(G, 'pos').values())
        T = list(nx.get_node_attributes(G, 'type').values())
        A = nx.adjacency_matrix(G, weight='path_length')

        data['pos'] = P
        data['type'] = T
        data['A'] = A
        savemat(path + file_name + '.mat', data)


