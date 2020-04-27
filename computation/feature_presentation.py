import copy
import random
import collections
import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from neurontree.utils import smooth_gaussian
from computation.persistence_functions import radial_distance


def get_persistence(neurontree=None, f=None):
    """
    Creates the persistence barcode for the graph G. The algorithm is taken from
    _Quantifying topological invariants of neuronal morphologies_ from Lida Kanari et al
    (https://arxiv.org/abs/1603.08432).
    changed for use with networkx v2 (works also in old version: list(G.neighbors()))
    :param neurontree: instance of a NeuronTree class which holds the data of the swc file
    :param f: user defined function for computing persitence (see persistence_functions.py)
    :return: pandas.DataFrame with entries node_id | birth | death . Where birth and death are defined in radial
        distance from soma.
    """
    # Initialization
    L = neurontree.get_tips()
    R = neurontree.get_root()
    G = neurontree.get_graph()
    D = dict(node_id=[], node_type=[], birth=[], death=[])  # holds persistence barcode
    v = dict()  # holds 'aging' function of visited nodes defined by f

    # active nodes
    A = list(copy.copy(L))

    if f is None:
        # radial distance function
        f = radial_distance

    # set the initial value for leaf nodes
    for l in L:
        v[l] = f(G, l, R)

    while R not in A:
        for l in A:
            p = list(G.predecessors(l))[0]
            C = list(G.successors(p))

            # if all children are active
            if all(c in A for c in C):
                # choose randomly from the oldest children
                age = np.array([np.abs(v[c]) for c in C])
                indices = np.where(age == age[np.argmax(age)])[0]
                c_m = C[random.choice(indices)]

                A.append(p)

                for c_i in C:
                    A.remove(c_i)
                    if c_i != c_m:
                        D['node_id'].append(c_i)
                        if neurontree._nxversion == 2:
                            D['node_type'].append(G.nodes[c_i]['type'])
                        else:
                            D['node_type'].append(G.node[c_i]['type'])
                        D['birth'].append(v[c_i])
                        D['death'].append(f(G, p, R))
                v[p] = v[c_m]
    D['node_id'].append(R)
    D['node_type'].append(1)
    D['birth'].append(v[R])
    D['death'].append(f(G, R, R))
    return pd.DataFrame(D)


def compute_morphometric_statistics(neurontree=None, format='wide'):
    """
    Compute various morphometric statistics of a NeuronTree which is passed as an object
    :param neurontree: NeuronTree instance, holds complete data of an swc file
    :param format: String (default='wide'), determines the data format of the returned statistics. Options are 'wide' and
    'long'. For more information read http://www.cookbook-r.com/Manipulating_data/Converting_data_between_wide_and_long_format/.
    :return: pandas dataframe object with dictionary of all statistics
    """
    if neurontree is None:
        return None

    z = dict()
    z['branch_points'] = neurontree.get_branchpoints().size

    extend = neurontree.get_extend()
    z['width'] = extend[0]
    z['depth'] = extend[1]
    z['height'] = extend[2]

    tips = neurontree.get_tips()
    z['tips'] = tips.size

    root = neurontree.get_root()
    z['stems'] = len(neurontree.edges(root))

    z['total_length'] = np.sum(list(neurontree.get_edge_attributes('path_length').values()))
    # get all radii
    radii = neurontree.get_node_attributes('radius')
    # delete the soma
    radii.pop(root)
    z['avg_thickness'] = np.mean(list(radii.values()))
    z['max_thickness'] = np.max(list(radii.values()))

    z['total_surface'] = np.sum(list(neurontree.get_surface().values()))
    z['total_volume'] = np.sum(list(neurontree.get_volume().values()))

    z['max_path_dist_to_soma'] = np.max(list(neurontree.get_path_length().values()))
    z['max_branch_order'] = np.max(list(neurontree.get_branch_order().values()))

    path_angles = list(neurontree.get_path_angles().values())

    z['max_path_angle'] = np.percentile(path_angles, 99.5)
    z['min_path_angle'] = np.min(path_angles)
    z['median_path_angle'] = np.median(path_angles)

    z['mean_soma_exit_angle'] = np.mean(neurontree.get_soma_angles())

    R = neurontree.get_topological_minor()
    segment_length = R.get_segment_length()
    terminal_segment_pl = [item[1] for item in segment_length.items() if item[0][1] in tips]
    intermediate_segment_pl = [item[1] for item in segment_length.items() if item[0][1] not in tips]

    z['max_segment_path_length'] = np.max(list(segment_length.values()))
    z['median_intermediate_segment_pl'] = np.median([0] + intermediate_segment_pl)
    z['median_terminal_segment_pl'] = np.median(terminal_segment_pl)

    tortuosity = [e[2]['path_length'] / e[2]['euclidean_dist'] for e in R.edges(data=True)]

    z['log_max_tortuosity'] = np.log(np.percentile(tortuosity, 99.5))
    z['log_min_tortuosity'] = np.log(np.min(tortuosity))
    z['log_median_tortuosity'] = np.log(np.median(tortuosity))

    branch_angles = list(R.get_branch_angles().values())
    z['max_branch_angle'] = np.max(branch_angles)
    z['min_branch_angle'] = np.min(branch_angles)
    z['mean_branch_angle'] = np.mean(branch_angles)

    # get maximal degree within data
    if neurontree._nxversion == 2:
        # changed for version 2.x of networkX
        z['max_degree'] = np.max([item[1] for item in R.get_graph().out_degree() if item[0] != R.get_root()])
    else:
        z['max_degree'] = np.max([item[1] for item in R.get_graph().out_degree().items() if item[0] != R.get_root()])
    
    # get tree asymmetry
    weights, psad = R.get_psad()
    if np.sum(list(weights.values())) != 0:
        z['tree_asymmetry'] = np.sum([weights[k] * psad[k] for k in psad.keys()]) / np.sum(list(weights.values()))
    else:
        z['tree_asymmetry'] = 0


    morph = pd.DataFrame.from_dict(z, orient='index').T
    if format == 'long':
        morph = morph.T.reset_index().rename(columns={0: 'value', 'index': 'statistic'})
    return morph


def compute_density_maps(neurontree=None, config_params=None):
    """
    function for computing density maps which can be specified by a config and passed with a neurontree
    several projections are computed: x,y,z,xy,xz,yz
    :param neurontree:      NeuronTree object wich holds a swc file data
    :param config_params:   configuration params passed as dictionary which was load from file
                            containing all customizable params for density maps
    :return:                returns a dictionary with all density maps computed with all projections (x,y,z,xy,xz,yz)
    """
    # read distance from config and set default if no config available:
    if config_params is None:
        distance = 1
    else:
        distance = config_params.get('distance', 1)

    # get the resampled point could along each neurite
    # at specific distance default: 1 micron.
    # pc is an array of 3D coordinates for each resampled node
    pc = neurontree.resample_nodes(d=distance)

    ###### PARAMETER ################

    # read all missing params from config and set default values if no config available:
    if config_params is None:
        density = True
        smooth = True
        sigma = 1
        min = np.min(pc, axis=0)
        max = np.max(pc, axis=0)

        bin_size = 20
        n_bins_x, n_bins_y, n_bins_z = np.ceil((max-min)/bin_size).astype(int)

    else:

        min_ = np.min(pc, axis=0)
        max_ = np.max(pc, axis=0)
        # if config available use params else default values
        density = config_params.get('density', True)
        smooth = config_params.get('smooth', True)
        sigma = config_params.get('sigma', 1)

        # normalization ranges
        r_min_x = config_params.get('r_min_x', min_[0])
        r_min_y = config_params.get('r_min_y', min_[1])
        r_min_z = config_params.get('r_min_z', min_[2])

        r_max_x = config_params.get('r_max_x', max_[0])
        r_max_y = config_params.get('r_max_y', max_[1])
        r_max_z = config_params.get('r_max_z', max_[2])

        min = np.array([r_min_x, r_min_y, r_min_z])
        max = np.array([r_max_x, r_max_y, r_max_z])


        bin_size = config_params.get('bin_size', 20)
        n_bins_x, n_bins_y, n_bins_z = np.ceil((max - min) / bin_size).astype(int)
        if 'n_bins_x' in config_params.keys():
            n_bins_x = config_params.get('n_bins_x')
        if 'n_bins_y' in config_params.keys():
            n_bins_y = config_params.get('n_bins_y')
        if 'n_bins_z' in config_params.keys():
            n_bins_z = config_params.get('n_bins_z')

    # dictionary for axes and all labels of projection in right order
    axes = collections.OrderedDict([('0', 'x'), ('1', 'y'), ('2', 'z'), ('01', 'xy'), ('02', 'xz'), ('12', 'yz')])

    # holds binning per projection
    bins = {'x': (n_bins_x,), 'y': (n_bins_y,), 'z': (n_bins_z,),
            'xy': (n_bins_x, n_bins_y), 'xz': (n_bins_x, n_bins_z), 'yz': (n_bins_y, n_bins_z)}

    ######## COMPUTATION ############
    # create ranges
    ranges = {'x': [[min[0], max[0]]], 'y': [[min[1], max[1]]], 'z': [[min[2], max[2]]],
              'xy': [[min[0], max[0]], [min[1], max[1]]],
              'xz': [[min[0], max[0]], [min[2], max[2]]],
              'yz': [[min[1], max[1]], [min[2], max[2]]]}

    # all computed density maps will be stored in a dictionary
    densities = collections.OrderedDict()
    # loop over all axes
    for p_ax, ax in axes.items():

        dim = len(p_ax)
        # holds the range for binning of the histogram. So far the cells are noramlized to be between max --> 1 and min --> 0
        # I can therefore know, that the point will lie between 0 and 1. However, the range could also be a parameter set
        # in the config file.
        data = _project_data(p_ax, pc)
        n_bins = bins[ax]
        ranges_ = ranges[ax]

        # compute histogram hence density map
        h, edges = np.histogramdd(data, bins=n_bins, range=ranges_, density=density)
        # perform smoothing
        if smooth:
            h = smooth_gaussian(h, dim=dim, sigma=sigma)
        densities['%s_proj' %ax] = {'data': h, 'edges': edges, 'bins':n_bins}

    return densities


def plot_density_maps(densities=None, figure=None):
    """
        functions to plot density maps from densities dictionary with data from x,y,z,xy,xz,yz projections

        :return:            figure will be returned with all plotted maps from densities
        :param densities:   dictionary which holds all projections for the plots
        :param figure:      you can pass a figure if you want to use a custom plot format
    """
    # if figure is not passed create a new one
    if figure is None:
        figure = plt.figure(figsize=(8,8))

    if densities is None:
        return figure

    # set the axes layout right
    ax_x = plt.subplot2grid((4, 4), (0, 1), rowspan=1, colspan=2)
    ax_y = plt.subplot2grid((4, 4), (1, 3), rowspan=2, colspan=1)
    ax_z = plt.subplot2grid((4, 4), (3, 3), rowspan=1, colspan=1)

    ax_xy = plt.subplot2grid((4, 4), (1, 1), rowspan=2, colspan=2, sharex=ax_x, sharey=ax_y)
    ax_xz = plt.subplot2grid((4, 4), (3, 1), rowspan=1, colspan=2, sharex=ax_xy, sharey=ax_z)
    ax_yz = plt.subplot2grid((4, 4), (1, 0), rowspan=2, colspan=1, sharey=ax_xy)
    axes = {'x': ax_x, 'y': ax_y, 'z': ax_z, 'xy': ax_xy, 'xz': ax_xz, 'yz': ax_yz}


    # loop over all densities, keys contain type of data
    for name, density in densities.items():
        # get name and split projection axes from it
        p_axes = name.split('_')[0]
        dim = len(p_axes)

        if dim > 1:
            ax = axes[p_axes]
            x_edges = density['edges'][0]
            y_edges = density['edges'][1]

            x_min = np.floor(np.min(x_edges))
            x_max = np.ceil(np.max(x_edges))
            y_min = np.floor(np.min(y_edges))
            y_max = np.ceil(np.max(y_edges))

            if p_axes == 'yz':
                ax.imshow(density['data'], extent=(y_min, y_max, x_max, x_min))
                ax.set_xlabel(p_axes[1].capitalize() + r' ($\mu$m)')
                ax.set_ylabel(p_axes[0].capitalize() + r' ($\mu$m)')
            else:
                ax.imshow(density['data'].T, extent=(x_min, x_max, y_max, y_min))
                ax.set_xlabel(p_axes[0].capitalize() + r' ($\mu$m)')
                ax.set_ylabel(p_axes[1].capitalize() + r' ($\mu$m)')
            ax.invert_yaxis()

        else:
            ax = axes[p_axes]
            if p_axes == 'x':
                ax.plot(density['edges'][0][:-1], density['data'])
                ax.set_xlabel(p_axes.capitalize() + r' ($\mu$m)')
                ax.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0))
            else:
                ax.plot(density['data'], density['edges'][0][:-1])
                ax.set_ylabel(p_axes.capitalize() + r' ($\mu$m)')
                ax.ticklabel_format(axis='x', style='scientific', scilimits=(0, 0))

            sns.despine()
        ax.set_aspect('auto')
    figure.tight_layout(rect=[0, 0.03, 1, 0.9])
    return figure


def _project_data(proj_axes, data):
    """
        Helper function to project data onto the the axes defined in proj_axes.

        :param proj_axes:   str that holds the axes that are projected to as number, e.g. '01' for projection onto xy
                            or '02' for projection onto xz.
        :param data:
    """
    p_a = proj_axes
    dim = len(proj_axes)

    if dim == 2:
        indices = '012'
        for ix in range(len(p_a)):
            indices = indices.replace(p_a[ix], '')
        deleted_axis = int(indices)
        ax = [0, 1, 2]
        ax.remove(deleted_axis)
        result = data[:, ax]

    elif dim == 1:

        ax = int(p_a)
        result = data[:, ax]
    else:
        result = data

    return result
