import numpy as np
import pandas as pd
import networkx as nx
from copy import deepcopy
import logging

from .summarize import *

def get_logger(loglevel):

    """
    Log out useful or debug infos.

    Parameters
    ---------
    loglevel: str
        'debug', 'info', 'warning', 'error', 'critical'. 

    Returns
    -------
    logger: logging.RootLogger
        a logging object for turning on and off the log.    
    
    """
    
    logger = logging.getLogger()
    
    LEVELS = {'debug': logging.DEBUG,
              'info': logging.INFO,
              'warning': logging.WARNING,
              'error': logging.ERROR,
              'critical': logging.CRITICAL}
    
    try:
        LEVEL = LEVELS[loglevel]
        logger.setLevel(LEVEL)
    except ValueError:
        logger.setLevel(logging.INFO)
        logging.info('  Please enter a valid logging mode (DEBUG, INFO, WARNING, ERROR, CRITICAL).')
        logger.setLevel(logging.ERROR)
        
    return logger

def read_swc(filepath):
    """
    Read swc into Pandas DataFrame

    Parameters
    ----------
    filepath: str
        path to swc file
    
    Returns
    -------
    G: networkx Graph object.
    """
    
    swc =  pd.read_csv(filepath, delim_whitespace=True, comment='#',
                          names=['n', 'type', 'x', 'y', 'z', 'radius', 'parent'], index_col=False)
    swc.index = df_swc.n.as_matrix()
    
    G = nx.DiGraph()
    
    # raw data
    n = swc['n'].tolist()
    pos = np.array([swc['x'].tolist(), swc['y'].tolist(), swc['z'].tolist()]).T
    radius = swc['radius'].tolist()
    t = swc['type'].tolist()
    pid = swc['parent'].tolist()
    t[pid == -1] = 1 # if soma is missing, the first point is soma
    
    # node
    node_keys = ['pos', 'type', 'radius']
    node_data = list(zip(n,
                        [dict(zip(node_keys, [pos[ix], t[ix], radius[ix]])) for ix in range(pos.shape[0])]))
    parent_idx = np.array([n.index(pid[ix]) for ix in range(1, len(pid))])

    # edge
    ec = np.sqrt(np.sum((pos[parent_idx] - pos[1:]) ** 2, axis=1)) 
    edge_keys =['euclidean_dist', 'path_length']
    edge_data = list(zip(pid[1:], n[1:],
                        [dict(zip(edge_keys, [ec[ix], ec[ix]])) for ix in range(ec.shape[0])]))

    G.add_nodes_from(node_data)
    G.add_edges_from(edge_data)
    
    return G, swc

def graph_to_path(G):
    
    edges_all = G.edge
    nodes_all = np.array(list(G.node.keys()))
    
    num_branch = np.array(list(G.out_degree().values()))
    branchpoints = nodes_all[np.logical_and(num_branch != 1, num_branch !=0)]
    
    path_all = {}
    i = 1
    for current_key in branchpoints: 
        next_keys = list(edges_all[current_key].keys())
        path = [current_key]
        for next_key in next_keys:
            while len(edges_all[next_key]) == 1:
                path.append(next_key)
                next_key += 1
            else:
                path.append(next_key)
                path_all[i] = np.array(path)
                path = [current_key]
            i+=1
            
    return path_all

def get_df_paths(G):
    
    """
    Split the original swc into paths (Soma, Dendrites, Axon, etc..). 
    
    Parameters
    ----------
    G : networkx Graph object.
    
    Returns
    -------
    df_paths: pandas.DataFrame
        A DataFrame with columns ['type', 'path', 'radius', 'n_index']
        * the first row (df.iloc[0]) is soma. 
        * the first point of each path should be the branch point.
    """ 

    path_idx_dict = graph_to_path(G)
    
    path_dict = {}
    type_dict = {}
    radius_dict = {}
    for key in path_idx_dict.keys():
        path_dict[key] = np.vstack([G.node[key]['pos'] for key in path_idx_dict[key]])
        type_dict[key] = np.vstack([G.node[key]['type'] for key in path_idx_dict[key]])[1]
        radius_dict[key] = np.vstack([G.node[key]['radius'] for key in path_idx_dict[key]])
    
    type_dict[0] = G.node[1]['type']
    path_dict[0] = G.node[1]['pos'].reshape(1,3)
    radius_dict[0] = [G.node[1]['radius']]
    path_idx_dict[0] = [1]
    
    
    df_paths = pd.DataFrame()
    df_paths['type'] = pd.Series(type_dict)
    df_paths['path'] = pd.Series(path_dict)
    df_paths['radius'] = pd.Series(radius_dict)
    df_paths['n_index'] = pd.Series(path_idx_dict)
    
    return df_paths

# def get_consecutive_pairs_of_elements_from_list(l, s=None, e=None):
    
#     """
#     Get pair of items in a list. 
    
#     Parameters
#     ----------
#     l : list or array-like 
#         e.g. [1,2,3,4]

#     s : int
#         An integer inserted to the front of the list or array.

#     e : int or None
#         An integer or None appended to the list or array.
#         `None` is always added unless other value is specified. 
#         This is for slicing the last part of the list. e.g. [4:None] is equivalent to [4:] 
    
#     Returns
#     -------
#     pair : list
#         a list of consecutive pairs of elements from the input list or array.
#         e.g. [(1,2), (2,3), (3,4), (4,None)]
#     """
#     if s is not None:
#         l = np.append(s, l)
#     l = np.append(l, e)
#     pair = list(zip(l[:-1], l[1:]))
    
#     return pair

# def preprocess_swc(df_swc):

#     """
    
#     Preprocess swc data into two parts: df_soma and df_neurites 
#     (although the soma point is also inclued in the neurites...)

#     The validity of soma points are also assessed. 

#     For 1-point soma, it will be retrieved without modification.

#     For 3 or multi-point soma, the AVERAGE of all soma points is used as the 1-point soma.
#     All points connected to soma points will be reconnected to point 1. 

#     """
    
#     df_swc = deepcopy(df_swc)
    
#     df_soma = df_swc[df_swc.type == 1]
#     if len(df_soma) > 1:
#         logging.info('  More then one points of soma are found. The average of all points is used as soma.\n')
        
#         # using the average of all soma points as THE soma point.
#         df_soma.set_value(1, 'x', df_soma['x'].mean())
#         df_soma.set_value(1, 'y', df_soma['y'].mean())
#         df_soma.set_value(1, 'z', df_soma['z'].mean())
#         df_soma.set_value(1, 'radius', df_soma['radius'].mean())
        
#         # reconnecting all other soma points to THE soma point
#         # also changing the parents of neurites which connected to soma to 1 
#         for i in range(2, len(df_soma)+1):
#             n_index = df_swc[df_swc.parent == i].index
#             for j in n_index:
#                 df_swc.set_value(j, 'parent', 1)
#         df_soma = df_swc.iloc[[0]]
        
#     elif len(df_soma) < 1:
#         logging.info('  No soma is founded. The first row is used as soma.\n')
#         df_swc.set_value(1, 'type', 1)
#         df_soma = df_swc[df_swc.type == 1]
# #     elif len(df_soma) == 1:
#     else:
#         logging.info('  One point soma is found.\n')
        
#     df_neurites = df_swc[df_swc.type != 1]
#     df_neurites = pd.concat([df_soma.iloc[[0]], df_neurites])
    
#     df_swc[['x', 'y', 'z']] = df_swc[['x', 'y', 'z']] - df_swc[['x', 'y', 'z']].iloc[0]

#     return df_swc, df_soma, df_neurites

# def get_df_paths(df_swc):
    
#     """
#     Split the original swc into paths (Soma, Dendrites, Axon, etc..). 
    
#     Parameters
#     ----------
#     df_swc : pandas.DataFrame
    
#     Returns
#     -------
#     df_paths: pandas.DataFrame
#         A DataFrame with columns ['type', 'path', 'radius', 'n_index']
#         * the first row (df.iloc[0]) is soma. 
#         * the first point of each path should be the branch point.
#     """ 

#     # df_swc, df_soma, df_neurites = preprocess_swc(df_swc)
#     df_soma = df_swc.iloc[[0]]
#     df_neurites = df_swc[df_swc.type != 1]
#     df_neurites = pd.concat([df_soma.iloc[[0]], df_neurites])
    
#     n = df_neurites.n.values
#     parent = df_neurites.parent.values
#     diff_n_parent = n - parent
#     diff_n_parent[1] = 1 # the first non-soma point should be included.
    
#     df_starting_points = df_neurites[diff_n_parent != 1] # starting point of each path, which is not the branch point.
#     branchpoint_index = np.unique(df_starting_points.parent.values[1:]) # branch point is the parent point of starting point.

#     path_dict = {}
#     type_dict = {}
#     radius_dict = {}
#     n_index_dict = {}
#     path_id = 1

#     starting_points_pairs = get_consecutive_pairs_of_elements_from_list(df_starting_points.n.values)

#     for s, e in starting_points_pairs: # s: start index, e: end index
#         logging.debug('{}, {}'.format(s, e))
#         if e is not None:
#             b = branchpoint_index[np.logical_and(branchpoint_index > s, branchpoint_index < e)] # b: list, branch index
#         else:
#             b = branchpoint_index[branchpoint_index > s]

#         logging.debug("  \t{}".format(b))
#         branchpoint_index_pairs = get_consecutive_pairs_of_elements_from_list(b, s, e)

#         for bs, be in branchpoint_index_pairs:
            
#             logging.debug('\t {}: {}, {}'.format(path_id, bs, be))
            
#             path = df_neurites.loc[bs:be][['x', 'y', 'z']].values

#             logging.debug('  Lenght of path is {}'.format(len(path)))

#             path_type = df_neurites.loc[bs:be][['type']].values
#             path_radius = df_neurites.loc[bs:be][['radius']].values
#             path_n_index = df_neurites.loc[bs:be][['n']].values
            
#             if bs == s:
#                 parent_idx = df_neurites.loc[[bs]].parent.values[0]    
#                 if parent_idx != -1:
#                     parent_coord = df_neurites.loc[[parent_idx]][['x', 'y', 'z']].values
#                     parent_radius = df_neurites.loc[[parent_idx]][['radius']].values
#                     if not (parent_coord == path).all(1).any():
#                         path = np.vstack([parent_coord, path])
#                         path_radius = np.vstack([parent_radius, path_radius])

#             if be == e:

#                 if len(path[:-1]) <= 2: continue # remove a path contains just one point 

#                 path_dict[path_id] = path[:-1]
#                 radius_dict[path_id] = path_radius[:-1]
#                 n_index_dict[path_id] = path_n_index[:-1]
#             else:
#                 path_dict[path_id] = path
#                 radius_dict[path_id] = path_radius
#                 n_index_dict[path_id] = path_n_index

#             type_dict[path_id] = max(path_type)[0]
#             path_id += 1
            
#     type_dict[0] = df_soma['type'].as_matrix()[0]
#     path_dict[0] = df_soma[['x', 'y', 'z']].as_matrix()
#     radius_dict[0] = df_soma['radius'].as_matrix()
#     n_index_dict[0] = df_soma['n'].as_matrix()
    
#     df_paths = pd.DataFrame()
#     df_paths['type'] = pd.Series(type_dict)
#     df_paths['path'] = pd.Series(path_dict)
#     df_paths['radius'] = pd.Series(radius_dict)
#     df_paths['n_index'] = pd.Series(n_index_dict)
    
#     return df_paths

def swc_to_linestack(df_swc, voxelsize=None):
    
    """
    Convert SWC to Line Stack (from real length to voxel coordinates).
    :param df_swc:
    :param unit:
    :param voxelsize:
    :return:
    """

    coords = df_swc[['x', 'y', 'z']].as_matrix()

    if voxelsize is None:
        logging.debug('  Not able to build linestack from real length coordinates.')
        return None
    else:
        # coords = np.round(coords / voxelsize).astype(int) # not able to handle the soma-centered swc.
        logging.debug('  Real length coordindates are converted back to pixel.')
        coords = coords - coords.min(0)
        coords = np.round(coords / voxelsize).astype(int)   

    imagesize = coords.max(0) + 1
    logging.debug('  Start: Creating linestack...')
    linestack = np.zeros(imagesize)
    for c in coords:
        linestack[tuple(c)] = 1
    
    reset_neg = coords.min(0)
    reset_neg[reset_neg > 0] = 0
    
    xyz = (coords - reset_neg).max(0) + 1
    xy = max(xyz[:2])
    xy = np.ceil(xy/100) * 100
    z = xyz[2]
    z = np.ceil(z / 10) * 10
    logging.debug("{}, {}".format([xy, xy, z], linestack.shape))
    padding_cp = np.ceil((np.array([xy, xy, z]) - linestack.shape) / 2).astype(int)
    padding_x = padding_cp.copy()
    padding_y = padding_cp.copy()
    
    odd = np.array(linestack.shape) % 2 == 1
    padding_y[odd] = padding_y[odd] - 1
    
    padding = np.vstack([padding_x, padding_y]).T

    npad = ((padding[0]), (padding[1]), (padding[2]))
    linestack = np.pad(linestack, pad_width=npad, mode='constant')
    soma_on_stack = coords[0] + padding_x

    logging.debug('  Finished.\n')
    
    return linestack, soma_on_stack, padding_x


def connect_to_soma(current_path, soma):
    """

    :param current_path:
    :param soma:
    :return:
    """
    # return (current_path[0] == soma).all(1).any()
    # logging.debug(current_path == soma)
    return (current_path == soma).all(1).any()

def find_connection(all_paths, soma, path_id, paths_to_ignore=[]):
    """

    :param all_paths:
    :param soma:
    :param path_id:
    :param paths_to_ignore:
    :return:
    """
    
    current_path = all_paths[path_id]
    
    
    if connect_to_soma(current_path, soma):

        connect_to = -1
        connect_to_at = soma
        
        return connect_to, connect_to_at
        
    sub_paths = deepcopy(all_paths)
    sub_paths.pop(path_id)
    
    for key in sub_paths.keys():
        
        
        if key in paths_to_ignore: continue
        
        target_path = sub_paths[key]
        connect_to_at_loc = np.where((current_path[0] == target_path).all(1))[0] 
        
        if len(connect_to_at_loc) != 0:
            connect_to = key
            connect_to_at = target_path[connect_to_at_loc[0]]
            return connect_to, connect_to_at
    
    logging.info("Path {} connects to no other path. Try fix it.".format(path_id))
    connect_to = -99 
    connect_to_at = np.nan
    
    return connect_to, connect_to_at

def back2soma(df_paths, path_id):
    """
    Given a path_id, find all the other paths which lead back to the soma.
    :param df_paths:
    :param path_id:
    :return:
    """

    path_id_original = path_id
    
    paths_to_soma = []
    
    counter = 0
        
    while df_paths.loc[path_id].connect_to != -1:
        if path_id in paths_to_soma:
            # logging.info("\tPath {} cannot trace back to soma: {}".format(path_id_original, paths_to_soma))
            logging.info("\tPath {} cannot trace back to soma.".format(path_id_original))
            break
        else:
            paths_to_soma.append(path_id)
            path_id = df_paths.loc[path_id].connect_to   

        if path_id == -99:
            break

    paths_to_soma.append(path_id)  

    return paths_to_soma

def check_path_connection(df_paths):
    """
    This function updates the df_paths object. It loops through all paths and adds the fours columns
    ['connect_to', 'connect_to_at', 'connected_by', 'connected_by_at'] to df_paths.

    :param df_paths: pandas.DataFrame containing paths.
    :return:
        df_paths: pandas.DataFrame
        Updated df_paths with added columns ['connect_to', 'connect_to_at', 'connected_by', 'connected_by_at'] for
        each path.
    """

    soma = df_paths[df_paths.type == 1]['path'][0]
    
    logging.debug(soma)

    all_paths = df_paths.path.to_dict()
    all_keys = list(all_paths.keys())

    # find which path the current path connects to.
    connect_to_dict = {}
    connect_to_at_dict = {}
    for path_id in all_keys:
        connect_to_dict[path_id], connect_to_at_dict[path_id] = find_connection(all_paths, soma, path_id, paths_to_ignore=[])
    df_paths['connect_to'] = pd.Series(connect_to_dict)
    df_paths['connect_to_at'] = pd.Series(connect_to_at_dict)

    # find all paths connect to current path.  
    connected_by_dict = {}
    connected_by_at_dict = {}
    for path_id in all_keys:
        connected_by_dict[path_id]    = df_paths[df_paths.connect_to == path_id].index.tolist()
        connected_by_at_dict[path_id] = df_paths[df_paths.connect_to == path_id].connect_to_at.tolist()
    df_paths['connected_by'] = pd.Series(connected_by_dict)
    df_paths['connected_by_at'] = pd.Series(connected_by_at_dict)


    # fix unexpected broken paths (e.g. branch points exist when there are no branching.)
    df_paths[df_paths.connected_by.apply(len) == 1].index

    # check if all paths can goes back to soma.
    back_to_soma_dict = {}
    for path_id in all_keys:
        back_to_soma_dict[path_id] = back2soma(df_paths, path_id)
    
        # logging.info('  All paths can be traced back to soma. It is a single tree.')

    df_paths['back_to_soma'] = pd.Series(back_to_soma_dict)

    return df_paths

def unique_row(a):

    """
    Returns an array of the ordered, unique set of rows for input array a.

    Parameters
    ----------
    a: array
        an array with replicated rows.

    returns
    -------
    unique_a: array
        an ordered array without replicated rows.

    example
    -------
    >>> a = np.array([[9,9],
                      [8,8],
                      [1,1],
                      [9,9]])
    >>> unique_row(a)
    >>> array([[1, 1],
               [8, 8],
               [9, 9]])
    """

    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    _, idx = np.unique(b, return_index=True)
    
    unique_a = a[idx]
    
    return unique_a

def get_path_statistics(df_paths):
    """

    Add path statistics (e.g. dendritic/euclidean lenght of each path, ordering, index of paths back to soma...)

    Parameters
    ==========
    df_paths

    Returns
    =======
    a updated df_paths
    """
    
    logging.info('  Start: Calculating path statistics (e.g. dendritic length, branching order...)')

    all_keys = df_paths.index
    
    dendritic_length_dict = {}
    euclidean_length_dict = {}
    back_to_soma_dict = {}
    corder_dict = {}
    
    for path_id in all_keys:
        
        path = df_paths.loc[path_id].path
        
        dendritic_length_dict[path_id] = get_path_dendritic_length(path)
        euclidean_length_dict[path_id] = get_path_euclidean_length(path)
        corder_dict[path_id] = len(df_paths.loc[path_id].back_to_soma)

    df_paths['dendritic_length'] = pd.Series(dendritic_length_dict)
    df_paths['euclidean_length'] = pd.Series(euclidean_length_dict)
    df_paths['corder'] = pd.Series(corder_dict)

    # # Strahler order
    # df_paths['sorder'] = np.ones(len(df_paths)) * np.nan
    # df_paths = df_paths.set_value(df_paths.connected_by.apply(len) == 0, 'sorder', 1)
    
    # while np.isnan(df_paths.sorder).any():
    
    #     df_sub = df_paths[np.isnan(df_paths.sorder)]

    #     for row in df_sub.iterrows():

    #         path_id = row[0]
    #         connected_by = row[1]['connected_by']

    #         sorder0 = df_paths.loc[connected_by[0]].sorder
    #         sorder1 = df_paths.loc[connected_by[1]].sorder

    #         if np.isnan(sorder0) or np.isnan(sorder1):
    #             continue
    #         else:
                
    #             if sorder0 == sorder1:
    #                 df_paths.set_value(path_id, 'sorder', sorder0+1)
    #             else:
    #                 df_paths.set_value(path_id, 'sorder', np.max([sorder0, sorder1]))

    # df_paths.sorder = df_paths['sorder'].astype(int)

    logging.info('  Done. \n')
    
    return df_paths

def calculate_density(linestack, voxelsize):
    """
    reference: A. Stepanyantsa & D.B. Chklovskiib (2005). Neurogeometry and potential synaptic connectivity. Trends in Neuroscience.
    :param linestack:
    :param voxelsize:
    :return:
    """

    import scipy.ndimage as ndimage

    logging.debug('  Start: Calculating dendritic density...')

    smoothed_layer_stack = []
    for i in range(linestack.shape[2]):

        layer = linestack[:,:,i]
        smoothed_layer_stack.append(ndimage.gaussian_filter(layer, sigma=25/voxelsize[0], order=0))

    density_stack = np.dstack(smoothed_layer_stack)
    center_of_mass = np.array(ndimage.measurements.center_of_mass(density_stack.sum(2)))

    logging.debug('  Finished. \n')

    return density_stack, center_of_mass

def get_path_on_stack(df_paths, voxelsize, coordinate_padding):
    """

    :param df_paths:
    :param voxelsize:
    :param coordinate_padding:
    :return:
    """

    if voxelsize is None:
        voxelsize = np.array([1,1,1])

    path_dict = df_paths.path.to_dict()
    path_stack_dict = {}

    all_keys = path_dict.keys()

    coords = np.vstack(df_paths.path)
    reset_neg = coords.min(0)
    reset_neg[reset_neg > 0] = 3 
    
    for path_id in all_keys:

        path = path_dict[path_id]
        path = path - reset_neg
        path = np.round(path / voxelsize).astype(int)    

        path_stack_dict[path_id] = path + coordinate_padding

    df_paths['path_stack'] = pd.Series(path_stack_dict)

    cols = list(df_paths.columns)
    cols.remove('path_stack')
    cols.insert(1, 'path_stack')
    
    return df_paths[cols]


def print_summary(summary):
    
    """
    Print out summary statistics of the cell.

    Parameters
    ----------
    summary: dict
        a nested dict that contains summary of the cell.
    """

    import logging
    
    num_dendritic_segments = summary['general']['number_of_dendritic_segments']
    num_branchpoints = summary['general']['number_of_branch_points']
    num_irreducible_nodes = summary['general']['number_of_irreducible_nodes']
    max_branch_order = summary['general']['max_branch_order']
    # max_strahler_order = summary['general'] ['max_strahler_order']
    average_nodal_angle_deg = summary['angle']['average_nodal_angle_in_degree']
    average_nodal_angle_rad = summary['angle']['average_nodal_angle_in_radian']
    average_local_angle_deg = summary['angle']['average_local_angle_in_degree']
    average_local_angle_rad = summary['angle']['average_local_angle_in_radian']
    average_tortuosity = summary['length']['tortuosity']
    dendritic_sum = summary['length']['dendritic']['sum']
    dendritic_mean = summary['length']['dendritic']['mean']
    dendritic_median = summary['length']['dendritic']['median']
    dendritic_min = summary['length']['dendritic']['min']
    dendritic_max = summary['length']['dendritic']['max']
    euclidean_sum = summary['length']['euclidean']['sum']
    euclidean_mean = summary['length']['euclidean']['mean']
    euclidean_median = summary['length']['euclidean']['median']
    euclidean_min = summary['length']['euclidean']['min']
    euclidean_max = summary['length']['euclidean']['max']
    
    logging.info('  Summary of the cell')
    logging.info('  ======================\n')

    logging.info('  # General Infomation\n')
    logging.info('    Number of dendritic arbor segment: {}'.format(num_dendritic_segments))
    logging.info('    Number of branch points: {}'.format(num_branchpoints))
    logging.info('    Number of irreducible nodes: {}\n'.format(num_irreducible_nodes))

    logging.info('    Max branching order: {}'.format(max_branch_order))
    # logging.info('    Max Strahler order: {}\n\n'.format(max_strahler_order))

    logging.info('  # Angle \n')
    logging.info('    Average nodal angle in degree: {:.3f}'.format(average_nodal_angle_deg))
    logging.info('    Average nodal angle in radian: {:.3f} \n'.format(average_nodal_angle_rad))
    logging.info('    Average local angle in degree: {:.3f}'.format(average_local_angle_deg))
    logging.info('    Average local angle in radian: {:.3f} \n'.format(average_local_angle_rad))

    logging.info('  # Average tortuosity: {:.3f}\n'.format(average_tortuosity))

    logging.info('  # Dendritic length (μm)\n')
    # logging.info('  ## Dendritic length\n')
    logging.info('    Sum: {:.3f}'.format(dendritic_sum))
    logging.info('    Mean: {:.3f}'.format(dendritic_mean))
    logging.info('    Median: {:.3f}'.format(dendritic_median))
    logging.info('    Min: {:.3f}'.format(dendritic_min))
    logging.info('    Max: {:.3f}\n'.format(dendritic_max))

    logging.info('  # Euclidean length (μm)\n')
    # logging.info('  ## Euclidean length\n')

    logging.info('    Sum: {:.3f}'.format(euclidean_sum))
    logging.info('    Mean: {:.3f}'.format(euclidean_mean))
    logging.info('    Median: {:.3f}'.format(euclidean_median))
    logging.info('    Min: {:.3f}'.format(euclidean_min))
    logging.info('    Max: {:.3f}\n'.format(euclidean_max))    
    
    if 'density' in summary.keys():
        
        asymmetry = summary['density']['asymmetry']
        outer_radius = summary['density']['outer_radius']
        typical_radius = summary['density']['typical_radius']
        dendritic_area = summary['density']['dendritic_area']
        
        logging.info('  # Density related (μm)\n')
        logging.info('    Asymmetry: {:.3f}'.format(asymmetry))
        logging.info('    Outer Radius: {:.3f}'.format(outer_radius))
        logging.info('    Typical Radius : {:.3f}\n'.format(typical_radius))
        logging.info('    Dendritic Area: {:.3f} ×10\u00b3 um\u00b2\n\n'.format(dendritic_area))