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
    swc.index = swc.n.as_matrix()

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

    """
    Turning graph into vectors for paths

    Parameters
    ----------
    G: nx.Graph()
        the graph representation of the cell morphology

    Returns
    -------
    path_all: dict
        a dict holding all paths, 
        which is the segment between two branchpoitns. 
    """

    edges_all = G.edge
    nodes_all = np.array(list(G.node.keys()))

    num_branch = np.array(list(G.out_degree().values()))
    branchpoints = nodes_all[np.logical_and(num_branch != 1, num_branch !=0)]
    if 1 not in branchpoints:
        branchpoints = np.hstack([1, branchpoints])

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
        type_dict[key] = np.vstack([G.node[key]['type'] for key in path_idx_dict[key]])[1][0]
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

# def connect_to_soma(current_path, soma):
#     """

#     :param current_path:
#     :param soma:
#     :return:
#     """

#     return (current_path == soma).all(1).any()

# def find_connection(all_paths, soma, path_id, paths_to_ignore=[]):
#     """

#     :param all_paths:
#     :param soma:
#     :param path_id:
#     :param paths_to_ignore:
#     :return:
#     """

#     current_path = all_paths[path_id]


#     if connect_to_soma(current_path, soma):

#         connect_to = -1
#         connect_to_at = soma

#         return connect_to, connect_to_at

#     sub_paths = deepcopy(all_paths)
#     sub_paths.pop(path_id)

#     for key in sub_paths.keys():


#         if key in paths_to_ignore: continue

#         target_path = sub_paths[key]
#         connect_to_at_loc = np.where((current_path[0] == target_path).all(1))[0]

#         if len(connect_to_at_loc) != 0:
#             connect_to = key
#             connect_to_at = target_path[connect_to_at_loc[0]]
#             return connect_to, connect_to_at

#     logging.info("Path {} connects to no other path. Try fix it.".format(path_id))
#     connect_to = -99
#     connect_to_at = np.nan

#     return connect_to, connect_to_at

# def back2soma(df_paths, path_id):
#     """
#     Given a path_id, find all the other paths which lead back to the soma.
#     :param df_paths:
#     :param path_id:
#     :return:
#     """

#     path_id_original = path_id

#     paths_to_soma = []

#     counter = 0

#     while df_paths.loc[path_id].connect_to != -1:
#         if path_id in paths_to_soma:
#             # logging.info("\tPath {} cannot trace back to soma: {}".format(path_id_original, paths_to_soma))
#             logging.info("\tPath {} cannot trace back to soma.".format(path_id_original))
#             break
#         else:
#             paths_to_soma.append(path_id)
#             path_id = df_paths.loc[path_id].connect_to

#         if path_id == -99:
#             break

#     paths_to_soma.append(path_id)

#     return paths_to_soma

# def check_path_connection(df_paths):
#     """
#     This function updates the df_paths object. It loops through all paths and adds the fours columns
#     ['connect_to', 'connect_to_at', 'connected_by', 'connected_by_at'] to df_paths.

#     :param df_paths: pandas.DataFrame containing paths.
#     :return:
#         df_paths: pandas.DataFrame
#         Updated df_paths with added columns ['connect_to', 'connect_to_at', 'connected_by', 'connected_by_at'] for
#         each path.
#     """

#     soma = df_paths[df_paths.type == 1]['path'][0]

#     logging.debug(soma)

#     all_paths = df_paths.path.to_dict()
#     all_keys = list(all_paths.keys())

#     # find which path the current path connects to.
#     connect_to_dict = {}
#     connect_to_at_dict = {}
#     for path_id in all_keys:
#         connect_to_dict[path_id], connect_to_at_dict[path_id] = find_connection(all_paths, soma, path_id, paths_to_ignore=[])
#     df_paths['connect_to'] = pd.Series(connect_to_dict)
#     df_paths['connect_to_at'] = pd.Series(connect_to_at_dict)

#     # find all paths connect to current path.
#     connected_by_dict = {}
#     connected_by_at_dict = {}
#     for path_id in all_keys:
#         connected_by_dict[path_id]    = df_paths[df_paths.connect_to == path_id].index.tolist()
#         connected_by_at_dict[path_id] = df_paths[df_paths.connect_to == path_id].connect_to_at.tolist()
#     df_paths['connected_by'] = pd.Series(connected_by_dict)
#     df_paths['connected_by_at'] = pd.Series(connected_by_at_dict)


#     # fix unexpected broken paths (e.g. branch points exist when there are no branching.)
#     df_paths[df_paths.connected_by.apply(len) == 1].index

#     # check if all paths can goes back to soma.
#     back_to_soma_dict = {}
#     for path_id in all_keys:
#         back_to_soma_dict[path_id] = back2soma(df_paths, path_id)

#         # logging.info('  All paths can be traced back to soma. It is a single tree.')

#     df_paths['back_to_soma'] = pd.Series(back_to_soma_dict)

#     return df_paths
def sort_path_direction(df_paths):
    
    df_paths = df_paths.copy()
    soma = df_paths.loc[0].path.flatten()    
    
        
    df_paths['connect_to'] = np.nan
    df_paths['connect_to_at'] = ''
    df_paths['connect_to_at'] = df_paths['connect_to_at'].apply(np.array)

    for row in df_paths.iterrows():

        path_id = row[0]
        path = row[1]['path']

        if (path[0] == soma).all():
            df_paths.set_value(path_id, 'connect_to', -1)
            df_paths.set_value(path_id, 'connect_to_at', soma)
            continue

        if (path[-1] == soma).all():
            df_paths.set_value(path_id, 'path', path[::-1])
            df_paths.set_value(path_id, 'connect_to', -1)
            df_paths.set_value(path_id, 'connect_to_at', soma)
    
    new_target_paths = list(df_paths[~np.isnan(df_paths.connect_to)].index) # seed the first round of paths to check
    
    logging.info('  num of paths connected to soma: {}'.format(len(new_target_paths)))
    
    while np.count_nonzero(~np.isnan(df_paths.connect_to)) != len(df_paths):
        
        all_checked_paths = list(df_paths[~np.isnan(df_paths.connect_to)].index)
        num_check_paths_before = len(all_checked_paths)
        logging.info("  num of paths checked: {}".format(num_check_paths_before))

        target_paths = new_target_paths
        new_target_paths = [] # empty the list to hold new target paths for next round
        
        for target_path_id in target_paths:

            target_path = df_paths.loc[target_path_id].path
            
            for row in df_paths.iterrows():

                path_id = row[0]
                path = row[1]['path']

                if path_id in all_checked_paths:
                    continue
                else:
    #                 print(path_id)
                    if (path[0] == target_path).all(1).any():

                        df_paths.set_value(path_id, 'connect_to', target_path_id)
                        df_paths.set_value(path_id, 'connect_to_at', target_path[np.where((path[0] == target_path).all(1))[0]])
                        new_target_paths.append(path_id)
                        continue
                    
                    if (path[-1] == target_path).all(1).any(): 
                        df_paths.set_value(path_id, 'path', path[::-1])
                        df_paths.set_value(path_id, 'connect_to', target_path_id)
                        df_paths.set_value(path_id, 'connect_to_at', target_path[np.where((path[-1] == target_path).all(1))[0]])
                        new_target_paths.append(path_id)
             
        num_check_paths_after = len(list(df_paths[~np.isnan(df_paths.connect_to)].index))
        
        if num_check_paths_before == num_check_paths_after:
            num_disconneted = len(df_paths) - num_check_paths_after
            logging.info('\tNumber of disconnected path(s): {}'.format(num_disconneted))
            break
    
    df_paths = df_paths.drop(df_paths[np.isnan(df_paths.connect_to)].index)
        
    # find all paths connect to current path.
    connected_by_dict = {}
    connected_by_at_dict = {}
    for path_id in df_paths.index:
        connected_by_dict[path_id]    = df_paths[df_paths.connect_to == path_id].index.tolist()
        connected_by_at_dict[path_id] = df_paths[df_paths.connect_to == path_id].connect_to_at.tolist()
    df_paths['connected_by'] = pd.Series(connected_by_dict)
    df_paths['connected_by_at'] = pd.Series(connected_by_at_dict)
    
    back_to_soma_dict = {}
    for path_id in df_paths.index:
        list_to_soma = [path_id]
        next_path_id = df_paths.loc[path_id].connect_to
        while next_path_id != -1:
            list_to_soma.append(next_path_id)
            next_path_id = df_paths.loc[next_path_id].connect_to
        back_to_soma_dict[path_id] = list_to_soma
    df_paths['back_to_soma'] = pd.Series(back_to_soma_dict)
    
    return df_paths

def get_path_statistics(df_paths):
    """

    Add path statistics (e.g. real/euclidean length of each path, ordering, index of paths back to soma...)

    Parameters
    ==========
    df_paths

    Returns
    =======
    a updated df_paths
    """

    logging.info('  Calculating path statistics (e.g. real length, branching order...)')

    df_paths = df_paths.copy()

    all_keys = df_paths.index

    real_length_dict = {}
    euclidean_length_dict = {}
    back_to_soma_dict = {}
    branch_order_dict = {}

    for path_id in all_keys:

        path = df_paths.loc[path_id].path

        real_length_dict[path_id] = get_path_real_length(path)
        euclidean_length_dict[path_id] = get_path_euclidean_length(path)
        branch_order_dict[path_id] = len(df_paths.loc[path_id].back_to_soma) - 1

    df_paths['real_length'] = pd.Series(real_length_dict)
    df_paths['euclidean_length'] = pd.Series(euclidean_length_dict)
    df_paths['branch_order'] = pd.Series(branch_order_dict)

    return df_paths

def get_density_data_of_type(neurites, soma):

    """
    A helper function to gether all summarized infomation

    Parameters
    ----------
    neurites: pandas.DataFrame
    soma: pandas.DataFrame

    Returns
    -------
    result tuple: 
        (type, asymmetry, neurites_radius, neurites_size)

    Z: numpy.array
        the density map of neurites, a (100, 100) matrix
        
    """
    
    if len(neurites) < 2:
        return None, None
    
    import cv2
    from scipy.stats import gaussian_kde
    from scipy.spatial import ConvexHull
    from scipy.ndimage.measurements import center_of_mass
    
    soma_coords = soma.path.as_matrix()[0].flatten()
    xy = (np.vstack(neurites.path)[:, :2] - soma_coords[:2]).T
    kernel = gaussian_kde(xy, bw_method='silverman')

    lim_max = int(np.ceil((xy.T).max() / 20) * 20)
    lim_min = int(np.floor((xy.T).min() / 20) * 20)
    lim = max(abs(lim_max), abs(lim_min))
    X, Y = np.mgrid[-lim:lim:100j, -lim:lim:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])

    Z = np.flipud(np.rot90(np.reshape(kernel(positions).T, X.shape)))

    density_center = np.array(center_of_mass(Z))
    density_center = density_center * (2 * lim) / Z.shape[0] - lim
    
    asymmetry = np.sqrt(np.sum(density_center ** 2))

    hull = ConvexHull(xy.T)
    outer_terminals = xy.T[hull.vertices]
    outer_terminals = np.vstack([outer_terminals, outer_terminals[0]])
    neurites_radius = np.mean(np.sqrt(np.sum((outer_terminals - density_center)**2, 1)))
    neurites_size = cv2.contourArea(outer_terminals.astype(np.float32))

    if neurites.iloc[0].type == 2:
        t = 'axon'
    elif neurites.iloc[0].type == 3:
        t = 'basal_dendrites'
    elif neurites.iloc[0].type == 4:
        t = 'apical_dendrites'
    else:
        t = 'undefined'

    return (t, asymmetry, neurites_radius, neurites_size), Z

def get_density_data(df_paths):

    """
    A helper function to gether all summarized infomation

    Parameters
    ----------
    df_paths: pandas.DataFrame

    Returns
    -------
    df_density: pandas.DataFrame

    density_maps: numpy.array
        a (3, 100, 100) matrix, each layer is a density map of one neurites type 
        (0: axon; 1: basal dendrites; 2: apical dendrites)

    """

    logging.info('  Calculating density data...')

    df_paths = df_paths.copy()

    soma = df_paths[df_paths.type == 1]
    axon = df_paths[df_paths.type == 2]
    dend_basal = df_paths[df_paths.type == 3]
    dend_apical = df_paths[df_paths.type == 4]

    axon_density_summary, axon_density_map = get_density_data_of_type(axon, soma)
    dend_basal_density_data, dend_basal_density_map = get_density_data_of_type(dend_basal, soma)
    dend_apical_density_data, dend_apical_density_map = get_density_data_of_type(dend_apical, soma)

    density_maps = np.zeros([3, 100, 100])
    density_maps[0] = axon_density_map
    density_maps[1] = dend_basal_density_map
    density_maps[2] = dend_apical_density_map

    labels = ['type', 'asymmetry', 'radius', 'size']
    neurites = [axon_density_summary,dend_basal_density_data,dend_apical_density_data]
    df_density = pd.DataFrame.from_records([n for n in neurites if n is not None], columns=labels)

    return df_density, density_maps