import numpy as np
import pandas as pd
from copy import deepcopy
import logging
import matplotlib.pyplot as plt


def get_logger(loglevel):

    """
    Log out useful or debug infos.

    Paramters
    ---------
    loglevel: str
        'debug', 'info', 'warning', 'error', 'critical'. 
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
    except:
        logger.setLevel(logging.INFO)
        logging.info('  Please enter a valid logging mode (DEBUG, INFO, WARNING, ERROR, CRITICAL).')
        logger.setLevel(logging.ERROR)
        
    return logger


def read_swc(filepath, unit, voxelsize):
    
    """
    Read swc into Pandas DataFrame

    Parameters
    ----------
    filepath: str
        path to swc file
    unit: str
        unit of the swc
    voxelsize: numpy.array, shape = (3,)
        the voxel separation (x,y,z)
    
    Returns
    -------
    df_swc: pandas.DataFrame
        A Pandas DataFrame, with columns ['n', 'type', 'x', 'y', 'z', 'radius', 'parent']
        * 'n' is the index of current row, 'parent' is the row current row connects to
        * 'x', 'y', 'z' are spatial the coordinates
        * 'type' is morphological structure identifier, 0-undefined, 1-soma, 2-axon, 3-dendrite, 4-apical dendrite, 5-custom.
    """
    
    df_swc =  pd.read_csv(filepath, delim_whitespace=True, comment='#',
                          names=['n', 'type', 'x', 'y', 'z', 'radius', 'parent'], index_col=False)
    df_swc.index = df_swc.n.as_matrix()
    if unit == 'pixel' and voxelsize is not None:
        df_swc[['x', 'y', 'z']] = df_swc[['x', 'y', 'z']] * voxelsize
        unit = 'um'
    
    return df_swc


def get_consecutive_pairs_of_elements_from_list(l, s=None, e=None):
    
    """
    Get pair of items in a list. 
    
    Argument `e=
    e.g. [4:None] is equivalent to [4:].    

    Argument
    
    Parameters
    ----------
    l : list
        e.g. [1,2,3,4]

    s : int
        An integer inserted to the front of a list.


    e : int or None
        An integer or None appended to the list.
        `None` is always added unless other value is specified. 
        This is for slicing the last part of the list. e.g. [4:None] is equivalent to [4:] 
    
    Returns
    -------
    pair : list
        e.g. [(1,2), (2,3), (3,4), (4,None)]
    """
    if s is not None:
        l = np.append(s, l)
    l = np.append(l, e)
    pair = list(zip(l[:-1], l[1:]))
    
    return pair


def get_df_paths(df_swc):
    
    """
    Split swc into paths. 
    
    Parameters
    ----------
    df_swc : pandas.DataFrame
    
    Returns
    -------
    df_paths: pandas.DataFrame
        A DataFrame with columns ['type', 'path', 'radius']
        * the first row (df.iloc[0]) is soma. 
        * the first point of each path should be the branch point.
    """
    
    n = df_swc.n.values
    parent = df_swc.parent.values
    df_starting_points = df_swc[n - parent != 1]
    branchpoint_index = np.unique(df_starting_points.parent.values[1:])

    path_dict = {}
    type_dict = {}
    radius_dict = {}
    path_id = 1

    starting_points_pairs = get_consecutive_pairs_of_elements_from_list(df_starting_points.n.values)

    for s, e in starting_points_pairs: # s: start index, e: end index
        logging.debug('{}, {}'.format(s, e))
        if e is not None:
            b = branchpoint_index[np.logical_and(branchpoint_index > s, branchpoint_index < e)] # b: list, branch index
        else:
            b = branchpoint_index[branchpoint_index >= s]

        branchpoint_index_pairs = get_consecutive_pairs_of_elements_from_list(b, s, e)

        for bs, be in branchpoint_index_pairs:
            logging.debug('\t {}: {}, {}'.format(path_id, bs, be))
            
            path = df_swc.loc[bs:be][['x', 'y', 'z']].values
            path_type = df_swc.loc[bs:be][['type']].values
            path_radius = df_swc.loc[bs:be][['radius']].values
            
            if bs == s:
                parent_idx = df_swc.loc[[bs]].parent.values[0]
                if parent_idx != -1:

                    parent_coord = df_swc.loc[[parent_idx]][['x', 'y', 'z']].values
                    parent_radius = df_swc.loc[[parent_idx]][['radius']].values

                    if not (parent_coord == path).all(1).any():

                        path = np.vstack([parent_coord, path])
                        path_radius = np.vstack([parent_radius, path_radius])

            if be == e:
                path_dict[path_id] = path[:-1]
                radius_dict[path_id] = path_radius[:-1]
            else:
                path_dict[path_id] = path
                radius_dict[path_id] = path_radius

            type_dict[path_id] = max(path_type)[0]
            path_id += 1

    path_dict[0] = df_swc[df_swc.type == 1][['x', 'y', 'z']].values
    type_dict[0] = max(df_swc[df_swc.type == 1][['type']].values)[0]
    radius_dict[0] = df_swc[df_swc.type == 1][['radius']].values

    
    df_paths = pd.DataFrame([type_dict, path_dict, radius_dict]).T
    df_paths.columns = [['type', 'path', 'radius']]
    
    return df_paths

def swc_to_linestack(filepath, unit, voxelsize=None):

    """
    Convert SWC to Line Stack (from real length to voxel coordinates).
    """

    coords = pd.read_csv(filepath, comment='#', sep=' ', header=None)[[2,3,4]].as_matrix()
    
    if unit == 'pixel':
        coords = np.round(coords).astype(int)

    else: # unit == 'um'

        if voxelsize is None:
            logging.debug('  Not able to build linestack from float coordinates.')
            return None
        else:
            # coords = np.round(coords / voxelsize).astype(int) # not able to handle the soma-centered swc.
            logging.debug('  coords in um are converted back to pixel.')
            coords = coords - coords.min(0)
            coords = np.round(coords / voxelsize).astype(int)   
            # coords = np.round(coords / voxelsize)
            # coords = coords - coords.min(0)
            # coords = coords.astype(int)
            # logging.debug('{}'.format(coords.max(0)))

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
    # return (current_path[0] == soma).all(1).any()
    return (current_path == soma).all(1).any()

def find_connection(all_paths, soma, path_id, paths_to_ignore=[]):
    
    current_path = all_paths[path_id]
    
    
    if connect_to_soma(current_path, soma):

        connect_to = 0 
        connect_to_at = soma[0]
        
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
    
    logging.debug("Path {} connects to no other path. Try fix it.".format(path_id))
    connect_to = -99 
    connect_to_at = np.nan
    
    return connect_to, connect_to_at

def back2soma(df_paths, path_id):
    
    '''
    Given a path_id, find all the other paths which leads back to soma.
    '''
    
    path_id_original = path_id
    
    paths_to_soma = []
    
    counter = 0
        
    while df_paths.loc[path_id].connect_to != 0:
        if path_id in paths_to_soma:
            logging.debug("\tPath {} cannot trace back to soma: {}".format(path_id_original, paths_to_soma))
            break
        else:
            paths_to_soma.append(path_id)
            path_id = df_paths.loc[path_id].connect_to   

        if path_id == -99:
            break

    paths_to_soma.append(path_id)  

    return paths_to_soma

def update_df_paths(df_paths):

    """
    Loop through all paths, add fours columns ['connect_to', 'connect_to_at', 'connected_by', 'connected_by_at'] to df_paths.
    """
    
    logging.debug('  Start: Updating `df_paths` with connectivity debug.')

    soma = df_paths.loc[0]['path']
    
    all_paths = df_paths.loc[1:].path.to_dict()
    all_keys = list(all_paths.keys())

    connect_to_dict = {}
    connect_to_at_dict = {}

    for path_id in all_keys:
        connect_to_dict[path_id], connect_to_at_dict[path_id] = find_connection(all_paths, soma, path_id, paths_to_ignore=[])

    connect_to_dict.update({0: -1})
    connect_to_at_dict.update({0: []})
    
    df_paths['connect_to'] = pd.Series(connect_to_dict)
    df_paths['connect_to_at'] = pd.Series(connect_to_at_dict)

    connected_by_dict = {}
    connected_by_at_dict = {}

    for path_id in all_keys:

        connected_by_dict[path_id]    = df_paths[df_paths.connect_to == path_id].index.tolist()
        connected_by_at_dict[path_id] = df_paths[df_paths.connect_to == path_id].connect_to_at.tolist()

    connected_by_dict.update({0: list(df_paths[df_paths['connect_to'] == 0].index)})
    connected_by_at_dict.update({0: df_paths.loc[0].path[0]})

    df_paths['connected_by'] = pd.Series(connected_by_dict)
    df_paths['connected_by_at'] = pd.Series(connected_by_at_dict)

#     df_paths = add_soma_to_some_paths(df_paths, soma)
    
    logging.debug('  Finished.\n')

    return df_paths

def get_sorder(df_paths):
    
    df_paths['sorder'] = np.ones(len(df_paths)) * np.nan
    df_paths = df_paths.set_value(df_paths.connected_by.apply(len) == 0, 'sorder', 1)
    
    while np.isnan(df_paths.sorder).any():
    
        df_sub = df_paths[np.isnan(df_paths.sorder)]

        for row in df_sub.iterrows():

            path_id = row[0]
            connected_by = row[1]['connected_by']

            sorder0 = df_paths.loc[connected_by[0]].sorder
            sorder1 = df_paths.loc[connected_by[1]].sorder

            if np.isnan(sorder0) or np.isnan(sorder1):
                continue
            else:
                
                if sorder0 == sorder1:
                    df_paths.set_value(path_id, 'sorder', sorder0+1)
                else:
                    df_paths.set_value(path_id, 'sorder', np.max([sorder0, sorder1]))

    df_paths.sorder = df_paths['sorder'].astype(int)
                    
    return df_paths

def get_path_dendritic_length(path):
    return np.sum(np.sqrt(np.sum((path[1:] - path[:-1])**2, 1)))

def get_path_euclidean_length(path):
    return np.sqrt(np.sum((path[0] - path[-1]) ** 2))

def unique_row(a):
    
    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    _, idx = np.unique(b, return_index=True)
    
    unique_a = a[idx]
    
    return unique_a

def get_outer_terminals(all_terminals):
    
    from scipy.spatial import ConvexHull
    hull = ConvexHull(all_terminals[:,:2])
    outer_terminals_3d = all_terminals[hull.vertices]
#     outer_terminals_2d = all_terminals[:, :2][hull.vertices]
    outer_terminals_3d = np.vstack([outer_terminals_3d, outer_terminals_3d[0]])
    
    return outer_terminals_3d

def get_angle(v0, v1):
    c = np.dot(v0, v1) / np.linalg.norm(v0) / np.linalg.norm(v1)
    return np.arccos(np.clip(c, -1, 1)), np.degrees(np.arccos(np.clip(c, -1, 1)))

def get_node_vector(df_paths, path_id):
    s = df_paths.loc[path_id].path[0]
    # e = df_paths.loc[path_id].connected_by_at
    # if np.isnan(e).any():
    e = df_paths.loc[path_id].path[-1]
    v= e-s
    return v/np.linalg.norm(v)

def get_nearest_recorded_point_vector(df_paths, path_id):

    v = np.array([0,0,0])
    s = df_paths.loc[path_id].path[0]
    i = 0
    while (v == 0).all():
        i+=1
        e = df_paths.loc[path_id].path[i]
        v = e-s
        if i>5:break
    
    return v/np.linalg.norm(v)

def get_path_statistics(df_paths):
    
    logging.debug('  Start: Calculating path statistics (e.g. dendritic length, branch order...)')


    all_keys = df_paths.index[1:]
    
    dendritic_length_dict = {}
    euclidean_length_dict = {}
    back_to_soma_dict = {}
    corder_dict = {}
    
    for path_id in all_keys:
        
        path = df_paths.loc[path_id].path
        
        dendritic_length_dict[path_id] = get_path_dendritic_length(path)
        euclidean_length_dict[path_id] = get_path_euclidean_length(path)
        back_to_soma_dict[path_id] = back2soma(df_paths, path_id)
        corder_dict[path_id] = len(back_to_soma_dict[path_id])

    corder_dict.update({0: 0})

    df_paths['dendritic_length'] = pd.Series(dendritic_length_dict)
    df_paths['euclidean_length'] = pd.Series(euclidean_length_dict)
    df_paths['back_to_soma'] = pd.Series(back_to_soma_dict)
    df_paths['corder'] = pd.Series(corder_dict)

    logging.debug('  Finished. \n')
    
    return df_paths

def calculate_density(linestack, voxelsize):
    '''
    reference: A. Stepanyantsa & D.B. Chklovskiib (2005). Neurogeometry and potential synaptic connectivity. Trends in Neuroscience.
    '''
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

def get_average_angles(df_paths):

    nodal_angles_deg = {}
    nodal_angles_rad = {}
    
    local_angles_deg = {}
    local_angles_rad = {}
    
    n = 0
    for i in np.unique(df_paths.connect_to):

        if i == -1:
            continue

        path_ids = df_paths[df_paths.connect_to == i].index.tolist()
        
        # print(i, len(path_ids))
        if len(path_ids) == 2:
            v00 = get_node_vector(df_paths, path_ids[0])
            v01 = get_node_vector(df_paths, path_ids[1])
            nodal_angles_rad[n], nodal_angles_deg[n] = get_angle(v00,v01)

            v10 = get_nearest_recorded_point_vector(df_paths, path_ids[0])
            v11 = get_nearest_recorded_point_vector(df_paths, path_ids[1])
            local_angles_rad[n], local_angles_deg[n] = get_angle(v10,v11)

            n+=1
        else:
            continue

    average_nodal_angle_deg = np.mean(list(nodal_angles_deg.values()))
    average_nodal_angle_rad = np.mean(list(nodal_angles_rad.values()))

    average_local_angle_deg = np.mean(list(local_angles_deg.values()))
    average_local_angle_rad = np.mean(list(local_angles_rad.values()))

    return average_nodal_angle_deg, average_nodal_angle_rad, average_local_angle_deg, average_local_angle_rad


def plot_skeleten(ax, df_paths, soma, axis0, axis1, order_type, lims):

    if order_type == 'c':
        colors = plt.cm.viridis.colors
        colors_idx = np.linspace(0, 255, max(df_paths.corder)+1).astype(int)
    elif order_type == 's':
        colors = plt.cm.viridis_r.colors
        colors_idx = np.linspace(0, 255, max(df_paths.sorder)+1).astype(int)
        
    ax.scatter(soma[axis0], soma[axis1], s=280, color='grey')
    for row in df_paths.iterrows():

        path_id = row[0]

        path = row[1]['path']
        bpt = path[0]
        if order_type == 'c':
            order = row[1]['corder']
            ax.plot(path[:, axis0], path[:, axis1], color=colors[colors_idx[int(order)]])
            ax.scatter(bpt[axis0], bpt[axis1], color=colors[colors_idx[int(order)]], zorder=1)

        elif order_type == 's':
            order = row[1]['sorder']
            ax.plot(path[:, axis0], path[:, axis1], color=colors[colors_idx[int(order)-1]])
            ax.scatter(bpt[axis0], bpt[axis1], color=colors[colors_idx[int(order)-1]], zorder=1)
    
    xylims, zlim = lims

    if axis0 == 2 and axis1 == 0: # ax2
        ax.set_xlim(zlim[0], zlim[1])
        ax.set_ylim(xylims[0], xylims[1])   

    elif axis0 == 1 and axis1 == 2: # ax3
        ax.set_xlim(xylims[0], xylims[1])
        ax.set_ylim(zlim[0], zlim[1])
     

    elif axis0 == 1 and axis1 == 0: # ax1
        ax.set_xlim(xylims[0], xylims[1])
        ax.set_ylim(xylims[0], xylims[1])

    ax.axis('off')

def find_lims(df_paths):
    points = np.vstack(df_paths.path)
    maxlims = np.max(points, 0)
    minlims = np.min(points, 0)
    xylims = np.hstack([maxlims[:2], minlims[:2]])
    zlims = np.hstack([maxlims[2], minlims[2]])
    if (xylims >= 0).all():
        xylims = np.array([0, max(xylims).astype(int) + 30])
    else:
        xylims = np.array([-max(abs(xylims)).astype(int) - 30, max(abs(xylims)).astype(int) + 30])

    if (zlims >= 0).all():
        zlims = np.array([-10, max(zlims).astype(int) + 30])
    else:
        zlims = np.array([-max(abs(zlims)).astype(int) - 30, max(abs(zlims)).astype(int) + 30])
        
    return xylims, zlims

def get_path_on_stack(df_paths, voxelsize, coordinate_padding):

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


def print_summary(summary, unit):
    
    """
    Print out summary statistics of the cell.

    Parameters
    ----------
    summary: dict
        a nested dict contains summary of the cell.

    unit: str
        the unit of swc file. Either 'um' or 'pixel'.
    """

    import logging
    
    num_dendritic_segments = summary['general']['number_of_dendritic_segments']
    num_branchpoints = summary['general']['number_of_branch_points']
    num_irreducible_nodes = summary['general']['number_of_irreducible_nodes']
    max_branch_order = summary['general']['max_branch_order']
    max_strahler_order = summary['general'] ['max_strahler_order']
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
    logging.info('    Max Strahler order: {}\n\n'.format(max_strahler_order))

    logging.info('  # Angle \n')
    logging.info('    Average nodal angle in degree: {:.3f}'.format(average_nodal_angle_deg))
    logging.info('    Average nodal angle in radian: {:.3f} \n'.format(average_nodal_angle_rad))
    logging.info('    Average local angle in degree: {:.3f}'.format(average_local_angle_deg))
    logging.info('    Average local angle in radian: {:.3f} \n'.format(average_local_angle_rad))

    logging.info('  # Average tortuosity: {:.3f}\n'.format(average_tortuosity))

    logging.info('  # Dendritic length ({})\n'.format(unit))
    # logging.info('  ## Dendritic length\n')
    logging.info('    Sum: {:.3f}'.format(dendritic_sum))
    logging.info('    Mean: {:.3f}'.format(dendritic_mean))
    logging.info('    Median: {:.3f}'.format(dendritic_median))
    logging.info('    Min: {:.3f}'.format(dendritic_min))
    logging.info('    Max: {:.3f}\n'.format(dendritic_max))

    logging.info('  # Euclidean length ({})\n'.format(unit))
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
        
        logging.info('  # Density related ({})\n'.format(unit))
        logging.info('    Asymmetry: {:.3f}'.format(asymmetry))
        logging.info('    Outer Radius: {:.3f}'.format(outer_radius))
        logging.info('    Typical Radius : {:.3f}\n'.format(typical_radius))
        logging.info('    Dendritic Area: {:.3f} ×10\u00b3 um\u00b2\n\n'.format(dendritic_area))
