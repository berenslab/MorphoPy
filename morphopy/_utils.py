import numpy as np
import pandas as pd
from copy import deepcopy
import logging
import matplotlib.pyplot as plt

def read_swc(filepath, unit, voxelsize):
    
    logging.info('  Start: Reading `.swc` file from \n\t\t{}\n\t\tinto Pandas DataFrame.'.format(filepath))
    df = pd.read_csv(filepath, comment='#', sep=' ', header=None)
    df = df.dropna(axis=1)
    df.columns = ['ID', 'Type', 'x', 'y', 'z', 'Radius', 'PID']
    df.index = df.ID.as_matrix()

    if unit=='pixel' and voxelsize is not None:
        df[['x', 'y', 'z']] = df[['x', 'y', 'z']] * voxelsize
        unit_of_df = 'um'
    else:
        unit_of_df = unit

    logging.info('  Finished.\n')
    return df, unit_of_df

def swc2linestack(filepath, unit, imagesize, voxelsize=None):

    coords = pd.read_csv(filepath, comment='#', sep=' ', header=None)[[2,3,4]].as_matrix()
    
    if unit == 'pixel':
        if imagesize is None:
            logging.info('  No `imagesize` is provided.')
            return None
        else:
            coords = np.round(coords).astype(int)

    else: # unit == 'um'
        if imagesize is None:
            logging.info('  No `imagesize` is provided.')
            return None
        else:
            if voxelsize is None:
                logging.info('  Not able to build linestack from float coordinates.')
                return None
            else:
                # coords = np.round(coords / voxelsize).astype(int) # not able to handle the soma-centered swc.
                logging.info('  coords in um are converted back to pixel.')
                coords = np.round(coords / voxelsize)
                coords = coords - coords.min(0) + 10
                coords = coords.astype(int)
                logging.info('{}'.format(coords.max(0)))

    logging.info('  Start: Creating linestack...')
    linestack = np.zeros(imagesize)
    for c in coords:
        linestack[tuple(c)] = 1
    logging.info('  Finished.\n')
    
    return linestack   


def get_soma(df_swc):
    
    soma = df_swc[df_swc.Type == 1]
    soma_pos = soma[['x','y','z']].as_matrix()
    neurites = df_swc[df_swc.Type != 1]
       
    return soma_pos, neurites

def connect_to_soma(current_path, soma):
    return (current_path[0] == soma).all(1).any()


def get_df_paths(df_swc):

    logging.info('  Start: Creating `df_paths` from `df_swc`.')
    
    dict_path = {}
    path_id = 0
    path = []
    
    soma, neurites = get_soma(df_swc)

    for counter, row in enumerate(neurites.iterrows()):

        idx = row[0]
        ID = row[1]['ID']
        PID = row[1]['PID']
        current_point = row[1][['x', 'y', 'z']].as_matrix()

        if counter == 0:
            path.append(soma[0])

        if ID - 1 == PID:
            path.append(current_point)
        else:
            path_id += 1
            path = []
            parent_point = df_swc.loc[np.where(df_swc.ID == PID)[0]+1][['x','y','z']].as_matrix()
            path.append(parent_point)
            path.append(current_point)

        dict_path[path_id] = np.vstack(path)

    df_paths = pd.DataFrame()
    df_paths['path'] = pd.Series(dict_path)
    
    logging.info('  Finished.\n')

    return df_paths

def find_connection(all_paths, soma, path_id, paths_to_ignore=[]):
    
    current_path = all_paths[path_id]
    
    
    if connect_to_soma(current_path, soma):

        connect_to = -1 
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
    
    logging.info("Path {} connects to no other path. Try fix it.".format(path_id))
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
        
    while df_paths.loc[path_id].connect_to != -1:

        if path_id in paths_to_soma:
            logging.info("\tPath {} cannot trace back to soma: {}".format(path_id_original, paths_to_soma))
            break
        else:
            paths_to_soma.append(path_id)
            path_id = df_paths.loc[path_id].connect_to   

        if path_id == -99:
            break
    
    paths_to_soma.append(path_id)  

    return paths_to_soma

def update_df_paths(df_paths, soma):
    
    logging.info('  Start: Updating `df_paths` with connectivity info.')

    all_paths = df_paths.path.to_dict()
    all_keys = list(all_paths.keys())

    connect_to_dict = {}
    connect_to_at_dict = {}

    for path_id in all_keys:
        connect_to_dict[path_id], connect_to_at_dict[path_id] = find_connection(all_paths, soma, path_id, paths_to_ignore=[])

    df_paths['connect_to'] = pd.Series(connect_to_dict)
    df_paths['connect_to_at'] = pd.Series(connect_to_at_dict)

    connected_by_dict = {}
    connected_by_at_dict = {}

    for path_id in all_keys:

        connected_by_dict[path_id]    = df_paths[df_paths.connect_to == path_id].index.tolist()
        connected_by_at_dict[path_id] = df_paths[df_paths.connect_to == path_id].connect_to_at.tolist()

    df_paths['connected_by'] = pd.Series(connected_by_dict)
    df_paths['connected_by_at'] = pd.Series(connected_by_at_dict)

    df_paths = add_soma_to_some_paths(df_paths)
    
    logging.info('  Finished.\n')

    return df_paths

def detect_messiness(df_paths, threshold):
    
    logging.info('  Start: Finding if there are disconnected paths. ')

    paths_to_fix = []
    
    for i, path in enumerate(df_paths.path):

        tip_distance = np.sqrt(np.sum((path[0] - path[1]) ** 2))

        if tip_distance > threshold:
            paths_to_fix.append(i+1)
    
    if len(paths_to_fix) > 0:
        logging.info('  Paths {} needs to be fixed'.format(paths_to_fix))
    else:
        logging.info('  All paths are fine. ')

    logging.info('  Finished.\n')
    
    return paths_to_fix

def get_distance_tip_to_path(point, path):
    
    distance_tip_to_all_points_in_path = np.sqrt(((point - path)**2).sum(1))
    
    closest_point = np.argmin(distance_tip_to_all_points_in_path)
    closest_distance = distance_tip_to_all_points_in_path[closest_point]
    
    return [closest_point, closest_distance]

def find_closest_paths(all_paths, current_path, path_id):
    
    closest_paths = []
    
    sub_paths = deepcopy(all_paths)
    sub_paths.pop(path_id)
    
    for key in sub_paths.keys():
        
        target_path = sub_paths[key]
        
        closest_paths.append([key] + get_distance_tip_to_path(current_path[0], target_path))
        
    closest_paths = np.vstack(closest_paths) 
    
    connect_to = closest_paths[:, 0][np.argmin(closest_paths[:, 2])]
    connect_to_at_loc = closest_paths[:, 1][np.argmin(closest_paths[:, 2])]
    connect_to_at = all_paths[connect_to][int(connect_to_at_loc)]
    
    return connect_to, connect_to_at

def clean_messiness(df_paths, paths_to_fix):
    
    if len(paths_to_fix)>0:
        logging.info('  Start: Fixing disconnected paths...')

        df_paths_fixed = deepcopy(df_paths)
        
        all_paths = df_paths.path.to_dict()

        connect_to_dict = {}
        connect_to_at_dict = {} 
        
        for path_id in paths_to_fix:
            
            current_path = all_paths[path_id][1:]
                
            connect_to, connect_to_at = find_closest_paths(all_paths, current_path, path_id)
         
            current_path = np.vstack([connect_to_at, current_path])

            df_paths_fixed.set_value(path_id, 'path', current_path)
            df_paths_fixed.set_value(path_id, 'connect_to', connect_to)
            df_paths_fixed.set_value(path_id, 'connect_to_at', connect_to_at)
            
            logging.info("Fixed Path {} connections.".format(path_id))

        logging.info('  Finished.\n')
        
        return df_paths_fixed
    else:
        return df_paths

def get_bpt_loc_pairs(bpt_loc):
    
    if len(bpt_loc) == 1:
        return None
    else:
        bpt_loc_pairs = [[bpt_loc[i-1], bpt_loc[i]] for i in range(1,len(bpt_loc))]
        return np.vstack(bpt_loc_pairs)

def add_soma_to_some_paths(df_paths):
    
    paths_connect_to_soma = df_paths[df_paths.connect_to == -1]
    soma = paths_connect_to_soma.iloc[0].connect_to_at    
    for row in paths_connect_to_soma.iterrows():

        path_id = row[0]
        path = row[1]['path']

        if (soma == path).all(1).any():
            continue
        else:
            path = np.vstack([soma, path])

        df_paths = df_paths.set_value(path_id, 'path', path)

    return df_paths 

def get_better_paths(df_paths, soma):

    logging.info('  Start: breaking up paths at branch points...')
        
    new_paths_dict = {}
    path_id_old = {}

    for row in df_paths.iterrows():

        current_number_of_paths= len(new_paths_dict)
        path_id = row[0]
        path = row[1]['path']

        bpts = row[1]['connected_by_at']

        if len(bpts) != 0:
            # bpts = np.vstack([path[0], bpts])
            if not (path[0] == bpts).all(1).any():
                bpts = np.vstack([path[0], bpts])
            else:
                bpts = np.vstack(bpts)
        else:
            bpts = path[0].reshape(1,3)

        logging.debug('  Processing Path {}'.format(path_id))

        if (bpts[-1] == path[-1]).all():
            logging.debug('\tDelete the last branchpoint,\n\tbecause it is the last point of the path.')
            bpts = np.delete(bpts, -1, axis=0)


        logging.debug('\t# bpts {}'.format(len(bpts)))

        bpt_loc = []
        for i in range(len(bpts)):

            loc_tuple = np.where((bpts[i] == path).all(1))
            if len(loc_tuple[0]) > 0:
                bpt_loc.append(loc_tuple[0][0])
        bpt_loc = np.sort(bpt_loc)

        bpt_loc_pairs = get_bpt_loc_pairs(bpt_loc)

        if bpt_loc_pairs is not None:

            num_skip = 0
            for key in range(len(bpt_loc_pairs)):

                new_path = path[bpt_loc_pairs[key][0]:bpt_loc_pairs[key][1]+1]
                if len(new_path)>1:
                    logging.debug('\tadding new path {}'.format(key+current_number_of_paths-num_skip))
                    new_paths_dict[key+current_number_of_paths-num_skip] = new_path
                else:
                    logging.debug('\tCurrent new path {} is just one point. Skipped.'.format(key+current_number_of_paths))
                    num_skip += 1

            new_path = path[bpt_loc_pairs[key][1]:]
            if len(new_path) > 1:
                logging.debug('\tadding new path {}\n'.format(key+current_number_of_paths-num_skip+1))
                new_paths_dict[key+current_number_of_paths-num_skip+1] = new_path
            else:
                logging.debug('\tCurrent new path {} is just one point. Skipped.'.format(key+current_number_of_paths))

        else:
            logging.debug('\tadding new path {}\n'.format(current_number_of_paths))
            new_paths_dict[current_number_of_paths] = path
            
    
    df_paths_updated = pd.DataFrame()
    df_paths_updated['path'] = pd.Series(new_paths_dict)
    df_paths_updated = update_df_paths(df_paths_updated, soma)

    logging.info('  Finished. \t')
    
    return df_paths_updated

def cleanup_better_paths(df_paths):

    logging.info('  Start: Checking and cleaning up messy paths (e.g. connecting the nonbifucated paths.)..')

    df = deepcopy(df_paths)
    
    df_paths_to_clean = df[df.connected_by.apply(len) == 1]
    
    path_ids_to_clean = list(df_paths_to_clean.index)
    
    for path_id in path_ids_to_clean:
        
        path_id_head = path_id
        path_id_tail = df.loc[path_id].connected_by[0]
        
        logging.info("Connecting Path {} and {}".format(path_id_head, path_id_tail))

        path_head = df.loc[path_id].path
        path_tail = df.loc[path_id_tail].path
        
        path = np.vstack([path_head, path_tail])
        
        connect_to = df.loc[path_id_head].connect_to
        connect_to_at = df.loc[path_id_head].connect_to_at
        
        df.set_value(path_id_tail, 'path', path)
        df.set_value(path_id_tail, 'connect_to', connect_to)
        df.set_value(path_id_tail, 'connect_to_at', connect_to_at)
        
        path_id_connect_to = df.loc[path_id].connect_to
        
        if path_id_connect_to != -1:
            connected_by = np.array(df.loc[path_id_connect_to].connected_by)
            connected_by[np.where(connected_by == path_id_head)[0]] = path_id_tail
            df.set_value(path_id_connect_to, 'connected_by', connected_by)

        logging.info('  Path {} is removed.'.format(path_id_head))
        df.drop(path_id_head, inplace=True)

    logging.info('  Finished. \n')

    return df

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
    
    logging.info('  Start: Calculating path statistics (e.g. dendritic length, branch order...)')

    all_keys = df_paths.index
    
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

    df_paths['dendritic_length'] = pd.Series(dendritic_length_dict)
    df_paths['euclidean_length'] = pd.Series(euclidean_length_dict)
    df_paths['back_to_soma'] = pd.Series(back_to_soma_dict)
    df_paths['corder'] = pd.Series(corder_dict)

    logging.info('  Finished. \n')
    
    return df_paths

def calculate_density(linestack, voxelsize):
    '''
    reference: A. Stepanyantsa & D.B. Chklovskiib (2005). Neurogeometry and potential synaptic connectivity. Trends in Neuroscience.
    '''
    import scipy.ndimage as ndimage

    logging.info('  Start: Calculating dendritic density...')

    smoothed_layer_stack = []
    for i in range(linestack.shape[2]):

        layer = linestack[:,:,i]
        smoothed_layer_stack.append(ndimage.gaussian_filter(layer, sigma=25/voxelsize[0], order=0))

    density_stack = np.dstack(smoothed_layer_stack)
    center_of_mass = np.array(ndimage.measurements.center_of_mass(density_stack.sum(2)))

    logging.info('  Finished. \n')

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
        
    ax.scatter(soma[0][axis0], soma[0][axis1], s=280, color='grey')
    for row in df_paths.iterrows():

        path_id = row[0]

        path = row[1]['path']
        bpt = path[0]
        if order_type == 'c':
            order = row[1]['corder']
            ax.plot(path[:, axis0], path[:, axis1], color=colors[colors_idx[order]])
            ax.scatter(bpt[axis0], bpt[axis1], color=colors[colors_idx[order]], zorder=1)

        elif order_type == 's':
            order = row[1]['sorder']
            ax.plot(path[:, axis0], path[:, axis1], color=colors[colors_idx[order-1]])
            ax.scatter(bpt[axis0], bpt[axis1], color=colors[colors_idx[order-1]], zorder=1)
    
    # maxlims, minlims = lims
    # ax.set_xlim(0, lims[axis0])
    # ax.set_ylim(0, lims[axis1])
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