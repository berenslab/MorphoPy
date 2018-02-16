import numpy as np
import pandas as pd
import logging

# from ._morph import *

def swc_to_df(filepath):
    df_swc =  pd.read_csv(filepath, delim_whitespace=True, comment='#',
                          names=['n', 'type', 'x', 'y', 'z', 'radius', 'parent'], index_col=False)
    df_swc.index = df_swc.n.as_matrix()
    return df_swc

def read_swc(filepath):

    logging.info('  Reading {}\n'.format(filepath))
    
    swc = swc_to_df(filepath)
    # raw data
    n = np.array(swc['n'].tolist())
    pos = np.array([swc['x'].tolist(), swc['y'].tolist(), swc['z'].tolist()]).T
    radius = np.array(swc['radius'].tolist())
    t = np.array(swc['type'].tolist())
    pid = np.array(swc['parent'].tolist())
    t[pid == -1] = 1 # if soma is missing, the first point is soma

    e = np.vstack(list(zip(pid[1:], n[1:])))
    e = remove_duplicate(pos, e)
    
    
    soma_loc = 1
    
    return {'n': n,
            'pos':pos,
            'radius':radius,
            't': t,
            'e': e,
            'soma_loc': soma_loc
           }

def read_imx(filepath):

    import numpy as np
    import xml.etree.ElementTree as ET
    
    from bs4 import BeautifulSoup
    
    logging.info('  Reading {}\n'.format(filepath))
    
    with open(filepath, 'rb') as f:
        xml_binary = f.read().decode('utf-8')
    
    soup = BeautifulSoup(xml_binary, 'lxml')
    
    vertex = soup.find('mfilamentgraphvertex').text
    vertex = vertex.replace(')', '')
    vertex = vertex.split('(')
    vertex = [i.strip() for i in vertex]
    vertex = [i for i in vertex if i != '']
    
    edges = soup.find('mfilamentgraphedge').text
    edges = edges.replace(')', '')
    edges = edges.split('(')
    edges = [i.strip() for i in edges]
    edges = [i for i in edges if i != '']
    
    vortex = np.array([list(map(float, sub.split(','))) for sub in vertex])
    # node position/coordinates(x,y,z)
    
    pos = vortex[:, :3] # n->pos    
    e = np.array([list(map(int, sub.split(','))) for sub in edges]) + 1
    e = remove_duplicate(pos, e)
    
    # node index
    n = np.arange(1, len(pos)+1) # n is the node index
    
    # radius
    radius = vortex[:, 3]
    
    #type
    t = np.ones(len(radius)) * 3
    soma_loc = int(soup.find('mfilamentgraph').attrs['mrootvertex']) + 1
    t[soma_loc-1] = 1
        
    return {'n': n,
            'pos':pos,
            'radius':radius,
            't': t,
            'e': e,
            'soma_loc': soma_loc
           }

def remove_duplicate(pos, e):

    logging.info('  Checking duplicate points in data.')
    
    pos_uniq, pos_count = np.unique(pos, return_counts=True, axis=0)

    logging.info('\t{} duplicate points are found.\n'.format(len(pos_uniq[pos_count > 1])))
    for duplicate_point in pos_uniq[pos_count > 1]:
        a, *rest = np.where((duplicate_point == pos).all(1))[0]+1
        for dupl in rest:
            for i in [i for i in np.where(e[:, 0] == dupl)[0]]:
                e[i][0] = a
            for i in [i for i in np.where(e[:, 1] == dupl)[0]]:
                e[i][1] = a
    
    
    e[(e[:, 1] - e[:, 0]) < 0] = e[(e[:, 1] - e[:, 0]) < 0][:, ::-1]
    
    e = np.delete(e, np.where(e[:, 0] == e[:, 1])[0], axis=0)
    
    e, index = np.unique(e, axis=0, return_index=True)
    
    return e[index.argsort()]

def get_edge_dict(n, e, soma_loc):

    logging.info("  Creating path arrays from voxel data points.\n")

    def _point_already_in_dict(point, edge_dict):
    
        for path_id, points_list in edge_dict.items():
            if point in points_list:
                return True
        else:
            return False

    edge_dict = {}
    path_id = 1
    
    branch_loc = [i for i in n if len(np.where(e[:, 0] == i)[0]) >=2 ]
    if 1 not in branch_loc:
        branch_loc = [1] + branch_loc
    
    for bpt in branch_loc:
        
        bpt_locs_on_e = np.where(e[:, 0] == bpt)[0]
        
        for edge in e[bpt_locs_on_e]:
            current_point = edge[0]
            next_point = edge[1]
            
            a_list = [current_point]
            
            if _point_already_in_dict(next_point, edge_dict):
                logging.debug('\tPoint {} is already used by other paths. Skip.'.format(next_point))
            else:   
                a_list.append(next_point)# a list for holding the index of point of one paths
            
            next_point_locs_on_e = np.where(e[:, 0] == next_point)[0]
            
            while len(next_point_locs_on_e) == 1:
                
                next_point = e[next_point_locs_on_e[0]][1]
                if _point_already_in_dict(next_point, edge_dict):
                    logging.debug('\tPoint {} is already used by other paths. Skip.'.format(next_point))
                else:   
                    a_list.append(next_point)# a list for holding the index of point of one paths
                next_point_locs_on_e = np.where(e[:, 0] == next_point)[0]
                
            if len(a_list) < 2:
                logging.debug('\t\tCurrent path has fewer than two points, skipped.')
                continue
            
            edge_dict[path_id] = np.array(a_list)
            path_id += 1
            
    if soma_loc not in branch_loc:
        paths_soma_on = [key for key, value in edge_dict.items() if soma_loc in value]


        for path_id in paths_soma_on:
            path = edge_dict[path_id]
            breakup_point = np.where(edge_dict[path_id] == soma_loc)[0] 
            path_0 = path[:breakup_point[0]+1][::-1]
            path_1 = path[breakup_point[0]:]
            edge_dict[path_id] = path_0
            edge_dict[len(edge_dict)] = path_1

    logging.debug('  Done.\n')
            
    return edge_dict

def get_path_dict(pos, radius, t, edge_dict, soma_loc):
    
    path_dict = {}
    radius_dict = {}
    type_dict = {}
    all_keys = edge_dict.keys()
    for key in all_keys:
        
        path_dict[key] = pos[edge_dict[key]-1]
        radius_dict[key] = radius[edge_dict[key]-1]
        type_dict[key] = np.unique(t[edge_dict[key]-1][1:])[0]
    
    path_dict.update({0: pos[soma_loc-1].reshape(1,3)})
    radius_dict.update({0:[radius[soma_loc-1]]})
    type_dict.update({0: 1})
    
    df_paths = pd.DataFrame()
    df_paths['type'] = pd.Series(type_dict)
    df_paths['path'] = pd.Series(path_dict)
    df_paths['radius'] = pd.Series(radius_dict)
    
    return df_paths

def sort_path_direction(df_paths):
    
    df_paths = df_paths.copy()
    soma = df_paths.loc[0].path.flatten()    
    
        
    df_paths['connect_to'] = np.nan
    df_paths['connect_to_at'] = ''
    df_paths['connect_to_at'] = df_paths['connect_to_at'].apply(np.array)

    path_ids_head = df_paths[df_paths.path.apply(lambda x: (x[0] == soma).all())].index

    if len(path_ids_head) > 0:
        df_paths.loc[path_ids_head, 'connect_to'] = -1
        df_paths.loc[path_ids_head, 'connect_to_at'] = pd.Series({path_id:soma for path_id in path_ids_head})

    path_ids_tail = df_paths[df_paths.path.apply(lambda x: (x[-1] == soma).all())].index
    if len(path_ids_tail) > 0:
        df_paths.loc[path_ids_tail, 'path'] = df_paths.loc[path_ids_tail].path.apply(lambda x: x[::-1])
        df_paths.loc[path_ids_tail, 'connect_to'] = -1
        df_paths.loc[path_ids_tail, 'connect_to_at'] = pd.Series({path_id:soma for path_id in path_ids_tail})

    new_target_paths = list(df_paths[~np.isnan(df_paths.connect_to)].index) # seed the first round of paths to check

    logging.info('  Checking path connection.')
    logging.info('\tTotal num of paths to check: {}\n'.format(len(df_paths)))    
    logging.debug('\tNumber of paths connected to soma: {}'.format(len(new_target_paths)-1))
    
    while np.count_nonzero(~np.isnan(df_paths.connect_to)) != len(df_paths):
        
        all_checked_paths = list(df_paths[~np.isnan(df_paths.connect_to)].index)
        num_check_paths_before = len(all_checked_paths)

        target_paths = new_target_paths
        new_target_paths = [] # empty the list to hold new target paths for next round
        logging.debug("  Current target paths {}".format(target_paths))
        logging.debug("  Next target paths: {}".format(new_target_paths))
        
        for target_path_id in target_paths:
            
            if target_path_id == 0: continue

            target_path = df_paths.loc[target_path_id].path

            path_ids_head = df_paths[df_paths.path.apply(lambda x: (x[0] == target_path[-1]).all())].index.tolist()
            path_ids_head = [i for i in path_ids_head if i not in all_checked_paths]
            logging.debug("\t{} paths connect to target path by their head: {}".format(len(path_ids_head), path_ids_head))
            if len(path_ids_head) > 0:
                df_paths.loc[path_ids_head, 'connect_to'] = target_path_id
                df_paths.loc[path_ids_head, 'connect_to_at'] = pd.Series({path_id:target_path[-1] for path_id in path_ids_head})
                new_target_paths = new_target_paths + path_ids_head
                logging.debug("  \tadding {} to next target paths".format(path_ids_head))
                
            path_ids_tail = df_paths[df_paths.path.apply(lambda x: (x[-1] == target_path[-1]).all())].index.tolist()
            path_ids_tail = [i for i in path_ids_tail if i not in all_checked_paths]
            logging.debug("\t{} paths connect to target path by their tail: {}".format(len(path_ids_tail), path_ids_tail))

            if len(path_ids_tail) > 0:
                df_paths.loc[path_ids_tail, 'path'] = df_paths.loc[path_ids_tail].path.apply(lambda x: x[::-1])
                df_paths.loc[path_ids_tail, 'connect_to'] = target_path_id
                df_paths.loc[path_ids_tail, 'connect_to_at'] = pd.Series({path_id:target_path[-1] for path_id in path_ids_tail})
                new_target_paths = new_target_paths + path_ids_tail
                logging.debug("  \tadding {} to next target paths".format(path_ids_tail))
                
        logging.debug('  Next target paths is: {}'.format(new_target_paths))
                
        num_check_paths_after = len(list(df_paths[~np.isnan(df_paths.connect_to)].index))
        logging.debug("\tNumber of paths checked: {}".format(num_check_paths_after))
        
        if num_check_paths_before == num_check_paths_after:
            num_disconneted = len(df_paths) - num_check_paths_after
            logging.info('\tNumber of disconnected path(s): {}'.format(num_disconneted))
            break

    
    df_paths_drop = df_paths[np.isnan(df_paths.connect_to)] 
    df_paths = df_paths.drop(df_paths[np.isnan(df_paths.connect_to)].index)

    return df_paths, df_paths_drop

def get_paths_nearest_to_tree(df_paths, df_drops, num_all_paths):
    
    """
    Get paths in df_drops which stay nearest to the connected paths. 
    
    Paremeters
    ----------
    df_paths: 
        DataFrame holding all connected paths
    
    df_drops:
        DataFrame holding all disconnected paths
        
    Return 
    ------
    res: dict
        {'p0': {'path': path,
                'path_id': path_id_drop,
                'radius': radius,
                'type': t},
         'target': {'path': path_tree,
                     'path_id': path_id_tree,
                     'radius': radius_tree,
                     'type': type_tree
                     } }

    """
    
    res = []
    for row in df_drops.iterrows():
        
        path_id_drop = row[0]
        path = row[1].path
        
        path_id_tree = df_paths.path.apply(lambda x: np.sqrt(((x[-1]-path)**2).sum(1).min())).argmin()
        distance_arr = np.sqrt(np.sum((path - df_paths.loc[path_id_tree].path[-1]) ** 2, 1))
        
        loc_nearest = distance_arr.argmin()
        dist_nearest = distance_arr.min()        
        
        res.append((path_id_drop, path_id_tree, loc_nearest, dist_nearest))
    
    res = np.array(res)            
    nearest_paths_data = res[np.where(res[:, 3] == res[:, 3].min())[0]]
    
    ### get all paths data
    
    num_paths_to_tree = len(nearest_paths_data)
    
    path_id_tree = nearest_paths_data[0][1]
    path_tree = df_paths.loc[path_id_tree].path
    radius_tree = df_paths.loc[path_id_tree].radius
    type_tree = df_paths.loc[path_id_tree].type
    
    loc_nearest = int(nearest_paths_data[0][2])
    distance_between = nearest_paths_data[0][3]
    
    res = {}

    for i, datum in enumerate(nearest_paths_data):
        
        path_id_drop = datum[0]
        path = df_drops.loc[path_id_drop].path
        radius = df_drops.loc[path_id_drop].radius
        
        t = df_drops.loc[path_id_drop].type
        
        if loc_nearest == 0:
            path = path 
            point_connect = path[0]
            radius_connect = radius[0]
            
            res['p{}'.format(i)] = {'path': path,
                                    'path_id': path_id_drop,
                                    'radius': radius,
                                    'type': t}
        elif loc_nearest == len(path)-1:
            path = path[::-1]
            point_connect = path[0]
            radius_connect = radius[0]
            
            res['p{}'.format(i)] = {'path': path,
                                    'path_id': path_id_drop,
                                    'radius': radius,
                                    'type': t}
        else:
            p0 = path[:loc_nearest+1][::-1] 
            p1 = path[loc_nearest:]
            
            r0 = radius[:loc_nearest+1]
            r1 = radius[loc_nearest:]
            
            point_connect = p0[0]
            radius_connect = r0[0]
            
            if len(p0) > 1:
                res['p0'.format(i)] = {'path': p0,
                                        'path_id': path_id_drop,
                                        'radius': r0,
                                        'type': t}
            if len(p1) > 1:
                res['p1'.format(i)] = {'path': p1,
                                    # 'path_id': len(df_paths) + len(df_drops),
                                    'path_id': num_all_paths,
                                    'radius': r1,
                                    'type': t}
    
    if distance_between > 0:
        logging.debug("\tPath {} appended point {} at the end".format(path_id_tree, point_connect))
        path_tree = np.vstack([path_tree, point_connect])
        radius_tree = np.hstack([radius_tree, radius_connect])
    
    res['target'] = {'path': path_tree,
                     'path_id': path_id_tree,
                     'radius': radius_tree,
                     'type': type_tree
                     }    
    
    return res

def reconnect_dropped_paths(df_paths, df_drops):
    
    if len(df_drops) > 0:
        logging.info('  Connecting disconnected paths.\n')
    else:
        logging.info('  No disconnected paths.\n')

    df_paths = df_paths.copy()
    df_drops = df_drops.copy()
    num_dropped_paths = len(df_drops)
    
    while len(df_drops) > 0:
        logging.debug('Lenght of df_paths: {}; df_drops: {}'.format(len(df_paths), len(df_drops)))
        num_all_paths = num_dropped_paths + len(df_paths)
        paths_data = get_paths_nearest_to_tree(df_paths, df_drops, num_all_paths)

        target = paths_data.pop('target')
        path_id_tree = target['path_id']
        path_tree = target['path']
        radius_tree = target['radius']
        tail_path_tree = path_tree[-1]
        
        df_paths.at[int(path_id_tree), 'path'] = path_tree
        df_paths.at[int(path_id_tree), 'radius'] = radius_tree

        if len(paths_data) == 1:

            p = paths_data.pop('p0')

            path_id_drop = p['path_id']
            path_drop = p['path']
            radius_drop = p['radius']
            t = p['type']
            logging.debug('  Append {} to {}'.format(path_id_drop, path_id_tree))
            df_paths.at[int(path_id_tree), 'path'] = np.vstack([path_tree[:-1], path_drop]) 
            df_paths.at[int(path_id_tree), 'radius'] = np.hstack([radius_tree[:-1], radius_drop] ) 
            df_drops.drop(path_id_drop, inplace=True)

        else:
            for key, values in paths_data.items():

                p = paths_data[key]

                path_id_drop = p['path_id']
                path_drop = p['path']
                logging.debug('  Add row: {}, len: {}'.format(int(path_id_drop), len(path_drop)))
                radius_drop = p['radius']
                t = p['type']
                df_paths.loc[int(path_id_drop)] = [t, path_drop, radius_drop, path_id_tree, tail_path_tree]            

                try:
                    logging.debug('\tdrop {}'.format(int(path_id_drop)))
                    df_drops.drop(path_id_drop, inplace=True)
                except:
                    pass

            logging.info('Lenght of df_drops: {}'.format(len(df_drops)))
    
    return df_paths.sort_index()

def find_connection(df_paths):
    
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
    
    logging.debug('  Done.\n')
    
    return df_paths

def write_swc(df_paths):

    path_checked = []
    swc_arr = []
    
    list_back_to_soma = (df_paths.sort_values(['connect_to']).back_to_soma).tolist()
    for i, back_to_soma in enumerate(list_back_to_soma):

        for path_id in back_to_soma[::-1]:

            if path_id in path_checked: 
                continue

            path_s = df_paths.loc[path_id]
            path = path_s['path']
            path_radius = path_s['radius']

            if len(path) > 1:
                path = path[1:]
                path_radius = path_s['radius'][1:]

            path_type = path_s['type']
    #         
            connect_to = path_s['connect_to']
            connect_to_at = path_s['connect_to_at']

            swc_path = np.column_stack([np.ones(len(path)) * path_type, path]) # type
            swc_path = np.column_stack([np.arange(len(swc_arr)+1, len(path)+len(swc_arr)+1), swc_path]) #ID
            swc_path = np.column_stack([swc_path, path_radius * np.ones(len(path))]) # radius
            swc_path = np.column_stack([swc_path, swc_path[:, 0]-1]) # placeholder for PID

            if len(swc_arr) == 0:
                swc_arr = swc_path
                swc_arr[0][-1] = -1

            else:
                pid = np.where((swc_arr[:, 2:5] == connect_to_at).all(1))[0] + 1
#                 logging.info(pid)
                if len(pid)>1:
                    swc_path[0][-1] = pid[0]
#                     logging.info(swc_arr[pid[0]][2:5], swc_arr[pid[1]][2:5])
                else:
                    swc_path[0][-1] = pid
                swc_arr = np.vstack([swc_arr, swc_path])

            path_checked.append(path_id)
            
    df_swc = pd.DataFrame(swc_arr)
    df_swc.index = np.arange(1, len(df_swc)+1)
    df_swc.columns = ['n', 'type', 'x', 'y', 'z', 'radius', 'parent']
    df_swc[['n', 'type', 'parent']] = df_swc[['n', 'type', 'parent']].astype(int)
    
    return df_swc

def data_preprocessing(filepath):

    filetype = filepath.split('/')[-1].split('.')[-1].lower()
    filename = filepath.split('/')[-1].split('.')[0].lower()

    if filetype == 'swc':
        data = read_swc(filepath)
    elif filetype == 'imx':
        data = read_imx(filepath)
    else:
        logging.info('  `.{}` is not supported yet.'.format(filetype))
        raise NotImplementedError

    e = data['e']
    n = data['n']
    t = data['t']
    pos = data['pos']
    radius = data['radius']
    soma_loc = data['soma_loc']
    
    edge_dict = get_edge_dict(n, e, soma_loc)
    df_paths = get_path_dict(pos, radius, t, edge_dict, soma_loc)
    df_paths, df_paths_drop = sort_path_direction(df_paths)
    df_paths = reconnect_dropped_paths(df_paths, df_paths_drop)
    df_paths = find_connection(df_paths)
    df_swc = write_swc(df_paths)

    return df_swc, df_paths

# class Preprocessing(Morph):

#     def __init__(self, filepath):

#         filetype = filepath.split('/')[-1].split('.')[-1].lower()
#         filename = filepath.split('/')[-1].split('.')[0].lower()

#         if filetype == 'swc':
#             data = read_swc(filepath)
#         elif filetype == 'imx':
#             data = read_imx(filepath)
#         else:
#             logging.info('  `.{}` is not supported yet.'.format(filetype))
#             raise NotImplementedError

#         e = data['e']
#         n = data['n']
#         t = data['t']
#         pos = data['pos']
#         radius = data['radius']
#         soma_loc = data['soma_loc']
        
#         edge_dict = get_edge_dict(n, e, soma_loc)
#         df_paths = get_path_dict(pos, radius, t, edge_dict, soma_loc)
#         df_paths, df_paths_drop = sort_path_direction(df_paths)
#         df_paths = reconnect_dropped_paths(df_paths, df_paths_drop)
#         df_paths = find_connection(df_paths)
#         df_swc = write_swc(df_paths)
        
#         self.df_paths = df_paths
#         self.df_paths_drop = df_paths_drop
#         self.df_swc = df_swc
#         self.filename = filename

#     def save_as_swc(self, save_to='./'):

#         self.df_swc.to_csv(save_to + '{}.swc'.format(self.filename), sep=' ', index=None, header=None)