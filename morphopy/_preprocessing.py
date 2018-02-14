import numpy as np
import pandas as pd
import logging

logger = logging.getLogger()
logger.setLevel(logging.INFO)

__all__ = ['Preprocessing']

def read_swc(filepath):

    logging.info('  Reading {}\n'.format(filepath))
    
    swc =  pd.read_csv(filepath, delim_whitespace=True, comment='#',
                          names=['n', 'type', 'x', 'y', 'z', 'radius', 'parent'], index_col=False)
    swc.index = swc.n.as_matrix()

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

    logging.info("  Creating path arrays from voxel data points.")

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
                logging.info('\tPoint {} is already used by other paths. Skip.'.format(next_point))
            else:   
                a_list.append(next_point)# a list for holding the index of point of one paths
            
            next_point_locs_on_e = np.where(e[:, 0] == next_point)[0]
            
            while len(next_point_locs_on_e) == 1:
                
                next_point = e[next_point_locs_on_e[0]][1]
                if _point_already_in_dict(next_point, edge_dict):
                    logging.info('\tPoint {} is already used by other paths. Skip.'.format(next_point))
                else:   
                    a_list.append(next_point)# a list for holding the index of point of one paths
                next_point_locs_on_e = np.where(e[:, 0] == next_point)[0]
                
            if len(a_list) < 2:
                logging.info('\t\tCurrent path has fewer than two points, skipped.')
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

    logging.info('  Done.\n')
            
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

    logging.info('  Checking path connection.')
    logging.info('\tTotal num of paths to check: {}'.format(len(df_paths)))    
    logging.info('\tNumber of paths connected to soma: {}'.format(len(new_target_paths)))
    
    while np.count_nonzero(~np.isnan(df_paths.connect_to)) != len(df_paths):
        
        all_checked_paths = list(df_paths[~np.isnan(df_paths.connect_to)].index)
        num_check_paths_before = len(all_checked_paths)

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
        logging.info("\tNumber of paths checked: {}".format(num_check_paths_after))
        if num_check_paths_before == num_check_paths_after:
            num_disconneted = len(df_paths) - num_check_paths_after
            logging.info('\tNumber of disconnected path(s): {}'.format(num_disconneted))
            break
    
    df_paths_drop = df_paths[np.isnan(df_paths.connect_to)] 
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
    
    logging.info('  Done.\n')

    return df_paths, df_paths_drop

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
#                 print(pid)
                if len(pid)>1:
                    swc_path[0][-1] = pid[0]
#                     print(swc_arr[pid[0]][2:5], swc_arr[pid[1]][2:5])
                else:
                    swc_path[0][-1] = pid
                swc_arr = np.vstack([swc_arr, swc_path])

            path_checked.append(path_id)
            
    df_swc = pd.DataFrame(swc_arr)
    df_swc.index = np.arange(1, len(df_swc)+1)
    df_swc.columns = [['ID', 'Type', 'x', 'y', 'z', 'Raidus', 'PID']]
    df_swc[['ID', 'Type', 'PID']] = df_swc[['ID', 'Type', 'PID']].astype(int)
    
    return df_swc

class Preprocessing:

    def __init__(self, filepath):

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
        df_swc = write_swc(df_paths)
        
        self.df_paths = df_paths
        self.df_swc = df_swc

    def save_as_swc(self, filename, save_to='./'):

        self.df_swc.to_csv(save_to + '{}.swc'.format(filename), sep=' ', index=None, header=None)