import logging
import numpy as np

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

def get_path_real_length(path):

    """
    Get the dendritic length of a path, which is the sum of the distance between each consecutive points.

    Parameters
    ----------
    path: array 
        a coordinate array with dim=(n, 3)
    
    Returns
    -------
    the dendritic length of this path: float

    """

    return np.sum(np.sqrt(np.sum((path[1:] - path[:-1])**2, 1)))


def get_path_euclidean_length(path):
    """
    get the euclidean length of a path, which is the distance between the first and last points.

    Parameters
    ----------
    path: array 
        a coordinate array with dim=(n, 3)
    
    Returns
    -------
    the euclidean length of this path: float

    """
    return np.sqrt(np.sum((path[0] - path[-1]) ** 2))


def get_outer_terminals(all_terminals):

    """
    Get terminal points which form the convex hull of the cell.

    Parameters
    ----------
    all_terminals: array
        The array contains all terminal points from terminal paths (no other paths connected to them) 
    
    Returns
    -------
    outer_terminals_3d: array
        The array contains all terminal points which found the convex hull of the cell.

    """
    
    from scipy.spatial import ConvexHull
    hull = ConvexHull(all_terminals[:,:2])
    outer_terminals_3d = all_terminals[hull.vertices]
    outer_terminals_3d = np.vstack([outer_terminals_3d, outer_terminals_3d[0]])
    
    return outer_terminals_3d

def get_angle(v0, v1):

    """
    Get angle (in both radian and degree) between two vectors.
    
    Parameters
    ----------
    v0: array
        vector zero.
    v1: array
        vector one.

    Returns
    -------
    Return a tuple, (angle in radian, angle in degree). 

    """
    v0 = np.array(v0)
    v1 = np.array(v1)

    if not v0.any() or not v1.any():
        return 0, 0

    c = np.dot(v0, v1) / np.linalg.norm(v0) / np.linalg.norm(v1)
    return np.arccos(np.clip(c, -1, 1)), np.degrees(np.arccos(np.clip(c, -1, 1)))


def get_remote_vector(path):

    """
    Get vector of certain path between the first and the last point.
    
    Parameters
    ----------
    df_paths: pandas.DataFrame

    path_id: int

    Returns
    -------
    normalized v: array
        returned a normalized vector.
    """

    s = path[0]
    e = path[-1]
    v= e-s

    return v/np.linalg.norm(v)

def get_local_vector(path):
    
    """
    Get vector of certain path between the first and the second point.
    
    Parameters
    ----------
    df_paths: pandas.DataFrame

    path_id: int

    Returns
    -------
    normalized v: array
        returned a normalized vector.
    """

    s = path[0]
    e = path[1]
    v= e-s
    
    return v/np.linalg.norm(v)

def get_average_angles(df_paths):
    """
    a helper function to get the average of all kinds of angles.

    Parameters
    ----------
    df_paths: pandas.DataFrame

    Returns
    -------
    average_nodal_angle_deg 
    average_nodal_angle_rad
    average_local_angle_deg
    average_local_angle_rad    

    """

    nodal_angles_deg = {}
    nodal_angles_rad = {}
    
    local_angles_deg = {}
    local_angles_rad = {}
    
    n = 0
    for i in np.unique(df_paths.connect_to):

        if i == -1:
            continue

        path_ids = df_paths[df_paths.connect_to == i].index.tolist()
        
        if len(path_ids) == 2:
            
            p0 = df_paths.loc[path_ids[0]].path
            p1 = df_paths.loc[path_ids[1]].path
            
            v00 = get_remote_vector(p0)
            v01 = get_remote_vector(p1)
            nodal_angles_rad[n], nodal_angles_deg[n] = get_angle(v00, v01)

            v10 = get_local_vector(p0)
            v11 = get_local_vector(p1)
            local_angles_rad[n], local_angles_deg[n] = get_angle(v10, v11)

            n+=1
        else:
            continue

    average_nodal_angle_deg = np.mean(list(nodal_angles_deg.values()))
    average_nodal_angle_rad = np.mean(list(nodal_angles_rad.values()))

    average_local_angle_deg = np.mean(list(local_angles_deg.values()))
    average_local_angle_rad = np.mean(list(local_angles_rad.values()))

    return average_nodal_angle_deg, average_nodal_angle_rad, average_local_angle_deg, average_local_angle_rad

def get_summary_of_paths(df_paths):

    """
    A helper function to gether all summarized infomation

    Parameters
    ----------
    df_paths: pandas.DataFrame

    Return
    ------
    a list of all summarized information.

    """
    
    if len(df_paths) < 1:
        return None
    
    branchpoints = np.vstack(df_paths.connect_to_at)
    branchpoints = unique_row(branchpoints)
    num_branchpoints = len(branchpoints)

    max_branch_order = max(df_paths.branch_order)

    terminalpaths = df_paths.path[df_paths.connected_by.apply(len) == 0].as_matrix()
    terminalpoints = np.vstack([p[-1] for p in terminalpaths])
    num_terminalpoints = len(terminalpoints)

    outerterminals = get_outer_terminals(terminalpoints)

    num_irreducible_nodes = num_branchpoints + num_terminalpoints

    num_dendritic_segments = len(df_paths)

    # path length

    reallength = df_paths['real_length']
    reallength_sum = reallength.sum()
    reallength_mean = reallength.mean()
    reallength_median = reallength.median()
    reallength_min = reallength.min()
    reallength_max = reallength.max()

    euclidean = df_paths['euclidean_length']
    euclidean_sum = euclidean.sum()
    euclidean_mean = euclidean.mean()
    euclidean_median = euclidean.median()
    euclidean_min = euclidean.min()
    euclidean_max = euclidean.max()

    tortuosity = reallength / euclidean
    average_tortuosity = np.mean(tortuosity)      

    # node angles
    average_nodal_angle_deg, average_nodal_angle_rad, average_local_angle_deg, average_local_angle_rad = get_average_angles(df_paths)    

    if df_paths.iloc[0].type == 2:
        t = 'Axon'
    elif df_paths.iloc[0].type == 3:
        t = 'Basal Dendrites'
    elif df_paths.iloc[0].type == 4:
        t = 'Apical Dendrites'
    else:
        t = 'Undefined'
    
    return (t,int(num_dendritic_segments),
            int(num_branchpoints),
            int(num_irreducible_nodes),
            int(max_branch_order),
            average_nodal_angle_deg,
            average_nodal_angle_rad,
            average_local_angle_deg,
            average_local_angle_rad,
            average_tortuosity,
            reallength_sum,
            reallength_mean,
            reallength_median,
            reallength_min,
            reallength_max,
            euclidean_sum,
            euclidean_mean,
            euclidean_median,
            euclidean_min,
            euclidean_max,)