import logging
import numpy as np

def get_path_dendritic_length(path):

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


def get_remote_vector(df_paths, path_id):

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

    s = df_paths.loc[path_id].path[0]
    e = df_paths.loc[path_id].path[-1]
    v= e-s
    return v/np.linalg.norm(v)


def get_local_vector(df_paths, path_id):
    
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

    s = df_paths.loc[path_id].path[0]
    e = df_paths.loc[path_id].path[1]
    v= e-s
    
    return v/np.linalg.norm(v)

def get_average_angles(df_paths):
    """

    :param df_paths:
    :return:
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
        
        # print(i, len(path_ids))
        if len(path_ids) == 2:
            v00 = get_remote_vector(df_paths, path_ids[0])
            v01 = get_remote_vector(df_paths, path_ids[1])
            nodal_angles_rad[n], nodal_angles_deg[n] = get_angle(v00, v01)

            v10 = get_local_vector(df_paths, path_ids[0])
            v11 = get_local_vector(df_paths, path_ids[1])
            local_angles_rad[n], local_angles_deg[n] = get_angle(v10, v11)

            n+=1
        else:
            continue

    average_nodal_angle_deg = np.mean(list(nodal_angles_deg.values()))
    average_nodal_angle_rad = np.mean(list(nodal_angles_rad.values()))

    average_local_angle_deg = np.mean(list(local_angles_deg.values()))
    average_local_angle_rad = np.mean(list(local_angles_rad.values()))

    return average_nodal_angle_deg, average_nodal_angle_rad, average_local_angle_deg, average_local_angle_rad

