import math
from sys import getsizeof

import numpy as np
import pandas as pd
import scipy
from scipy.ndimage import gaussian_filter
from scipy.signal import convolve2d, gaussian
from scipy.sparse import coo_matrix
from scipy.spatial import ConvexHull
from sklearn.decomposition import PCA


def unit_vector(vector):
    """Returns the unit vector of the vector.

    :param vector: numpy.array d-dimensional vector
    :return: u numpy.array d-dimensional unit vector of 'vector'
    """
    if np.linalg.norm(vector) > 0:
        return vector / np.linalg.norm(vector)
    else:
        return vector


def angle_between(v1, v2):
    """Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793

    :param v1: numpy.array d-dimensional vector
    :param v2: numpy.array d-dimensional vector
    :return: theta - float angle btw v1 and v2 in radians

    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arctan2(np.linalg.norm(np.cross(v1_u, v2_u)), np.dot(v1_u, v2_u))


def get_axis(v1, v2):
    """Returns the axis between two vectors v1 and v2

    :param v1: numpy.array d-dimensional
    :param v2: numpy.array d-dimensional
    :return: axis - numpy.array d-dimensional axis btw v1 and v2
    """
    return (np.cross(v1, v2)) / (np.linalg.norm(np.cross(v1, v2)))


def get_A(x):
    """
    Returns matrix A. Used for rotation matrix calculation

    :param x: 3D numpy.array rotation axis
    :return: A - 3x3 numpy.array
    """

    if np.isnan(x).any():
        A = np.zeros((3, 3))
    else:
        A = np.array([[0, -x[2], x[1]], [x[2], 0, -x[0]], [-x[1], x[0], 0]])

    return A


def get_rotation_matrix(a, b):
    """
    Returns the rotation matrix to rotate vector a onto vector b.

    :param a: numpy.array (2 or 3 dimensional)
    :param b: numpy.array (2 or 3 dimensional)
    :return: R (2x2 or 3x3) rotation matrix to rotate a onto b
    """

    n = np.max(a.shape)
    R = np.eye(n)
    a_ = unit_vector(a)
    b_ = unit_vector(b)
    if not np.allclose(np.abs(a_), np.abs(b_)):
        x = get_axis(a_, b_)
        theta = angle_between(a_, b_)
        if n == 3:
            A = get_A(x)
            R += np.sin(theta) * A + (1 - np.cos(theta)) * np.linalg.matrix_power(A, 2)
        if n == 2:
            R = [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
    return np.array(R)


def isRotationMatrix(R):
    """Checks if R is a valid rotation matrix.

    :param R: [3x3] rotation matrix
    :return: boolean. Returns True when R is a valid rotation matrix, otherwise False.
    """
    Rt = np.transpose(R)
    shouldBeIdentity = np.dot(Rt, R)
    I = np.identity(3, dtype=R.dtype)
    n = np.linalg.norm(I - shouldBeIdentity)
    return n < 1e-6


def rotationMatrixToEulerAngles(R):
    """Calculates rotation matrix to euler angles

    :param R: rotation matrix [3x3]
    :return: vector of euler angles (rotations around x,y and z axis)
    """
    assert isRotationMatrix(R)

    sy = math.sqrt(R[0, 0] * R[0, 0] + R[1, 0] * R[1, 0])

    singular = sy < 1e-6

    if not singular:
        x = math.atan2(R[2, 1], R[2, 2])
        y = math.atan2(-R[2, 0], sy)
        z = math.atan2(R[1, 0], R[0, 0])
    else:
        x = math.atan2(-R[1, 2], R[1, 1])
        y = math.atan2(-R[2, 0], sy)
        z = 0

    return np.array([x, y, z])


def eulerAnglesToRotationMatrix(theta):
    """Calculates the rotation matrix from given euler angles

    :param theta: 3D vector with eulerangles for x,y and z axis
    :return: rotation matrix R [3x3]
    """
    R_x = np.array(
        [
            [1, 0, 0],
            [0, math.cos(theta[0]), -math.sin(theta[0])],
            [0, math.sin(theta[0]), math.cos(theta[0])],
        ]
    )

    R_y = np.array(
        [
            [math.cos(theta[1]), 0, math.sin(theta[1])],
            [0, 1, 0],
            [-math.sin(theta[1]), 0, math.cos(theta[1])],
        ]
    )

    R_z = np.array(
        [
            [math.cos(theta[2]), -math.sin(theta[2]), 0],
            [math.sin(theta[2]), math.cos(theta[2]), 0],
            [0, 0, 1],
        ]
    )

    R = np.dot(R_z, np.dot(R_y, R_x))

    return R


def sphereFit(spX, spY, spZ):
    """
    fit a sphere to X,Y, and Z data points returns the radius and center points of the best fit sphere.\n
    Implementation by http://jekel.me/2015/Least-Squares-Sphere-Fit/

    :param spX:
    :param spY:
    :param spZ:
    :return:
    """
    #   Assemble the A matrix
    spX = np.array(spX)
    spY = np.array(spY)
    spZ = np.array(spZ)
    A = np.zeros((len(spX), 4))
    A[:, 0] = spX * 2
    A[:, 1] = spY * 2
    A[:, 2] = spZ * 2
    A[:, 3] = 1

    #   Assemble the f matrix
    f = np.zeros((len(spX), 1))
    f[:, 0] = (spX * spX) + (spY * spY) + (spZ * spZ)
    C, residules, rank, singval = np.linalg.lstsq(A, f)

    #   solve for the radius
    t = (C[0] * C[0]) + (C[1] * C[1]) + (C[2] * C[2]) + C[3]
    radius = math.sqrt(t)

    return radius, C[0], C[1], C[2]


def shortpathFW(A):
    """
    Returns matrix B of shortest path distances from each node to each other node in the graph.
    If the graph is undirected, B will be symmetric.
    B[:,0] corresponds to all shortest paths from reached from node 0.
    Algorithm after Floyd and Warshall

    :param A: adjacency matrix of a graph (NxN)
    :param gpu: default True. If True, calculation is performed on GPU.
    :return: B (NxN) matrix of shortest paths
    """
    print("Calculating shortest path from every node to every other...")
    n = A.shape[0]
    sparse_ones = coo_matrix(np.ones((n, n)), dtype=np.float32)
    I = np.inf * (sparse_ones - np.matlib.eye(n))
    I = np.nan_to_num(I)
    I[A != 0] = np.squeeze(np.array(A[A != 0]))
    B = np.array(I, dtype=np.float32)
    print(getsizeof(A) * 1e-09)

    for k in range(n):
        C1 = np.tile(B[:, k], (n, 1))
        C2 = np.tile(B[k, :], (n, 1))
        C = C1 + C2
        B = np.matlib.minimum(B, C)

    print("done.")
    return B


def commuteDist(A):
    """
    Returns the commute distance within a graph given by adjacency matrix A

    :param A: adjacency matrix (NxN)
    :return: B (NxN) containing the commute distances from each node to each other
    """
    n = A.shape[0]
    if n == 1:
        return 0
    m = np.mean(A)
    A = scipy.linalg.expm(A / m)
    L = np.diag(np.sum(A, axis=0)) - A  # unnormalised Laplacian
    L = np.linalg.pinv(L)
    d = L.diagonal()
    print(d.shape)
    B = np.matlib.repmat(d, 1, n) + np.matlib.repmat(d.T, n, 1) - 2 * L
    return B


def computeStat(statType, W, d, maxDist):
    """computes a graph key on (sub)graph given by W

    :param statType: string
        'maxDist' : maximal distance
        'maxDist_norm' : normalized maximal distance
        'meanEdgeLength' : mean edge length
        'meanEdgeLength_norm' : normalized mean edge length
        '4starMotif' : 4 star motif
        'branchPoints' : number of branch points
    :param W: (KxK)
        adjacency matrix of graph
    :param d: (Kx1)
        shortest path distances of nodes in W
    :param maxDist: float
        maximal shortest path distance in graph of which W is a sub-graph of.
    :return: stat: float
        key calculated according to statType
    """

    numNode = d.size

    if statType == "maxDist":
        stat = np.max(d)

    elif statType == "maxDist_norm":
        stat = np.max(d) / maxDist

    elif statType == "meanEdgeLength":
        stat = np.mean(W)

    elif statType == "meanEdgeLength_norm":
        stat = np.mean(W) / numNode

    elif statType == "4starMotif":
        if W.size > 1:
            deg = np.array(np.sum(W != 0, axis=1)).flatten()
        else:
            deg = np.sum(W != 0)
        stat = np.sum((deg - 1) * (deg - 2) / 2)

    elif statType == "branchPoints":

        stat = np.sum(np.sum(W != 0, axis=1) > 2)

    else:
        raise ValueError("Calculation of key {0} is not implemented.".format(statType))

    return stat


def smooth_gaussian(data, dim, sigma=2):
    """
    Smooths the given data using a gaussian. This method only works for stacked one or two dimensional data so
    far! Smoothing in 3D is not implemented.

    :param data: (X,Y,N) numpy.array 1,2 or 3 dimensional.
    :param dim: int
        Dimension of the passed data. Used to determine if data is a single image or stacked.
    :param sigma: int
        The standard deviation of the smoothing gaussian used.
    :return: Xb: same dimension as the input array.
        Smoothed data.
    """
    N = 1
    if dim == 2:
        try:
            pX, pY, N = data.shape
        except ValueError:
            pX, pY = data.shape
            data = data.reshape(pX, pY, 1)

        # gaussian window
        win = gaussian(11, sigma)
        win = win * win.reshape(-1, 1)

        # add blur
        Xb = np.zeros((pX, pY, N))

        for k in range(N):
            Xb[:, :, k] = convolve2d(data[:, :, k], np.rot90(win), mode="same")

    elif dim == 1:

        try:
            pX, N = data.shape
        except ValueError:
            pX = data.shape[0]
            data = data.reshape(pX, 1)

        # add blur
        Xb = np.zeros((pX, N))

        for k in range(N):
            Xb[:, k] = gaussian_filter(data[:, k], sigma=sigma)
    else:
        raise NotImplementedError(
            "There is no gaussian smoothing implemented for {0} dimensions".format(dim)
        )

    if N == 1:
        Xb = np.squeeze(Xb)

    return Xb


def get_standardized_swc(
    swc, scaling=1.0, soma_radius=None, soma_center=True, pca_rot=False
):
    """
    This function collapses all soma points to a single node located at the centroid of the convex hull of the original
    soma nodes. It can also scale the coordinates, merge nodes into the soma that have a bigger radius than soma_radius
    and return the xyz coordinates of the swc file in their PCA rotation.

    :param swc: swc file, as pandas.DataFrame.
    :param scaling: float, (default=1), allows for a uniform scaling of x, y and z
    :param soma_radius: float (default=None), if set, then all nodes with a radius greater or equal to soma radius are
     set to be somatic (type=1). In a subsequent step these nodes will be merged to one. Careful! If this radius is set
     too small this leads to faulty skeletons.
    :param soma_center: bool (default=True), if True, x,y,z are soma centered.
    :param pca_rot: bool (default=False), if True, the x,y,z coordinates in the given swc file are rotated into their
     PCA frame. Then x corresponds to the direction of highest and z to the direction of lowest variance.
    :return: pandas.DataFrame
    """

    swc.update(swc["x"] / scaling)
    swc.update(swc["y"] / scaling)
    swc.update(swc["z"] / scaling)
    swc.update(swc["radius"] / scaling)

    if pca_rot:
        print("Rotating x and y into their frame of maximal extent...")
        pca = PCA(copy=True)
        pc = np.vstack((swc["x"], swc["y"])).T
        pca.fit(pc)
        result = np.matmul(pc, pca.components_.T)

        swc.update(pd.DataFrame(result, columns=["x", "y"]))

    if soma_radius:
        print(
            "Setting all nodes to type soma that have a larger radius than %s microns..."
            % soma_radius
        )

        d = np.vstack((swc["radius"], swc["type"])).T
        d[d[:, 0] >= soma_radius, 1] = 1
        swc.update(pd.DataFrame(d[:, 1].astype(int), columns=["type"]))

    # create one point soma when there is more than three soma points
    sp = swc[swc["type"] == 1]

    if sp.shape[0] > 1:
        root_id = np.min(sp["n"].values)

        if sp.shape[0] > 3:
            print(
                "There are more than 3 soma points. The location and the radius of the soma is estimated based on its"
                " convex hull..."
            )

            # calculate the convex hull of soma points
            convex_hull = ConvexHull(sp[["x", "y", "z"]].values, qhull_options="QJ")

            hull_points = convex_hull.points

            centroid = np.mean(hull_points, axis=0)

            distances_to_centroid = np.linalg.norm(hull_points - centroid, axis=1)
            rad = np.max(distances_to_centroid)
        else:
            print(
                "There are 2-3 soma points. The location and the radius of the soma is estimated based on their mean."
            )

            centroid = np.mean(sp[["x", "y", "z"]].values, axis=0)
            rad = np.mean(sp[["radius"]].values)

        # fix the parent connections
        connection_locations = [
            row.n
            for k, row in swc.iterrows()
            if row["parent"] in sp["n"].values and row["n"] not in sp["n"].values
        ]
        connected_points = pd.concat([swc[swc["n"] == n] for n in connection_locations])
        connected_points["parent"] = root_id
        swc.update(connected_points)

        # delete old soma points
        to_delete = [swc[swc["n"] == n].index[0] for n in sp["n"].values]
        swc = swc.drop(swc.index[to_delete])

        # add new soma point
        soma_dict = dict(
            zip(
                ["n", "type", "x", "y", "z", "radius", "parent"],
                [
                    int(root_id),
                    int(1),
                    centroid[0],
                    centroid[1],
                    centroid[2],
                    rad,
                    int(-1),
                ],
            )
        )

        swc = pd.concat([swc, pd.DataFrame(soma_dict, index=[0])])
        swc = swc.sort_index()
    else:
        # if no soma was assigned use the node that has no parent with the smallest id
        root_id = np.min(swc[swc["parent"] == -1]["n"].values)

    # soma center on first entry
    if soma_center:
        centroid = swc[swc["n"] == root_id][["x", "y", "z"]].values.reshape(-1)
        swc.update(swc[["x", "y", "z"]] - centroid)
    return swc
