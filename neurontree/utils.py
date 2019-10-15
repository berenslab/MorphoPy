import numpy as np
import scipy
import math

from scipy.sparse import coo_matrix
from sys import getsizeof

from scipy.signal import gaussian, convolve2d
from scipy.ndimage.filters import gaussian_filter


def unit_vector(vector):
    """ Returns the unit vector of the vector.

    Parameter
    --------
    vector : numpy.array
        d-dimensional vector

    Returns
    -------
    u : numpy.array
        d-dimensional unit vector of 'vector'

    """
    if np.linalg.norm(vector) > 0:
        return vector / np.linalg.norm(vector)
    else:
        return vector


def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    Parameters
    ----------
    v1 : numpy.array
        d-dimensional vector
    v2 : numpy.array
        d-dimensional vector

    Returns
    -------
    theta : float
        angle btw v1 and v2 in radians

    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arctan2(np.linalg.norm(np.cross(v1_u, v2_u)), np.dot(v1_u, v2_u))


def get_axis(v1, v2):
    """ Returns the axis between two vectors v1 and v2

    Parameters
    ---------
    v1: numpy.array
        d-dimensional
    v2: numpy.array
        d-dimensional

    Returns
    -------
    axis: numpy.array
        d-dimensional axis btw v1 and v2
    """
    return (np.cross(v1, v2)) / (np.linalg.norm(np.cross(v1, v2)))


def get_A(x):
    """
    Returns matrix A. Used for rotation matrix calculation
    :param x: 3D numpy.array
            rotation axis

    :return:
        A: 3x3 numpy.array
    """

    if np.isnan(x).any():
        A = np.zeros((3, 3))
    else:
        A = np.array([[0, -x[2], x[1]],
                  [x[2], 0, -x[0]],
                  [-x[1], x[0], 0]])

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
        theta = angle_between(a_,b_)
        if n == 3:
            A = get_A(x)
            R += np.sin(theta) * A + (1 - np.cos(theta)) * np.linalg.matrix_power(A, 2)
        if n == 2:
            R = [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
    return np.array(R)


def isRotationMatrix(R):
    """
    Checks if R is a valid rotation matrix.
    :param R: [3x3] rotation matrix
    :return: boolean. Returns True when R is a valid rotation matrix, otherwise False.
    """
    Rt = np.transpose(R)
    shouldBeIdentity = np.dot(Rt, R)
    I = np.identity(3, dtype=R.dtype)
    n = np.linalg.norm(I - shouldBeIdentity)
    return n < 1e-6


def rotationMatrixToEulerAngles(R):
    """
    Calculates rotation matrix to euler angles
    :param R: rotation matrix [3x3]
    :return: vector of euler angles (rotations around x,y and z axis)
    """
    assert (isRotationMatrix(R))

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
    """
    Calculates the rotation matrix from given euelr angles
    :param theta: 3D vector with eulerangles for x,y and z axis
    :return: rotation matrix R [3x3]
    """
    R_x = np.array([[1, 0, 0],
                    [0, math.cos(theta[0]), -math.sin(theta[0])],
                    [0, math.sin(theta[0]), math.cos(theta[0])]
                    ])

    R_y = np.array([[math.cos(theta[1]), 0, math.sin(theta[1])],
                    [0, 1, 0],
                    [-math.sin(theta[1]), 0, math.cos(theta[1])]
                    ])

    R_z = np.array([[math.cos(theta[2]), -math.sin(theta[2]), 0],
                    [math.sin(theta[2]), math.cos(theta[2]), 0],
                    [0, 0, 1]
                    ])

    R = np.dot(R_z, np.dot(R_y, R_x))

    return R


def sphereFit(spX, spY, spZ):
    """
    fit a sphere to X,Y, and Z data points returns the radius and center points of the best fit sphere. Implementation
    by http://jekel.me/2015/Least-Squares-Sphere-Fit/
    :param spX:
    :param spY:
    :param spZ:
    :return:
    """
    #   Assemble the A matrix
    spX = np.array(spX)
    spY = np.array(spY)
    spZ = np.array(spZ)
    A = np.zeros((len(spX),4))
    A[:,0] = spX*2
    A[:,1] = spY*2
    A[:,2] = spZ*2
    A[:,3] = 1

    #   Assemble the f matrix
    f = np.zeros((len(spX),1))
    f[:,0] = (spX*spX) + (spY*spY) + (spZ*spZ)
    C, residules, rank, singval = np.linalg.lstsq(A,f)

    #   solve for the radius
    t = (C[0]*C[0])+(C[1]*C[1])+(C[2]*C[2])+C[3]
    radius = math.sqrt(t)

    return radius, C[0], C[1], C[2]


def shortpathFW(A, gpu=True):
    """
    Returns matrix B of shortest path distances from each node to each other node in the graph.
    If the graph is undirected, B will be symmetric.
    B[:,0] corresponds to all shortest paths from reached from node 0.
    Algorithm after Floyd and Warshall

    :param A: adjacency matrix of a graph (NxN)
    :param gpu: default True. If True, calculation is performed on GPU.
    :return: B (NxN) matrix of shortest paths
    """
    print('Calculating shortest path from every node to every other...')
    n = A.shape[0]
    sparse_ones = coo_matrix(np.ones((n, n)), dtype=np.float32)
    I = np.inf * (sparse_ones - np.matlib.eye(n))
    I = np.nan_to_num(I)
    I[A != 0] = np.squeeze(np.array(A[A != 0]))
    B = np.array(I, dtype=np.float32)
    print(getsizeof(A) * 1e-09)

    if getsizeof(A) * 1e-09 * 2 > 11.5:
        print('input is too large to be computed on GPU. Switching to CPU instead.')
        gpu = False

        #if gpu:
        # g = tf.Graph()
        # with g.as_default():
        #
        #     X = tf.placeholder(dtype=tf.float32, shape=(n, n), name='X')
        #     B_ = tf.Variable(X, name='B_')
        #     k = tf.Variable(0, dtype=tf.int32, name='k')
        #
        #     C1 = tf.transpose(tf.reshape(tf.tile(B_[:, k], [n]), (n, -1)), name='C1')
        #     C2 = tf.reshape(tf.tile(B_[k, :], [n]), (-1, n), name='C2')
        #     C = tf.add(C1, C2, name='C')
        #     D = tf.minimum(B_, C)
        #     step = tf.group(B_.assign(D), tf.assign_add(k, 1))
        #
        # # ensure that only as much gpu memory is needed as necessary
        # with tf.Session(graph=g) as sess:
        #     # Initialize state to initial conditions
        #     tf.global_variables_initializer().run(feed_dict={X: B})
        #
        #     for i in range(n - 1):  # one step less to not exceed array range
        #         sess.run(step, feed_dict={X: B})
        #         B = B_.eval()
        #else:
    for k in range(n):
        C1 = np.tile(B[:, k], (n, 1))
        C2 = np.tile(B[k, :], (n, 1))
        C = C1 + C2
        B = np.matlib.minimum(B, C)

    print('done.')
    return B


def commuteDist(A):
    """
    Returns the commute distance within a graph given by adjacency matrix A
    :param A: adjacency matrix (NxN)
    :return: B (NxN) containing the commute distances from each node to each other
    """
    n = A.shape[0]
    if n==1:
        return 0
    m = np.mean(A)
    A = scipy.linalg.expm(A/m)
    L = np.diag(np.sum(A, axis=0)) - A # unnormalised Laplacian
    L = np.linalg.pinv(L)
    d = L.diagonal()
    print(d.shape)
    B = np.matlib.repmat(d,1,n) + np.matlib.repmat(d.T, n,1) - 2*L
    return B


def computeStat(statType, W, d, maxDist):
    """
    computes a graph statistic on (sub)graph given by W
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
        statistic calculated according to statType
    """

    numNode = d.size

    if statType == 'maxDist':
        stat = np.max(d)

    elif statType == 'maxDist_norm':
        stat = np.max(d) / maxDist

    elif statType == 'meanEdgeLength':
        stat = np.mean(W)

    elif statType == 'meanEdgeLength_norm':
        stat = np.mean(W) / numNode

    elif statType == '4starMotif':
        if W.size > 1:
            deg = np.array(np.sum(W != 0, axis=1)).flatten()
        else:
            deg = np.sum(W != 0)
        stat = np.sum((deg - 1) * (deg - 2) / 2)

    elif statType == 'branchPoints':

        stat = np.sum(np.sum(W != 0, axis=1) > 2)

    else:
        raise ValueError("Calculation of statistic {0} is not implemented.".format(statType))

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
    N=1
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
            Xb[:, :, k] = convolve2d(data[:, :, k], np.rot90(win), mode='same')

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
        raise NotImplementedError("There is no gaussian smoothing implemented for {0} dimensions".format(dim))

    if N == 1:
        Xb = np.squeeze(Xb)

    return Xb
