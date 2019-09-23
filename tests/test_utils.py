import numpy as np
import sys

sys.path.append('..')

#### TEST GET_ANGLE #####
from morphopy._utils.summarize import get_angle


def test_get_angle_with_orthogonal_vectors():
    v0 = np.array([0, 0, 1])
    v1 = np.array([0, 1, 0])

    r, d = get_angle(v0, v1)
    assert(r == 90*np.pi/180), "returned angle should be pi/2"
    assert (d == 90), "returned angle should be 90 degree"


def test_get_angle_with_opposite_vectors():
    v0 = np.array([0, 0, 1])
    v1 = np.array([0, 0, -1])

    r, d = get_angle(v0, v1)
    assert (r == np.pi), "returned angle should be pi"
    assert (d == 180), "returned angle should be 180 degree"


def test_get_angle_with_same_vector():
    v0 = np.array([0, 0, 1])

    r, d = get_angle(v0, v0)
    assert (r == 0), "returned angle should be 0"
    assert (d == 0), "returned angle should be 0"


def test_get_angle_with_unnormalized_vector():
    v0 = np.array([0, 0, 1])
    v1 = np.array([0, 2, 0])

    r, d = get_angle(v0, v1)
    assert (r == 90 * np.pi / 180), "returned angle should be pi/2"
    assert (d == 90), "returned angle should be 90 degree"


def test_get_angle_btw_zero_and_v1():
    v0 = np.array([0, 0, 0])
    v1 = np.array([0, 1, 0])

    r, d = get_angle(v0, v1)
    assert (r == 0), "returned angle should be 0"
    assert (d == 0), "returned angle should be 0"


def test_get_angle_returns_float():
    v0 = np.array([0, 0, 1])
    v1 = np.array([0, 1, 1])

    r, d = get_angle(v0, v1)
    assert (isinstance(r, np.float)), "get_angle() should return float"
    assert (isinstance(d, np.float)), "get_angle() should return float"


# ### TEST READING METHODS ####

# import networkx as nx
# from morphopy._utils.utils import read_swc


# def test_read_swc_returned_fileformat():

#     import pandas as pd
#     filepath = 'data/Image001-005-01.CNG.swc'
#     G, swc = read_swc(filepath)

#     assert(isinstance(G, nx.DiGraph)), "read_swc() should return a graph as networkx.DiGraph"
#     assert(isinstance(swc, pd.DataFrame)), "read_swc() should return a swc as pandas.DataFrame"


# def test_read_swc_all_variables_are_in():

#     filepath = 'data/Image001-005-01.CNG.swc'
#     G, swc = read_swc(filepath)

#     assert 'n' in swc.keys(), "column 'n' should be in pandas.DataFrame"
#     assert 'x' in swc.keys(), "column 'x' should be in pandas.DataFrame"
#     assert 'y' in swc.keys(), "column 'y' should be in pandas.DataFrame"
#     assert 'z' in swc.keys(), "column 'z' should be in pandas.DataFrame"
#     assert 'type' in swc.keys(), "column 'type' should be in pandas.DataFrame"
#     assert 'parent' in swc.keys(), "column 'parent' should be in pandas.DataFrame"



# ### TEST FUNCTIONAL METHODS ###


# from morphopy._utils.utils import get_df_paths


# def test_get_df_paths_creates_dataFrame():
#     import pandas as pd
#     filepath = 'data/Image001-005-01.CNG.swc'
#     G, swc = read_swc(filepath)

#     paths = get_df_paths(G)
#     assert (isinstance(paths, pd.DataFrame)), "get_df_paths() should return a pandas.DataFrame"



from morphopy._utils.utils import unique_row


def test_unique_row_pairs_of_same_value():
    a = np.array([[9, 9], [8, 8], [1, 1], [9, 9]])

    a_ = unique_row(a)

    assert (a_ == np.array([[1, 1], [8, 8], [9, 9]])).all()


def test_unique_row_pairs_of_different_values():
    a = np.array([[1, 2], [2, 3], [2, 3], [9, 8]])

    a_ = unique_row(a)

    assert (a_ == np.array([[1, 2], [2, 3], [9, 8]])).all()


def test_unique_row_higher_number_first():
    a = np.array([[2, 1], [3, 6], [7, 4], [3, 6]])

    a_ = unique_row(a)

    assert (a_ == np.array([[2, 1], [3, 6], [7, 4]])).all()