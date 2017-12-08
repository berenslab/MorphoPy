import numpy as np
import sys

sys.path.append('..')

#### TEST GET_ANGLE #####
from morphopy._utils import get_angle


def test_get_angle_with_orthogonal_vectors():
    v0 = np.array([0,0,1])
    v1 = np.array([0, 1, 0])

    r, d = get_angle(v0,v1)
    assert(r == 90*np.pi/180)
    assert (d == 90)


def test_get_angle_with_opposite_vectors():
    v0 = np.array([0, 0, 1])
    v1 = np.array([0, 0, -1])

    r, d = get_angle(v0, v1)
    assert (r == np.pi)
    assert (d == 180)


def test_get_angle_with_same_vector():
    v0 = np.array([0, 0, 1])

    r, d = get_angle(v0, v0)
    assert (r == 0)
    assert (d == 0)


def test_get_angle_with_unnormalized_vector():
    v0 = np.array([0, 0, 1])
    v1 = np.array([0, 2, 0])

    r, d = get_angle(v0, v1)
    assert (r == 90 * np.pi / 180)
    assert (d == 90)


def test_angle_btw_zero_and_v1():
    v0 = np.array([0, 0, 0])
    v1 = np.array([0, 1, 0])

    r, d = get_angle(v0, v1)
    assert (r == 0)
    assert (d == 0)


### TEST READING METHODS ####

from morphopy._utils import read_swc


def test_read_swc_returned_fileformat():

    import pandas as pd
    filepath = './data/Image001-005-01.CNG.swc'
    df = read_swc(filepath)

    assert(isinstance(df, pd.DataFrame))


def test_read_swc_all_variables_are_in():

    filepath = './data/Image001-005-01.CNG.swc'
    swc = read_swc(filepath)

    assert 'n' in swc.keys()
    assert 'x' in swc.keys()
    assert 'y' in swc.keys()
    assert 'z' in swc.keys()
    assert 'type' in swc.keys()
    assert 'parent' in swc.keys()

### TEST FUNCTIONAL METHODS ###

from morphopy._utils import get_consecutive_pairs_of_elements_from_list


def test_get_consecutive_pairs_of_elements_from_list_default():


    l = [1, 2, 3, 4]

    result = get_consecutive_pairs_of_elements_from_list(l, s=None, e=None)

    assert result == [(1, 2), (2, 3), (3, 4)]


def test_get_consecutive_pairs_of_elements_from_list_with_start_node():

    l = [1, 2, 3, 4]

    result = get_consecutive_pairs_of_elements_from_list(l, s=0, e=None)

    assert result == [(0, 1), (1, 2), (2, 3), (3, 4)]


def test_get_consecutive_pairs_of_elements_from_list_with_end_node():

    l = [1, 2, 3, 4]

    result = get_consecutive_pairs_of_elements_from_list(l, s=None, e=5)

    assert result == [(1, 2), (2, 3), (3, 4), (4, 5)]


def test_get_consecutive_pairs_of_elements_from_list_with_start_and_end_node():

    l = [1, 2, 3, 4]

    result = get_consecutive_pairs_of_elements_from_list(l, s=0, e=5)

    assert result == [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)]

