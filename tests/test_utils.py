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