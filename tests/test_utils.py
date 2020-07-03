import numpy as np

#### TEST GET_ANGLE #####
from morphopy.neurontree.utils import angle_between


def test_get_angle_with_orthogonal_vectors():
    v0 = np.array([0, 0, 1])
    v1 = np.array([0, 1, 0])

    r = angle_between(v0, v1)
    assert(r == 90*np.pi/180), "returned angle should be pi/2"


def test_get_angle_with_opposite_vectors():
    v0 = np.array([0, 0, 1])
    v1 = np.array([0, 0, -1])

    r= angle_between(v0, v1)
    assert (r == np.pi), "returned angle should be pi"


def test_get_angle_with_same_vector():
    v0 = np.array([0, 0, 1])

    r = angle_between(v0, v0)
    assert (r == 0), "returned angle should be 0"


def test_get_angle_with_unnormalized_vector():
    v0 = np.array([0, 0, 1])
    v1 = np.array([0, 2, 0])

    r = angle_between(v0, v1)
    assert (r == 90 * np.pi / 180), "returned angle should be pi/2"


def test_get_angle_btw_zero_and_v1():
    v0 = np.array([0, 0, 0])
    v1 = np.array([0, 1, 0])

    r = angle_between(v0, v1)
    assert (r == 0), "returned angle should be 0"
    
    
def test_get_angle_returns_float():
    v0 = np.array([0, 0, 1])
    v1 = np.array([0, 1, 1])

    r  = angle_between(v0, v1)
    assert (isinstance(r, np.float)), "get_angle() should return float"




  

