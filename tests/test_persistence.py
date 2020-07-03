from morphopy.computation.file_manager import load_swc_file
from morphopy.computation.feature_presentation import get_persistence
from morphopy.computation.persistence_functions import radial_distance, height, path_length, branch_order
import os
import pandas as pd

dir_path = os.path.dirname(os.path.realpath(__file__))

file_names = ['Image001-005_01_CNenhance_latest_LXS.swc', 'C4.swc',
              'V1_Layer23_Chat-IRES-Cre-neo_Ai14-299537.04.02.01_614430666_m.swc', 'EC3-80604.CNG.rotated.swc']

precomputed_data_radial = pd.read_csv(dir_path + "/test_data/precomputed_data_persistence_radial.csv")
precomputed_data_height = pd.read_csv(dir_path + "/test_data/precomputed_data_persistence_height.csv")
precomputed_data_path = pd.read_csv(dir_path + "/test_data/precomputed_data_persistence_path.csv")
precomputed_data_branch_order = pd.read_csv(dir_path + "/test_data/precomputed_data_persistence_branch_order.csv")


def compare_old_and_current_persistence_implementation(filename, old_data, distance_function):
    """
    HELPER function
    :param filename: sting name of the file that is loaded
    :param old_data: DataFrame holding the original outcome of the persistence  computation
    :param distance_function: distance function used for persistence data.
    """

    neuron = load_swc_file("../data/%s" % filename)
    columns = ['birth', 'death', 'node_type']
    old_persistence = old_data[old_data['filename'] == filename].sort_values(columns)[columns].reset_index()
    del old_persistence['index']

    current_persistence = get_persistence(neuron, f=distance_function).sort_values(columns)[columns].reset_index()
    del current_persistence['index']

    assert old_persistence.round(4).equals(current_persistence.round(4)), \
        "Computation of persistence with distance function %s has changed for file %s." % (filename, distance_function.__name__)

### TEST RADIAL PERSISTENCE ####


def test_radial_persistence_data_file_1():

    k = 0
    filename = file_names[k]
    compare_old_and_current_persistence_implementation(filename, precomputed_data_radial, radial_distance)


def test_radial_persistence_data_file_2():

    k = 1
    filename = file_names[k]
    compare_old_and_current_persistence_implementation(filename, precomputed_data_radial, radial_distance)


def test_radial_persistence_data_file_3():

    k = 2
    filename = file_names[k]
    compare_old_and_current_persistence_implementation(filename, precomputed_data_radial, radial_distance)

    
def test_radial_persistence_data_file_4():

    k = 3
    filename = file_names[k]
    compare_old_and_current_persistence_implementation(filename, precomputed_data_radial, radial_distance)

### TEST HEIGHT PERSISTENCE ###


def test_height_persistence_data_file_1():

    k = 0
    filename = file_names[k]
    compare_old_and_current_persistence_implementation(filename, precomputed_data_height, height)


def test_height_persistence_data_file_2():

    k = 1
    filename = file_names[k]
    compare_old_and_current_persistence_implementation(filename, precomputed_data_height, height)


def test_height_persistence_data_file_3():

    k = 2
    filename = file_names[k]
    compare_old_and_current_persistence_implementation(filename, precomputed_data_height, height)


def test_height_persistence_data_file_4():

    k = 0
    filename = file_names[k]
    compare_old_and_current_persistence_implementation(filename, precomputed_data_height, height)


### TEST PATH LENGTH PERSISTENCE ###


def test_path_length_persistence_data_file_1():

    k = 0
    filename = file_names[k]
    compare_old_and_current_persistence_implementation(filename, precomputed_data_path, path_length)

    
def test_path_length_persistence_data_file_2():

    k = 1
    filename = file_names[k]
    compare_old_and_current_persistence_implementation(filename, precomputed_data_path, path_length)


def test_path_length_persistence_data_file_3():

    k = 2
    filename = file_names[k]
    compare_old_and_current_persistence_implementation(filename, precomputed_data_path, path_length)


def test_path_length_persistence_data_file_4():

    k = 3
    filename = file_names[k]
    compare_old_and_current_persistence_implementation(filename, precomputed_data_path, path_length)


### TEST BRANCH ORDER PERSISTENCE ###


def test_branch_order_persistence_data_file_1():

    k = 0
    filename = file_names[k]
    compare_old_and_current_persistence_implementation(filename, precomputed_data_branch_order, branch_order)


def test_branch_order_persistence_data_file_2():

    k = 1
    filename = file_names[k]
    compare_old_and_current_persistence_implementation(filename, precomputed_data_branch_order, branch_order)


def test_branch_order_persistence_data_file_3():

    k = 2
    filename = file_names[k]
    compare_old_and_current_persistence_implementation(filename, precomputed_data_branch_order, branch_order)


def test_branch_order_persistence_data_file_4():

    k = 3
    filename = file_names[k]
    compare_old_and_current_persistence_implementation(filename, precomputed_data_branch_order, branch_order)
