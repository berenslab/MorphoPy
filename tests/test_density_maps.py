from morphopy.computation.file_manager import load_swc_file
from morphopy.computation.feature_presentation import compute_density_maps
from numpy import allclose, hstack
from scipy.io import loadmat
import os

dir_path = os.path.dirname(os.path.realpath(__file__))

file_names = ['Image001-005_01_CNenhance_latest_LXS.swc', 'C4.swc',
              'V1_Layer23_Chat-IRES-Cre-neo_Ai14-299537.04.02.01_614430666_m.swc', 'EC3-80604.CNG.rotated.swc']



def test_density_map_data_file_1():

    k=0
    neuron = load_swc_file('../data/'+file_names[k])
    old_density_maps = loadmat(dir_path+"/test_data/%s_density.mat"%file_names[k])

    current_density_map = compute_density_maps(neuron)

    for proj in ['x_proj', 'y_proj', 'z_proj']:

        old_data_values = old_density_maps[proj]['data'][0][0]
        current_data_values = current_density_map[proj]['data']

        assert allclose(old_data_values, current_data_values), \
            "Histogram values in %s have changed for file %s"%(proj, file_names[k])

        old_edge_values = old_density_maps[proj]['edges'][0][0]
        current_edge_values = current_density_map[proj]['edges'][0]

        assert allclose(old_edge_values, current_edge_values), \
            "Edge values in %s have changed for file %s" % (proj, file_names[k])

    for proj in ['xy_proj', 'xz_proj', 'yz_proj']:
        
        old_data_values = old_density_maps[proj]['data'][0][0]
        current_data_values = current_density_map[proj]['data']

        assert allclose(old_data_values, current_data_values), \
            "Histogram values in %s have changed for file %s" % (proj, file_names[k])

        old_edge_values = hstack(old_density_maps[proj]['edges'][0][0][0])
        current_edge_values = hstack(current_density_map[proj]['edges'])

        assert allclose(old_edge_values, current_edge_values), \
            "Edge values in %s have changed for file %s" % (proj, file_names[k])


def test_density_map_data_file_2():
    k = 1
    neuron = load_swc_file('../data/' + file_names[k])
    old_density_maps = loadmat(dir_path + "/test_data/%s_density.mat" % file_names[k])

    current_density_map = compute_density_maps(neuron)

    for proj in ['x_proj', 'y_proj', 'z_proj']:
        old_data_values = old_density_maps[proj]['data'][0][0]
        current_data_values = current_density_map[proj]['data']

        assert allclose(old_data_values, current_data_values), \
            "Histogram values in %s have changed for file %s" % (proj, file_names[k])

        old_edge_values = old_density_maps[proj]['edges'][0][0]
        current_edge_values = current_density_map[proj]['edges'][0]

        assert allclose(old_edge_values, current_edge_values), \
            "Edge values in %s have changed for file %s" % (proj, file_names[k])

    for proj in ['xy_proj', 'xz_proj', 'yz_proj']:
        old_data_values = old_density_maps[proj]['data'][0][0]
        current_data_values = current_density_map[proj]['data']

        assert allclose(old_data_values, current_data_values), \
            "Histogram values in %s have changed for file %s" % (proj, file_names[k])

        old_edge_values = hstack(old_density_maps[proj]['edges'][0][0][0])
        current_edge_values = hstack(current_density_map[proj]['edges'])

        assert allclose(old_edge_values, current_edge_values), \
            "Edge values in %s have changed for file %s" % (proj, file_names[k])


def test_density_map_data_file_3():
    k = 2
    neuron = load_swc_file('../data/' + file_names[k])
    old_density_maps = loadmat(dir_path + "/test_data/%s_density.mat" % file_names[k])

    current_density_map = compute_density_maps(neuron)

    for proj in ['x_proj', 'y_proj', 'z_proj']:
        old_data_values = old_density_maps[proj]['data'][0][0]
        current_data_values = current_density_map[proj]['data']

        assert allclose(old_data_values, current_data_values), \
            "Histogram values in %s have changed for file %s" % (proj, file_names[k])

        old_edge_values = old_density_maps[proj]['edges'][0][0]
        current_edge_values = current_density_map[proj]['edges'][0]

        assert allclose(old_edge_values, current_edge_values), \
            "Edge values in %s have changed for file %s" % (proj, file_names[k])

    for proj in ['xy_proj', 'xz_proj', 'yz_proj']:
        old_data_values = old_density_maps[proj]['data'][0][0]
        current_data_values = current_density_map[proj]['data']

        assert allclose(old_data_values, current_data_values), \
            "Histogram values in %s have changed for file %s" % (proj, file_names[k])

        old_edge_values = hstack(old_density_maps[proj]['edges'][0][0][0])
        current_edge_values = hstack(current_density_map[proj]['edges'])

        assert allclose(old_edge_values, current_edge_values), \
            "Edge values in %s have changed for file %s" % (proj, file_names[k])


def test_density_map_data_file_4():

    k = 3
    neuron = load_swc_file('../data/' + file_names[k])
    old_density_maps = loadmat(dir_path + "/test_data/%s_density.mat" % file_names[k])

    current_density_map = compute_density_maps(neuron)

    for proj in ['x_proj', 'y_proj', 'z_proj']:
        old_data_values = old_density_maps[proj]['data'][0][0]
        current_data_values = current_density_map[proj]['data']

        assert allclose(old_data_values, current_data_values), \
            "Histogram values in %s have changed for file %s" % (proj, file_names[k])

        old_edge_values = old_density_maps[proj]['edges'][0][0]
        current_edge_values = current_density_map[proj]['edges'][0]

        assert allclose(old_edge_values, current_edge_values), \
            "Edge values in %s have changed for file %s" % (proj, file_names[k])

    for proj in ['xy_proj', 'xz_proj', 'yz_proj']:
        old_data_values = old_density_maps[proj]['data'][0][0]
        current_data_values = current_density_map[proj]['data']

        assert allclose(old_data_values, current_data_values), \
            "Histogram values in %s have changed for file %s" % (proj, file_names[k])

        old_edge_values = hstack(old_density_maps[proj]['edges'][0][0][0])
        current_edge_values = hstack(current_density_map[proj]['edges'])

        assert allclose(old_edge_values, current_edge_values), \
            "Edge values in %s have changed for file %s" % (proj, file_names[k])