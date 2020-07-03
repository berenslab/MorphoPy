from morphopy.computation.file_manager import load_swc_file
from morphopy.computation import feature_presentation as fp

import pandas as pd
from numpy import round
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
precomputed_morphometric_statistics = pd.read_csv(dir_path + "/test_data/precomputed_data_stats.csv", index_col='filename')

file_names = ['Image001-005_01_CNenhance_latest_LXS.swc', 'C4.swc',
              'V1_Layer23_Chat-IRES-Cre-neo_Ai14-299537.04.02.01_614430666_m.swc', 'EC3-80604.CNG.rotated.swc']

morphometric_statistic_names = ['avg_thickness', 'max_segment_path_length', 'max_path_dist_to_soma',
       'max_thickness', 'log_max_tortuosity', 'tree_asymmetry', 'depth',
       'width', 'log_median_tortuosity', 'total_length', 'min_path_angle',
       'log_min_tortuosity', 'median_path_angle', 'branch_points',
       'median_intermediate_segment_pl', 'max_branch_angle', 'stems',
       'median_terminal_segment_pl', 'min_branch_angle', 'max_branch_order',
       'mean_branch_angle', 'total_volume', 'max_degree', 'tips',
       'mean_soma_exit_angle', 'total_surface', 'height', 'max_path_angle']


def test_morphometrics_for_file_1():
    k = 0
    neuron = load_swc_file('../data/'+file_names[k])
    neuron_morphometrics = fp.compute_morphometric_statistics(neuron)

    for m in morphometric_statistic_names:
        old_value = round(precomputed_morphometric_statistics.loc[file_names[k]][m], 3)
        current_value = round(neuron_morphometrics[m].values[0], 3)
        assert (old_value == current_value), "Computation of %s changed in file %s!" % (m, file_names[k])


def test_morphometrics_for_file_2():
    k = 1
    neuron = load_swc_file('../data/'+file_names[k])
    neuron_morphometrics = fp.compute_morphometric_statistics(neuron)

    for m in morphometric_statistic_names:
        old_value = round(precomputed_morphometric_statistics.loc[file_names[k]][m], 3)
        current_value = round(neuron_morphometrics[m].values[0], 3)
        assert (old_value == current_value), "Computation of %s changed in file %s!" % (m, file_names[k])


def test_morphometrics_for_file_3():
    k = 2
    neuron = load_swc_file('../data/'+file_names[k])
    neuron_morphometrics = fp.compute_morphometric_statistics(neuron)

    for m in morphometric_statistic_names:
        old_value = round(precomputed_morphometric_statistics.loc[file_names[k]][m], 3)
        current_value = round(neuron_morphometrics[m].values[0], 3)
        assert (old_value == current_value), "Computation of %s changed in file %s!" % (m, file_names[k])


def test_morphometrics_for_file_4():
    k = 3
    neuron = load_swc_file('../data/'+file_names[k])
    neuron_morphometrics = fp.compute_morphometric_statistics(neuron)

    for m in morphometric_statistic_names:
        old_value = round(precomputed_morphometric_statistics.loc[file_names[k]][m], 3)
        current_value = round(neuron_morphometrics[m].values[0], 3)
        assert (old_value == current_value), "Computation of %s changed in file %s!" % (m, file_names[k])

