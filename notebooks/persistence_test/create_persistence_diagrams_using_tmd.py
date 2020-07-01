#!bin python2

import tmd
import os

#files = ['eNGC-j140908b_cell1.swc', 'C4.swc', 'ds_1_cell_390.swc']
#files = ['C010398B-P2.CNG.swc']

# get files
root, _ , files = list(os.walk("./data/"))[0]

# define filter functions
features = ['radial_distances', 'path_distances', 'projection', 'section_branch_orders']

for file in files:
    n = tmd.io.load_neuron("./data/" + file)

    for f in features:
        ph = tmd.methods.get_ph_neuron(n, feature=f)
        tmd.methods.write_ph(ph, "./diagrams/%s_%s.txt"% (file.split(".")[0], f))
