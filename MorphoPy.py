import sys
import getopt
import os
import fnmatch
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import neurontree.NeuronTree as nt
import computation.feature_presentation as fp
import computation.persistence_functions as pf


def help():
    print(
        'Please use: MorphoPy.py -c <compute_feature> [--func <persistence_function>] [-f <swc_file> | -d <directory>]')
    sys.exit(2)


def main(argv):
    try:
        opts, args = getopt.getopt(argv, "c:f:d:", ["compute=", "file=", "dir=", "func="])
    except getopt.GetoptError:
        help()

    # default values:
    compute = ''       # no compute mode selected
    directory = './'   # default working directory if no file and dir specified
    file = ""          # default no file -> directory is used
    function = None    # default function none

    # Check arguments
    for opt, arg in opts:
        if opt in ('-c', '--compute'):
            compute = arg
        elif opt in ('-f', '--file'):
            file = arg
        elif opt in ('-d', '--dir'):
            directory = arg
        elif opt in '--func':
            # check if valid function selected and set pointer
            # if unknown, computing without special function
            if arg in pf.functions:
                function = getattr(pf, arg)

    # if single file or directory, fill array with all files
    allfiles = []
    if len(file) > 1:
        allfiles.append(file)
        directory = ""
    else:
        allfiles = os.listdir(directory)

    # test if files have valid extension
    files = []
    swc_ext = "*.swc"
    for f in allfiles:
        if fnmatch.fnmatch(f, swc_ext):
            files.append(f)

    # no valid files found
    if len(files) < 1:
        print("No valid file is specified or no file found in current directory!")
        help()

    # set version of networkX
    nxversion = 1
    if float(nx.__version__) >= 2:
        nxversion = 2

    ##### Compute morphometric statistics #####
    if compute == 'stats':
        print('##### Morphometric statistics #####')
        # process all files
        for file in files:
            print(" Process File: {} ".format(file))
            # compute morphometric stats
            try:
                swc = pd.read_csv(directory+file, delim_whitespace=True, comment='#',
                                  names=['n', 'type', 'x', 'y', 'z', 'radius', 'parent'], index_col=False)
                mytree = nt.NeuronTree(swc=swc, nxversion=nxversion)
                morpho_stats_table = fp.compute_Morphometric_Statistics(mytree)
                print(morpho_stats_table)
                # Export of pandas data frame to csv-file in same directory
                morpho_stats_table.to_csv(directory+file+"_morpho_stats.csv")
            except:
                print("Failure in computing morphometric statistics!")

    ##### Compute persistence data #####
    elif compute == 'persistence':
        print('##### Persistence barcodes #####')
        # process all files
        for file in files:
            print(" Process File: {} ".format(file))
            # compute persistence data
            try:
                swc = pd.read_csv(directory+file, delim_whitespace=True, comment='#',
                                  names=['n', 'type', 'x', 'y', 'z', 'radius', 'parent'], index_col=False)
                mytree = nt.NeuronTree(swc=swc, nxversion=nxversion)
                persistence_table = fp.get_persistence(mytree.get_mst(), f=function)
                print(persistence_table)

                # create diagram and save png-file to same directory
                x = persistence_table['birth']
                y = persistence_table['death']
                plt.figure()
                plt.scatter(x, y, alpha=0.5)
                if function is None:
                    plt.title('Persistence Diagram')
                else:
                    plt.title('Persistence Diagram ({})'.format(function.__name__))
                plt.xlabel('birth')
                plt.ylabel('death')
                plt.savefig(directory+file+"_persistence.png")
                plt.close()
            except:
                print("Failure in computing persistence data!")

    ##### Compute density map #####
    elif compute == 'density':
        print('##### Density Map #####')
        # process all files
        for file in files:
            print(" Process File: {} ".format(file))
            # compute density map
            try:
                swc = pd.read_csv(directory+file, delim_whitespace=True, comment='#',
                                  names=['n', 'type', 'x', 'y', 'z', 'radius', 'parent'], index_col=False)
                mytree = nt.NeuronTree(swc=swc, nxversion=nxversion)

                dist = 1  # in microns
                # get the resampled point could along each neurite at distance 1 micron.
                # pc is an array of 3D coordinates for each resampled node
                pc = mytree.resample_nodes(mytree.get_graph(), dist)
                plt.figure()
                plt.scatter(pc[:, 0], pc[:, 2], s=1)
                plt.title('Density Map')

                plt.savefig(directory+file+"_density.png")
                plt.close()
            except:
                print("Failure in computing density map!")
    else:
        print('Unknown compute mode. Use a valid compute parameter')
        help()


if __name__ == '__main__':
    '''
    callable from command line
    Can open data of the following file format:
    - .swc
    '''
    print(sys.argv)
    main(sys.argv[1:])
