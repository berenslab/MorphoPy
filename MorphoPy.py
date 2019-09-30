import sys
import getopt
import os
import fnmatch
import pandas as pd
import networkx as nx
import neurontree.NeuronTree as nt
import computation.feature_presentation as fp
import computation.persistence_functions as pf


def main(argv):
    try:
        opts, args = getopt.getopt(argv, "c:f:d:", ["compute=", "file=", "dir=", "func="])
    except getopt.GetoptError:
        print('Please use: MorphoPy.py -c <compute_feature> [-f <swc_file> | -d <directory>]')
        sys.exit(2)

    # default values:
    compute = 'stats'  # default mode
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
        print('Please use: MorphoPy.py -c <compute_feature> [-f <swc_file> | -d <directory>]')
        sys.exit(2)

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
                print(fp.compute_Morphometric_Statistics(mytree))
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
                print(fp.get_persistence(mytree, f=function))
            except:
                print("Failure in computing persistence data!")

    else:
        print('Unknown mode. Use a valid compute parameter')
        sys.exit(2)


if __name__ == '__main__':
    '''
    callable from command line
    Can open data of the following file format:
    - .swc
    '''
    print(sys.argv)
    main(sys.argv[1:])
