import sys
import getopt
import os
import fnmatch
import pandas as pd
import neurontree.NeuronTree as nt
import computation.feature_presentation as fp
import computation.persistence_functions as pf


def main(argv):
    try:
        opts, args = getopt.getopt(argv, "c:f:d:", ["compute=", "file=", "dir=", "func="])
    except getopt.GetoptError:
        print('MorphoPy.py -c <compute_feature> [-f <swc_file> | -d <directory>]')
        sys.exit(2)

    # default values:
    compute = 'stats'  # default mode
    file = ''          # default no file -> directory is used
    swcdir = '.'       # default working directory
    function = None    # default function none

    # Arguments
    for opt, arg in opts:
        if opt in ('-c', '--compute'):
            compute = arg
        elif opt in ('-f', '--file'):
            file = arg
        elif opt in ('-d', '--dir'):
            swcdir = arg
        elif opt in ('--func'):
            function = arg

    # test if single file or directory and fill array with all files
    files = []
    if len(file) > 1:
        files.append(file)
        swcdir = "./"
    else:
        allfiles = os.listdir(swcdir)
        for f in allfiles:
            swcpattern = "*.swc"
            if fnmatch.fnmatch(f, swcpattern):
                files.append(f)

    ##### Compute morphometric statistics #####
    if compute == 'stats':
        print('##### Morphometric statistics #####')
        # process all files
        for file in files:
            print(" Process File: {} ".format(file))
            # compute morphometric stats
            try:
                swc = pd.read_csv(swcdir+file, delim_whitespace=True, comment='#',
                                  names=['n', 'type', 'x', 'y', 'z', 'radius', 'parent'], index_col=False)
                mytree = nt.NeuronTree(swc)
                print(fp.compute_Morphometric_Statistics(mytree))
            except:
                print("Failure in computing morphometric statistics!")

    ##### Compute persistence data #####
    elif compute == 'persistence':
        print('##### Persistence barcodes #####')
        # process all files
        for file in files:
            print(" Process File: {} ".format(file))
            # check if valid function selected
            if function in pf.functions:
                function = getattr(pf, function)
            else:
                function = None
            # compute persistence data
            try:
                swc = pd.read_csv(swcdir+file, delim_whitespace=True, comment='#',
                                  names=['n', 'type', 'x', 'y', 'z', 'radius', 'parent'], index_col=False)
                mytree = nt.NeuronTree(swc)
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