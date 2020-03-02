import sys
import getopt
import os
import fnmatch
import pandas as pd
import networkx as nx
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import neurontree.NeuronTree as nt
import computation.feature_presentation as fp
import computation.persistence_functions as pf


def help():
    print('Usage: MorphoPy.py -c <compute_feature> [--func <persistence_function> | --conf <config_file>]')
    print('                   [-f <swc_file> | -d <directory>] [-o <output directory>]')
    print('')
    print('Options:')
    print('   -c, --compute                parameter for selecting the computing feature:')
    print('                                persistence: Compute persistence data         ')
    print('                                stats      : Compute Morphometric statistics  ')
    print('                                density    : Create density maps              ')
    print('       Persistence Options:                                                   ')
    print('       --func                   if persistence is selected as feature, you can')
    print('                                specify with this option a method function you')
    print('                                want to use at computing the persistence.     ')
    print('       Density Map Options:                                                   ')
    print('       --conf                   if density map is selected, you can pass a    ')
    print('                                config file with more parameters for creating ')
    print('                                the density maps. (optional)                  ')
    print('   -f, --file                   specifies a swc-file as input for Morphopy,   ')
    print('                                if no file or directory is selected, working  ')
    print('                                directory is used as default.                 ')
    print('   -d, --directory              specifies a directory as input for swc-files. ')
    print('                                (default: working directory)                  ')
    print('   -o, --output                 specifies the output directory for saving the ')
    print('                                results in. (default: same as source)         ')

    sys.exit(2)


def main(argv):
    try:
        opts, args = getopt.getopt(argv, "c:f:d:o:", ["compute=", "func=", "conf=", "file=", "dir=", "output="])
    except getopt.GetoptError:
        print("Wrong options are specified!")
        help()

    # check if arguments are empty
    if len(argv) < 1:
      print("No arguments are used! At least the compute mode has to be passed.")
      help()

    # default values:
    compute = ''       # no compute mode selected
    directory = './'   # default working directory if no file and dir specified
    file = ""          # default no file -> directory is used
    function = None    # default function none
    output = None      # default output directory is none
    conf = None        # default value, no config file is used

    # Check arguments
    for opt, arg in opts:
        # argument is missing because next option is in argument
        if arg.startswith("-"):
            print("Wrong argument in option: "+opt)
            help()
        if opt in ('-c', '--compute'):
            compute = arg
        elif opt in ('-f', '--file'):
            file = arg
        elif opt in ('-d', '--dir'):
            directory = arg
        elif opt in ('-o', '--output'):
            output = arg
        elif opt in '--conf':
            conf = arg
        elif opt in '--func':
            # check if valid function selected and set pointer
            # if unknown, computing without special function
            if arg in pf.functions:
                function = getattr(pf, arg)

    # if single file or directory, fill array with all files
    allfiles = []
    if len(file) > 1:
        allfiles.append(os.path.basename(file))
        directory = os.path.dirname(file)+"/"
    else:
        allfiles = os.listdir(directory)

    # if no output directory is specified use source directory
    if output is None:
        output = directory

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
                morpho_stats_table.to_csv(output+file+"_morpho_stats.csv")
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
                plt.savefig(output+file+"_persistence.png")
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

                plots = fp.compute_Density_Maps(mytree, conf=conf)

                for idx, plot in enumerate(plots, start=1):
                    plot.savefig(output+file+"_density_%s.png" % idx)
                    plot.close()
            except FileNotFoundError:
                print("Failure in computing density map: Config file not found or not readable!")
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
