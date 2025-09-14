#!/usr/bin/env python3
import fnmatch
import getopt
import os
import sys
import traceback

import matplotlib as mpl
import pandas as pd
import scipy.io as sio

mpl.use('Agg')
import matplotlib.pyplot as plt
from morphopy import about
from morphopy.computation import file_manager
import morphopy.computation.feature_presentation as fp
import morphopy.computation.persistence_functions as pf
# load global about variables
# set global version
__version__ = about.__version__

def help(exitcode=0):
    """
    Print help page and exit application depending on error passed

    :param exitcode: Errorcode which will be returned at Exit or 0 if no error occured
    """
    if exitcode == 0:
        version()
    print('')
    print('Usage: morphopy -c <compute_feature> -i <swc_file>|<directory> [-o <output directory>]')
    print('                   [--long | --func <persistence_function> | --conf <config_file>] [--help]')
    print('                                                                                           ')
    if exitcode == 0:
        print('Options:                                                                      ')
        print('   -h, --help                   show this help.                               ')
        print('                                                                              ')
        print('   -c, --compute                parameter for selecting the computing feature:')
        print('                                persistence: compute persistence data         ')
        print('                                stats      : compute morphometric statistics  ')
        print('                                density    : create density maps              ')
        print('       statistics options:                                                    ')
        print('       --long                   you can change your output format, in long    ')
        print('                                format you get all values in a separate row.  ')
        print('                                (default: all values in one row)              ')
        print('       persistence options:                                                   ')
        print('       --func                   if persistence is selected as feature, you can')
        print('                                specify with this option a method function you')
        print('                                want to use at computing the persistence.     ')
        print('                                (default: radial distance function)           ')
        print('       density map options:                                                   ')
        print('       --conf                   if density map is selected, you can pass a    ')
        print('                                config file with more parameters for creating ')
        print('                                the density maps. (optional)                  ')
        print('                                                                              ')
        print('   -i, --input                  specifies a swc-file or a directory as input  ')
        print('                                for morphopy (mandatory)                      ')
        print('                                                                              ')
        print('   -o, --output                 specifies the output directory for saving the ')
        print('                                results in. (default: same as source)         ')
    sys.exit(exitcode)


def version():
    """
    Print version of Morphopy to command line - nothing else will be done
    """
    print('MorphoPy version %s' % about.__version__)
    print('(c) %s' % about.__copyright__)
    print(about.__url__)


def printException(message="Unknown Error"):
    """
    Helper function for processing Exceptions in main function

    :param message: pass message for user defined Exception
    """
    tb = traceback.format_exc()
    tb = tb.split("\n", 3)[3];
    print()
    print(message)
    print()
    print("Exception occured:")
    print(tb)


def morphopy(argv):
    """
    Callable from command line or imported with several arguments
    see help page for more information

    :type argv: arguments passed for processing swc-files
    """
    try:
        opts, args = getopt.gnu_getopt(argv, "c:i:o:hv",
                                   ["compute=", "long", "func=", "conf=", "input=", "output=", "help", "version"])
    except getopt.GetoptError:
        print("Error: Wrong options are specified!")
        help(1)

    # check if arguments are empty
    if len(argv) < 1:
        print("Error: No arguments are used! At least the compute mode and input is needed.")
        help(1)

    # default values:
    compute = ''  # no compute mode selected
    input = ""  # input data as file or directory
    format = "wide"  # default is wide output format for stats
    function = None  # default function none
    output = None  # default output directory is none
    configfile = None  # default value, no config file is used

    for arg in args:
        if arg in ['stats', 'density', 'persistence']:
            compute = arg

    # Check arguments
    for opt, arg in opts:
        # argument is missing because next option is in argument
        if arg.startswith("-"):
            print('Error: Wrong argument in option: %s'%opt)
            help(1)
        if opt in ('-c', '--compute'):
            compute = arg
        elif opt in ('-i', '--input'):
            input = arg
        elif opt in ('-o', '--output'):
            output = arg
        elif opt in '--long':
            format = "long"
        elif opt in '--conf':
            configfile = arg
        elif opt in '--func':
            # check if valid function selected and set pointer
            # if unknown, computing without special function
            if arg in pf.functions:
                function = getattr(pf, arg)
        elif opt in ('-v', '--version'):
            version()
            exit(0)
        elif opt in ('-h', '--help'):
            help(0)

    # if single file or directory, fill array with all files and set directory of data
    allfiles = []
    directory = ""
    if os.path.isfile(input):
        allfiles.append(os.path.basename(input))
        directory = os.path.dirname(input)
        if len(directory) == 0:
            directory = "./"
    elif os.path.isdir(input):
        directory = input
        allfiles = os.listdir(input)

    # check directory for trailing slash
    if not directory.endswith("/"):
        directory = "%s/" % directory

    # if no output directory is specified use source directory and check for trailing slash
    if output is None:
        output = directory
    elif not output.endswith("/"):
        output = "%s/" % output

    # test if files have valid extension
    files = []
    swc_ext = "*.swc"
    for f in allfiles:
        if fnmatch.fnmatch(f, swc_ext):
            files.append(f)

    # no valid files found
    if len(files) < 1:
        print('Error: No input is specified! Please use a valid file or directory.')
        help(1)

    # save output data in a dataframe with all proceeded files for export
    output_data = pd.DataFrame()

    ##### Compute morphometric statistics #####
    if compute == 'stats':
        print('##### Morphometric statistics #####')
        # process all files
        for file in files:
            print('... Process File: ', file)
            # compute morphometric stats
            try:
                # import swc file, compute statistics and create an output
                mytree = file_manager.load_swc_file(directory + file)
                morpho_stats_frame = fp.compute_morphometric_statistics(mytree, format=format)
                morpho_stats_frame['filename'] = file
                print(morpho_stats_frame)
                print("")
                # Append stats to global dataframe
                output_data = output_data.append(morpho_stats_frame)
            except:
                printException("Failure in computing morphometric statistics!")

    ##### Compute persistence data #####
    elif compute == 'persistence':
        print('##### Persistence barcodes #####')
        # process all files
        for file in files:
            print('... Process File: ', file)
            # compute persistence data
            try:
                # import swc file, compute persistence table and give an output
                mytree = file_manager.load_swc_file(directory + file)
                persistence_frame = fp.get_persistence(mytree.get_topological_minor(), f=function)
                persistence_frame['filename'] = file
                print(persistence_frame)
                print()

                # Append persistence to global dataframe
                output_data = output_data.append(persistence_frame)

                # create diagram and save png-file to same directory
                x = persistence_frame['birth']
                y = persistence_frame['death']
                plt.figure()
                plt.scatter(x, y, alpha=0.5)
                if function is None:
                    plt.title('Persistence Diagram')
                else:
                    plt.title('Persistence Diagram (%s)' % function.__name__)
                plt.xlabel('Birth')
                plt.ylabel('Death')
                outputfile = '%s%s_persistence.png' % (output, file)
                plt.savefig(outputfile)
                print('Persistence saved to: ', outputfile)
                plt.close()
            except:
                printException("Failure in computing persistence data!")

    ##### Compute density map #####
    elif compute == 'density':
        print('##### Density Map #####')
        # process all files
        for file in files:
            print('... Process File: %s' % file)
            # compute density map
            try:
                # import swc file and initialize NeuronTree
                mytree = file_manager.load_swc_file(directory + file)
                # read configfile if available, else set dicts None
                config_params = file_manager.read_config(configfile)
                # compute density maps and use config if available
                densities = fp.compute_density_maps(mytree, config_params)
                # plot densities
                plot = fp.plot_density_maps(densities)
                # get bins from density keys
                x_bins = densities['x_proj']['bins'][0]
                y_bins = densities['y_proj']['bins'][0]
                z_bins = densities['z_proj']['bins'][0]
                bins = 'x%iy%iz%i' % (x_bins, y_bins, z_bins)
                # build output path and save plots there
                outputfile = '%s%s_density_%s.png' % (output, file, bins)
                # plot.savefig('%s_%s.png' % (outputfile, bins))
                plot.savefig(outputfile)
                plt.close()
                print('Density maps plotted to %s' % outputfile)
                # build output path and save density data
                outputfile = '%s%s_density.mat' % (output, file)
                sio.savemat(outputfile, densities)
                print('Density map data saved to %s' % outputfile)
            except:
                printException("Failure in computing density map!")
    else:
        print('Error: Unknown compute mode. Use a valid compute parameter')
        help(1)

    if len(output_data) > 0:
        # output all data in one file named after directory(multiple files) or single file
        if len(files) > 1:
            filename = os.path.basename(os.path.normpath(directory))
        else:
            filename = files[0]
        output_filename = '%s%s_%s.csv' % (output, filename, compute)
        output_data.to_csv(output_filename, index=False)
        print('Data saved to: %s' % output_filename)

    # program end
    sys.exit(0)


def main():
    """
    callable from console script without arguments
    """
    argv=sys.argv[1:]
    print(argv)
    morphopy(argv)


if __name__ == '__main__':
    """
    callable from command line
    """
    main()

