# Morphopy
(c) 2020 by Sophie Laturnus, Adam von Daranyi, Ziwei Huang and Philipp Berens.
********

## <a name="content">Content</a> 
- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Uninstallation](#uninstallation)
- [Usage](#usage)
- [Module description](#module)

## <a name="overview">Overview</a> 

MorphoPy is a Python3 package that uses networkX to compute and show information about neurites.
The input can be passed with single swc-files or it can handle whole directories with multiple files at once.
You can use MorphoPy imported in Python or from command line as batch-tool as well.

The current working build:
 
 **version 0.6**
 
All builds are tested on Linux (Debian and CentOS) and Windows 10.

You can find all working builds at [pypi.org](https://pypi.org/project/morphopy/).

[back to start](#content)

## <a name="requirements">Software Requirements</a>

In the following, all commands written in boxes need to be entered into your terminal.

### Mac:

 - homebrew (to install latest version of python3)
 - python >3.4: without homebrew you can find python [here](https://www.python.org/downloads/mac-osx/)
 - pip
 
**1)** If you want the latest python3 version you need to install homebrew,
       for that just run this command in your terminal:

    ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    
**2)** Then you can easily install the newest version of python3:

	brew install python

**3)** Now you should have python3 and pip3 installed. You can check the version like this:

	python3 --version
	pip3 -V

If both are enough you can continue with installing MorphoPy, see below.

### Linux:

 - python >3.4
 - pip 3 

Install Python3 on **Debian/Ubuntu** systems:

	apt-get install python3 python3-pip

Install Python3 on **RedHat/CentOS/Fedora systems** (perhaps python version differs):

	yum install python36 python3-pip

That's it. All requirements are met and you can continue with the MorphoPy installation. See below. 

### Windows:

 - python >3.5 (x64): download from [here](https://www.python.org/downloads/windows/)
 - pip : be sure you selected it during installation with the base python package
 - Microsoft Build Tools >14.0 you can download them [here](https://visualstudio.microsoft.com/visual-cpp-build-tools/) 

**1)** Install python with pip by executing the downloaded installation file and
       be sure to check the option to add python paths to enviroment variable at the first step!
       
**2)** Install Microsoft Build Tools with the automatic installation tool from the website above.
       Select C/C++ Compiler Tools at the selection page

All requirements are met now, see below for continue with MorphoPy installation.

[back to start](#content)

## <a name="installation">Installation:</a>

###  Install the MorphoPy package (all platforms):

Install the python package with pip3 and get the latest build:

	pip3 install morphopy

## <a name="uninstallation">Uninstallation:</a>

You can simply uninstall the package with pip3:

	pip3 uninstall morphopy

	
[back to start](#content)

## <a name="usage">Usage</a> 
Just call everywhere on the command line:

	MorphoPy.py <options>
Help:

	Usage: MorphoPy.py -c <compute_feature> [--wide | --func <persistence_function> | --conf <config_file>]
	                   [-f <swc_file> | -d <directory>] [-o <output directory>]
	
	Options:
	-c, --compute               parameter for selecting the computing feature:
	                            persistence: compute persistence data
                                stats      : compute morphometric statistics
                                density    : create density maps
       statistics options:
       --long                   you can change your output format, in long
                                format you get all values in a separate row.
                                (default: all values in one row)   
       persistence options:
       --func                   if persistence is selected as feature, you can
                                specify with this option a method function you
                                want to use at computing the persistence.
                                (default: radial distance function)
       density map options:
       --conf                   if density map is selected, you can pass a
                                config file with more parameters for creating
                                the density maps. (optional)
    -f, --file                   specifies a swc-file as input for Morphopy,
                                if no file is selected, directory is used
    -d, --directory              specifies a directory as input for swc-files.
                                (default: working directory)
    -o, --output                 specifies the output directory for saving the
                                results in. (default: same as source)


Available functions for persistence at the moment are:
 - radial_distance (default function)
 - height
 - path_length
 - branch_order

A sample config file for density maps looks like this (stored in a text file):

	[global]
	# specific distance for resampling nodes:
    distance: 1
    # width of each bin in microns across all dimensions
    #bin_size: 20
    # number of bins for each dimension (only if you don't use bin_size)
    n_bins_x: 20
    n_bins_y: 20
    n_bins_z: 20
    # if true: probabilty density is returned, count histogram otherwise
    density: True
    # smoothing the density data
    smooth: True
    # sigma used at smoothing
    sigma: 2
    # normalization bounds for density map:
    [norm_bound]
    r_max_x: 238.85
    r_max_y: 140.95
    r_max_z: 285.97
    r_min_x: -236.17
    r_min_y: -24.2
    r_min_z: -173.72

[back to start](#content)

## <a name="contributing">Contributing to MorphoPy </a>

We tested MorphoPy to the best of our knowledge and abilities in the scope of several projects. If you still find a bug
or you are missing a feature, please do not hesitate to contact us via [GitHub issues](https://github.com/berenslab/MorphoPy/issues).
Please try to provide a minimal example that reproduces the bug you want to be fixed.
If you want to develop the code base further, you can work with git pull requests. Please make sure that you document
the code and add tests and examples of how to use your code.


## <a name="citation"> Citing MorphoPy </a>

If you are using MorphoPy for your research or your work project please make sure to cite us and this repository.

## <a name="module">Module description</a> 

**Important:** MorphoPy requires the soma to be one single point. If several soma points are present in the
reconstruction file they will be automatically collapsed to the centroid of their convex hull on file loading (also see utils.get_standardized_swc).
MorphoPy currently only supports neurites that connect back to the soma. This means, axons that emerge from dendritic
structures can not be handled.

All data is stored in the [tidy data format](http://vita.had.co.nz/papers/tidy-data.pdf).

Please also refer to our [tutorial](https://github.com/berenslab/MorphoPy/blob/master/notebooks/MORPHOPY%20Tutorial.ipynb).

### Density maps
Density maps are marginal histograms over the neural mass. MorphoPy allows you to create density maps of different projections through the function compute_denisty_maps(). Per default it computes x, y, z, xy, xz and yz density maps from the point cloud of the original reconstruction. The point cloud is constructed through resampling along all neurites with a default distance of 1 micron. The resulting point clous is then binned into bins of 20 microns and smoothed using Gaussian smoothing with std of 1.

However, you can customize all these parameters by passing a config file to the function (see LINK).


### Morphometric statistics
Available morphometric statistics are:
- branch orders
- Strahler order
- branch angles
- path angles
- root angles
- thickness
- segment lengths
- path length to soma
- radial distance


Changing the Code

### Morphometric distributions
These features are not (yet) available via the command line tool.
Frequency histograms or Gaussian kernel density estimates can be queried for all single value statistics using the
methods `get_histogram(key)` or `get_kde_distribution(key)`. If you provide a distance measure (e.g. branch order,
path distance from soma or radial distance) the returned distribution will be two-dimensional and allows to investigate
a features' development across space.

### Persistence
changing distance function. Adding a distance function


### Not enough? ###

You want to compute your own features? Go for it! We recommend you to check out `networkx` and `shapely` for more options.


[back to start](#content)

