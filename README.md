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
 
 **version 0.4**
 
All builds are tested on Linux (Debian and CentOS) and Windows 10.

You can always find the latest working build in the "/dist" directory of this repository.

[back to start](#content)


In the following, all commands written in boxes need to be entered into your terminal.
## <a name="requirements">Software Requirements</a>
### Mac:
You need python with pip installed and graphviz/pygraphviz:

 - python >3.4
 - pip
 - graphviz
 - pygraphviz
 
**1)** Use homebrew to install, if it is not installed you can do it from your terminal with this command:

    ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    
**2)** Then you can easily install python and graphviz like this:

	brew install python graphviz
	
	# if you get an error that graphviz can't be found you can use this to install:
	pip install pygraphviz –install-option=”–include-path=/usr/local/lib/graphviz/” –install-option=”–library-path=/usr/local/lib/graphviz”
	
**3)** In the last step you need to build pygraphyviz with the installed graphviz library:

	pip install pygraphviz
    
Now you should be able to continue with the MophoPy installation. See below.
### Linux:
Debian/Ubuntu packages:

 - python 3
 - pip 3
 - graphviz development
 - geometry library development

You can install them with this command:

	apt-get install python3 python3-pip libgeos-dev libgraphviz-dev

RedHat/CentOS/Fedora packages:

 - python >3.6
 - python development >3.6
 - gcc-c++ compiler
 - graphviz development
 - geometry library development

You can install them with this command (perhaps python version differs):

	yum install gcc-c++ python36 python36-devel geos-devel graphviz-devel

That's it. All requirements are met and you can continue with the Morphopy installation. See below. 
### Windows:

 - python >3.5 (x64): download from [here](https://www.python.org/downloads/windows/)
 - pip >18 : be sure you selected it during installation with the base python package
 - Microsoft Build Tools >14.0 you can download them [here](https://visualstudio.microsoft.com/visual-cpp-build-tools/)
 - Graphviz (x64) [download here](https://github.com/mahkoCosmo/GraphViz_x64/)
 - Pygraphviz (x64) [from here](https://github.com/pygraphviz/pygraphviz/releases)
 - patched pygraphviz files from [Kagami@Pygraphviz](https://github.com/Kagami/pygraphviz/commit/fe442dc16accb629c3feaf157af75f67ccabbd6e)

**1)** Install python with pip and be sure to check the option to add python paths to enviroment variable.

**2)** Install Microsoft Build Tools with the automatic installation tool from the website above

**3)** Now download and unpack graphviz(x64), pygraphviz(x64) and the patched files from Kagami to any folder.

**4)** Then replace in pygraphviz the files with the two patched ones:

       ...\Kagami\pygraphviz\graphviz.i -> ...\pygraphviz-1.x\pygraphviz\graphviz.i
       ...\Kagami\pygraphviz\graphviz_wrap.c -> ...\pygraphviz-1.x\pygraphviz\graphviz_wrap.c

**5)** Now you can use this command in the base directory of pygraphviz-1.x to build it (use the full path to the directories of the unpacked graphviz files):

	python setup.py install --user --include-path="C:\\...\\graphviz\\include" --library-path="C:\\...\\graphviz\\lib"

Pygraphviz is now build without any errors and all requirements are met, see below for continue with MorphoPy installation.

[back to start](#content)
## <a name="installation">Installation:</a>

###  Install the MorphoPy package (all platforms):

Install the python package with pip/pip3 and get the latest build:

	pip install morphopy

## <a name="uninstallation">Uninstallation:</a>

You can simply uninstall the package with pip/pip3:

	pip uninstall morphopy

	
[back to start](#content)

## <a name="usage">Usage</a> 
Just call everywhere on the command line:

	MorphoPy.py <options>
Help:

	Usage: MorphoPy.py -c <compute_feature> [--func <persistence_function> | --conf <config_file>]
	                   [-f <swc_file> | -d <directory>] [-o <output directory>]

	Options:
	   -c, --compute                parameter for selecting the computing feature:
	                                persistence: Compute persistence data
	                                stats      : Compute Morphometric statistics
	                                density    : Create density maps
	       Persistence Options:
	       --func                   if persistence is selected as feature, you can
	                                specify with this option a method function you
	                                want to use at computing the persistence.
	       Density Map Options:
	       --conf                   if density map is selected, you can pass a
	                                config file with more parameters for creating
	                                the density maps. (optional)
	   -f, --file                   specifies a swc-file as input for Morphopy,
	                                if no file or directory is selected, working
	                                directory is used as default.
	   -d, --directory              specifies a directory as input for swc-files.
	                                (default: working directory)
	   -o, --output                 specifies the output directory for saving the
	                                results in. (default: same as source)



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

All data is stored in the tidy data format (http://vita.had.co.nz/papers/tidy-data.pdf).

Please also refer to our [tutorial](./notebooks/MORPHOPY Tutorial.ipynb).

### Density maps
Changing config

### Persistence
changing distance function. Adding a distance function

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


### Not enough? ###

You want to compute your own features? Go for it! We recommend you to check out `networkx` and `shapely` for more options.


[back to start](#content)

