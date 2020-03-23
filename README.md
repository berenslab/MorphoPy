# Morphopy
(c) 2020
********

## <a name="content">Content</a> 
- [Overview](#overview)
- [Installation](#installation)
- [Uninstallation](#uninstallation)
- [Usage](#usage)
- [Module description](#module)

## <a name="overview">Overview</a> 

[back to start](#content)

## <a name="installation">Installation:</a>
In the following, all commands written in boxes need to be entered into your terminal.
### Software Requirements
#### Mac:

#### Linux:
Debian/Ubuntu packages:

 - git
 - python 3
 - pip 3
 - graphviz development
 - geometry library development

You can install them with this command:

	apt-get install python3 python3-pip git libgeos-dev libgraphviz-dev

RedHat/CentOS/Fedora packages:

 - python 3.6
 - pip 3.6
 - python 3.6 development
 - gcc-c++ compiler
 - graphviz development
 - geometry library development

You can install them with this command:

	yum install git gcc-c++ python36 python36-pip python36-devel geos-devel graphviz-devel
#### Windows:

 - python 3.7
 - pip 3.7
 - graphviz development ???
 - geometry library development ???
###  Install the python package (Mac & Linux)

**1)** Clone the GitHub repository into your local directory:

	git clone https://github.com/berenslab/MorphoPy
Note: If the git command does not work under Mac OS, you have to install the missing package via xcode.

**2)** Now the last step: install the package from the new created GitHub folder. And make sure you are using pip3 (install it if it is missing):

	pip3 install MorphoPy/dist/morphopy-*.whl
	
The installation is finished and MorphoPy should be available from everywhere on the command prompt.
	
[back to start](#content)

## <a name="uninstallation">Uninstallation:</a>

You can simply uninstall the package with pip3:

	pip3 uninstall morphopy

	
[back to start](#content)

## <a name="usage">Usage</a> 
Just call everywhere on the command line:

	MorphoPy <options>
Help:

	Usage: MorphoPy -c <compute_feature> [--func <persistence_function> | --conf <config_file>]
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

## <a name="module">Module description</a> 

[back to start](#content)
