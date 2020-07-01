#!/usr/bin/python3

import setuptools
### load global MorphoPy variables
from morphopy import about

# get long description from readme
with open("README.md", "r", encoding='utf8') as f:
    long_description = f.read()

# get all requirements from txt file
with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setuptools.setup(
    name=about.__title__,
    version=about.__version__,
    entry_points={
        'console_scripts': ['morphopy=morphopy.MorphoPy:main'],
    },
    author=about.__author__,
    author_email=about.__email__,
    description=about.__summary__,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url=about.__url__,
    download_url="https://github.com/berenslab/MorphoPy/archive/v%s.tar.gz" % about.__version__,
    packages=setuptools.find_packages(),
    classifiers=[
     "Programming Language :: Python :: 3",
     "License :: OSI Approved :: MIT License",
     "Operating System :: OS Independent",
    ],
    install_requires=requirements,
)
