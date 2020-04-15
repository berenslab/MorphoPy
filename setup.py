#!/usr/bin/python3

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
      name='morphopy',
      version='0.2',
      scripts=['MorphoPy.py'] ,
      author="Adam von Daranyi",
      author_email="adam.von-daranyi@uni-tuebingen.de",
      description="morphology command line tool",
      long_description=long_description,
      long_description_content_type="text/markdown",
      url="https://github.com/berenslab/MorphoPy",
      packages=setuptools.find_packages(),
      classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
      ],
      install_requires=[
         'matplotlib',
         'scipy',
         'numpy',
         'pandas',
         'networkx',
         'seaborn',
         'sklearn',
         'shapely',
         'pygraphviz'
      ],
)
