---
title: 'MorphoPy: A python package for feature extraction of neural morphologies.'
tags:
    - Python
    - neurons
    - morphologies
    - anatomy
authors:
    - name: Sophie Laturnus
      orcid: 0000-0001-9532-788X
      affiliation: "1,3,4"

    - name: Adam von Daranyi
      affiliation: 4

    - name: Ziwei Huang
      affiliation: "1,3,4"

    - name: Philipp Berens
     orcid: 0000-0002-0199-4727
     affiliation: "1,2,3,4"

affiliations:
    - name: Institute for Ophthalmic Research, University of T\"ubingen, Germany
    index: 1
    - name: Institute for Bioinformatics and Medical Informatics, University of T\"ubingen, Germany
    index: 2
    - name: Bernstein Center for Computational Neuroscience, University of T\"ubingen, Germany
    index: 3
    - name: Center for Integrative Neuroscience, University of T\"ubingen, Germany
    index: 4
date: \today
bibliography: morphopy.bib
---


# Summary

Summarizing the anatomy of a neural cell type has been one of the oldest ways to describe differences in neural populations.
However, the quantitative analysis of neuronal morphologies persists to be a hard problem. It usually begins with choosing a
particular feature representation in order to make individual morphologies amenable to classical data analysis. Many
different feature representations have been suggested in the literature, such as density maps [jefferis:2007], single valued summary
statistics (morphometrics) [scorcioni:2008; @NeuroM] or, more recently, persistence images [@li:2017; @kanari:2018].
The tools for extracting them, however, are often focused on solely one such representation and scattered across various
programming languages.

Our software package `MorphoPy` is meant to easily extract different feature representations from neural morphologies for
downstream statistical analysis. It bundles common representations such as density maps, morphometrics, morphometric distributions
and persistence images in one simple open source framework to make them accessible to a large community.
`MorphoPy` can be used either as a standalone command line tool or as a Python package within a scientific workflow.

![Neural reconstructions are represented as direct acylic graphs with node and edge attributes.\label{fig:attributes}](./figures/Fig1_attributes_small.png)

`MorphoPy` builds on the functionality of the networkx package [@hagberg:2008] and represents each neuron as a directed
acylic tree graph with node and edge attributes \autoref{fig:attributes}. The package supports 2D plotting routines and the
generation of different feature representations.
```python
import MorphoPy.NeuronTree as nt
from MorphoPy.computation import file_manager as fm

N = fm.load_swc_file("../data/EC3-80604.CNG.swc")
Dendrites = N.get_dendritic_tree()

from neurontree.plotting import show_threeview
fig = plt.figure(figsize=(10,10))
show_threeview(N, fig)
```

![Plotting reconstructions in 2D. \label{fig:plot}](./figures/threeview_dendrites.png)

As shown in the code snippet above, it is also possible to split the reconstruction into its different parts (axon or dendrites only)
and operate on each neurite type separately.

Density maps are computed on the basis of a configuration file (or dictionary) that controls parameters such as bin size
and normalization ranges. Additionally, users can specify whether and to which degree
they want to smooth each density map \autoref{fig:dms}.
![XY-density map of the dendrite plotted above with different degrees of Gaussian smoothing. \label{fig:dms}](./figures/density_map_smoothing.png)

A variety of statistics can be computed on the nodes and edges of each reconstruction (\autoref{fig:morphometrics}).
The `get_morphometric_statistics()`-method offers a precompiled single valued selection of these statistics, but in principle,
they can be adjusted to the user's personal preference.
![Node and edge related morphometric statistics. \label{fig:morphometrics}](./figures/fig_morphometrics.png)
Additionally, it is possible to query the entire distribution of each statistic either in form of a histogram or as a
Gaussian kernel density estimate (kde). Fig. \autoref{fig:morphdist}, for example, shows the kde of radial distances, branch angles and their
combination for the dendrites shown in \autoref{fig:plot}.
![\label{fig:morphdist}](./figures/2D_morph_dist.png)

Furthermore, `MorphoPy` supports the generation of 2D persistence diagrams. Persistence diagrams are a fairly recent way of describing
the branching of neural morphologies [@li:2017; @kanari:2018] with respect to a specified distance function. By default,
`MorphoPy` computes a persistence diagram based on the radial distance from the soma, but users can choose from four different
pre-implemented distance functions (radial distance, path length, height or branch order) or provide their own.

```python
from computation.persistence_functions import radial_distance
from computation.feature_presentation import get_persistence

import numpy as np
def custom_distance(G, u, v):
    """
    Returns a distance between nodes u and v, which both are part of the graph given in G.
    """
    n = G.node[u]['pos']
    r = G.node[v]['pos']
    return np.dot(n, r)

df = get_persistence(Dendrites.get_topological_minor(), f=radial_distance)
df_custom = get_persistence(Dendrites.get_topological_minor(), f=custom_distance)
```

To make `MorphoPy` accessible to a non-programming audience it can be called from the command line to operate on
single files or entire batches.
```bash
MorphoPy.py -c [density|persistence|stats] -f ['path_to_file'] | -d ['path_to_folder']
```
For a full description of `MorphoPy`'s functionality please refer to our documentation and tutorial on our [gitHub page](https://github.com/berenslab/MorphoPy).

SUMMARIZING SENTENCE??

`MorphoPy` has been developed in the context of a benchmarking study for cortical interneuron cell type classification
based on their morphology [@laturnus:2019]. It has already been used in a series of scientific publications that tried
to relate transcriptome, electrophysiology and morphology of cortical interneurons in V1 and M1 [@scala:2019; @scala:2020].


# Acknowledgements

This work was funded by the German Ministry of Education and Research(FKZ 01GQ1601), the German Research Foundation (DFG)
under Germany’s Excellence Strategy (EXC2064/1 – 390727645; BE5601/4-1, SFB 1233 “Robust Vision”, Project number 276693517),
and the National Institute of Mental Health (U19MH114830).

# Reference