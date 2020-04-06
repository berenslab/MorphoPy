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
The tools for extracting them, however, are often focused on solely one such representation and written in various
programming languages.
Our software package `MorphoPy` is meant to easily extract different feature representations from neural morphologies for
downstream statistical analysis. It bundles common representations such as density maps, morphometrics, morphometric distributions
and persistence images in one simple open source framework to make them accessible to a large community. `MorphoPy`
can be used either as a standalone command line tool or as a Python package within a scientific workflow.

TODO
Add brief summary of functionality
Write about command line tool

```bash
MorphoPy -c [density|persistence|stats] -f ./data/C4.swc
```

`MorphoPy` builds on the functionality of the networkx package [CITATION] and represents each neuron internally as a
directed acylic tree graph with node and edge attributes (see Fig. 1a). One can also query dendrites, the axon or all
individual neurites extending from the soma separately.
The package allows you to query a variety of morphometric values and distributions (see Fig. 1),


Each node within this graph is assigned a type-label that determines whether a node belongs
to the soma, the axon or the (apical) dendrites. Additionally each node is


For a full description of the packages functionality please refer to our documentation and tutorial

`MorphoPy` has been developed in the context of a benchmarking study for cortical interneuron cell type classification
based on their morphology [@laturnus:2019]. It has already been used in a series of scientific publications that tried
to relate transcriptome, electrophysiology and morphology of cortical interneurons in V1 and M1 [@scala:2019; @scala:2020].


# Acknowledgements
We want to acknowledge contributions from Ziwei Huang in the beginning of the project.

This work was funded by the German Ministry of Education and Research(FKZ 01GQ1601), the German Research Foundation (DFG)
under Germany’s Excellence Strategy (EXC2064/1 – 390727645; BE5601/4-1, SFB 1233 “Robust Vision”, Project number 276693517),
and the National Institute of Mental Health (U19MH114830).

# Reference