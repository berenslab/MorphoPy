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
- ADD INTRO

`MorphoPy` is meant to easily extract different feature representations from neural morphologies for downstream
statistical analysis. It features representations such as density maps, single valued summary statistics (morphometrics),
morphometric distributions and persistence images[@li:2017; @kanari:2018]. For some
of the individual representations there has been code published already [@cuntz:2011; @li:2017; @kanari:2018; @NeuroM]
but `MorphoPy` bundles these most common representations in one simple open source framework to make them accessible
to a large community. It can be used twofold either as a standalone command line tool or as a Python package within a
scientific workflow.

Persistence images can also be
been written in different languages (MATLAB, C, Python2, Python3) and only they cover one or two types of representation
each.

`MorphoPy` has been developed in the context of a benchmarking study for cortical interneuron cell type classification
based on their morphology [@laturnus:2019]. It has already been used in a series of scientific publications that tried
to relate transcriptome, electrophysiology and morphology of cortical interneurons in V1 and M1 [@scala:2019; @scala:2020].

TODO
Add brief summary of functionality

# Acknowledgements
We want to acknowledge contributions from Ziwei Huang in the beginning of the project.

This work was funded by the German Ministry of Education and Research(FKZ 01GQ1601), the German Research Foundation (DFG)
under Germany’s Excellence Strategy (EXC2064/1 – 390727645; BE5601/4-1, SFB 1233 “Robust Vision”, Project number 276693517),
and the National Institute of Mental Health (U19MH114830).

# References