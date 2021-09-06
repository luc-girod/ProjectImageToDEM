---
title: 'Projimdem: A Python Package to Georefence and Project Images on Digital Elevation Models'
tags:
  - Python
  - photogrammetry
  - image processing
  - geosciences
  - glaciology
authors:
  - name: Simon Filhol^[co-first author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0003-1282-7307
    affiliation: 1
  - name: Luc Girod^[co-first author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0003-3627-5885
    affiliation: 1
affiliations:
 - name: Department of Geosciences, University of Oslo, Oslo, Norway
   index: 1
date: 6 Sepember 2021
bibliography: paper.bib
---

# Summary

[What,Why,How]


# Statement of need

Decribe why do we need this software: 
- Processing time-lapse imagery for studying snow cover extent, glaciers, etc.
- Need of a fully OS method as all equivalent software are Matlab dependent
- Flexibility of a Python implementation for extending the package as well as improving it within an OS environemnt


# Software Description

Describe the software steps and underlying mathematics

Exisiting software packages: 
- [PRACTISE](https://github.com/shaerer/PRACTISE) by `@gmd-9-307-2016` -> "H\"arer, S. et al., (2016)"
- [ImGRAFT](http://imgraft.glaciology.net/) by `@gi-4-23-2015` -> "Messerli, A. and Grinsted, A., (2015)"

# Quality Evaluation

Include short accuracy assement using a comparative test with ImGRAFT

# Software dependencies
- [gdal](https://gdal.org/python/)
- adapted and modified version of[photogrammetry-resection](https://github.com/jeffwalton/photogrammetry-resection)
- [Pandas](https://pandas.pydata.org/)
- [numpy](https://numpy.org)
- [scipy](https://docs.scipy.org)
- [matplotlib](https://matplotlib.org/)

Find citation for each software

# Acknowledgements

# References