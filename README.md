# BallMapper for Knot Theory data repository

This repository contains all data and scripts used in Dłotko, Gurnari and Sazdanovic "Mapper-type algorithms for complex data and relations".

The datasets we consider were originally created by J.S. Levitt [2].  

This repository contains the preprocessed data. Because of GitHub's limits on file sizes, the complete data, including original files from [2], can be found on [Zenodo](https://zenodo.org/records/7670819).


## Code/Software

Scripts are provided in the form of Jupyter notebooks. They are bases on the [python](https://github.com/dgurnari/pyBallMapper) implementation of the BallMapper algorithm that can be installed via

```
pip install pyballmapper
```

In addition to standard scientific packages, we use [Bokeh](https://bokeh.org) to generate interactive plots. 

## References

[1] Dłotko, Gurnari and Sazdanovic "Mapper-type algorithms for complex data and relations" 2023

[2] Levitt, Hajij and Sazdanovic "Big data approaches to knot theory: Understanding the structure of the Jones polynomial", 2022 https://doi.org/10.1142/S021821652250095X 
