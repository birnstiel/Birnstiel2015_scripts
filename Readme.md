# Data and Plotting Scripts for Birnstiel et al. 2015

## General

This script loads the data and (re-)creates the figures for the paper

> T\. Birnstiel, S. M. Andrews, P. Pinilla, and M. Kama  
> *Dust Evolution Can Produce Scattered Light Gaps in Protoplanetary Disks*  
> Astrophysical Journal Letters, 2015

View the IPython notebook

- [html version with images](https://htmlpreview.github.io/?https://github.com/birnstiel/Birnstiel2015_scripts/blob/master/make_figures.html)
- [github rendered (without images)](https://github.com/birnstiel/Birnstiel2015_scripts/blob/master/make_figures.ipynb).

The data contained in this repository include output of the dust evolution code of [Birnstiel et al. 2010](http://adsabs.harvard.edu/abs/2010A%26A...513A..79B) and a [`RADMC-3D`](http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/) setup based upon the dust data. To run the entire notebook yourself, you need the `radmc3d` executable in your path, or you change `images_done = False` to `True`, in which case the images are not recalculated on the fly, but the pre-calculated images contained in the repository are used.

## Images and Figures

The figures can be downloaded directly as pdf:

- [Figure 1](https://raw.githubusercontent.com/birnstiel/Birnstiel2015_scripts/master/fig1.pdf)
- [Figure 2](https://raw.githubusercontent.com/birnstiel/Birnstiel2015_scripts/master/fig2.pdf)
- [Figure 3](https://raw.githubusercontent.com/birnstiel/Birnstiel2015_scripts/master/fig3.pdf)
- [Figure 4](https://raw.githubusercontent.com/birnstiel/Birnstiel2015_scripts/master/fig4.pdf)

The `fits` images can also be directly downloaded below:

- [Image at 1.65 µm (~33 MB)](https://raw.githubusercontent.com/birnstiel/Birnstiel2015_scripts/master/radmc_data/data_disklifetime_mstar07_mratio1_rc200_vf10_alpha3_static_243_1e%2B06/paperimage_1.65micron.fits)
- [Image at 1.3 mm (~33 MB)](https://raw.githubusercontent.com/birnstiel/Birnstiel2015_scripts/master/radmc_data/data_disklifetime_mstar07_mratio1_rc200_vf10_alpha3_static_243_1e%2B06/paperimage_1300micron.fits)

[![DOI](https://zenodo.org/badge/1015/birnstiel/Birnstiel2015_scripts.svg)](https://zenodo.org/badge/latestdoi/1015/birnstiel/Birnstiel2015_scripts) [![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/) 

## File Descriptions

| Name | Description |
| -----|-----|
| `make_figures.ipynb`              | the ipython notebook to create the figures |
| `make_figures.html`               | a html rendered version of the ipython notebook | 
| `Readme.md`                       | description of the repository
| `aux_functions.py`                | constants and auxiliary functions for plotting or image processing |
| `data.hdf5`                       | dust simulation parameters and results |
| `radmc_data/`                     | folder containing the `RADMC-3D` setup and images |
| `distribution_reconstruction.py`  | function to reconstruct the size distribution according to the paper |
| `fig1.py`                         | plotting function to create Figure 1
| `fig*.pdf`                        | resulting figures used in the paper |
