# Data and Plotting Scripts for Birnstiel et al. 2015
This script loads the data and (re-)creates the figures for the paper

> T\. Birnstiel, S. M. Andrews, P. Pinilla, and M. Kama  
> *Dust Evolution Can Produce Scattered Light Gaps in Protoplanetary Disks*  
> Astrophysical Journal Letters, 2015

View the IPython notebook

- [html version with images](https://htmlpreview.github.io/?https://github.com/birnstiel/Birnstiel2015_scripts/blob/master/make_figures.html)
- [github rendered (without images)](https://github.com/birnstiel/Birnstiel2015_scripts/blob/master/make_figures.ipynb).

The data contained in this repository include output of the dust evolution code of [Birnstiel et al. 2010](http://adsabs.harvard.edu/abs/2010A%26A...513A..79B) and a [`RADMC-3D`](http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/) setup based upon the dust data. To run the entire notebook yourself, you need the `radmc3d` executable in your path, or you change `images_done = False` to `True`, in which case the images are not recalculated on the fly, but the pre-calculated images contained in the repository are used.

The `fits` images can also be directly downloaded below:

- [Image at 1.65 µm (~33 MB)](https://raw.githubusercontent.com/birnstiel/Birnstiel2015_scripts/master/radmc_data/data_disklifetime_mstar07_mratio1_rc200_vf10_alpha3_static_243_1e%2B06/paperimage_1.65micron.fits)
- [Image at 1.3 mm (~33 MB)](https://raw.githubusercontent.com/birnstiel/Birnstiel2015_scripts/master/radmc_data/data_disklifetime_mstar07_mratio1_rc200_vf10_alpha3_static_243_1e%2B06/paperimage_1300micron.fits)

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)