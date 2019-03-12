
<!-- README.md is generated from README.Rmd. Please edit that file -->
MGPMSimulations
===============

The package MGPMSimulations contains all scripts and data that have been created or generated for the simulation-based study in the article "Automatic Generation of Evolutionary Hypotheses using Mixed Gaussian Phylogenetic Models". Due to size constraints, this R-package is not intended for a release on CRAN. Instead, the readers of the article are encouraged to install or clone the package from its code-repository (currently github).

Installation
------------

This package contains large data files and is not intended for release on CRAN. The recommended way to use this package is to follow these steps:

-   Clone the package repository from github using the command (type the following three commands in your terminal / shell):

        git clone https://github.com/venelin/MGPMSimulations.git

-   Enable git large file system (LFS) in your working copy:

        git lfs install

-   Pull the binary data from the LFS

        git pull

-   \[optional: This step is only needed if you need to inspect the row result-files from the MGPM simulations.\] Concatenate the binary files ResultsDirs.tar.gz.part\* into one file. This binary file was split in parts using the command (`split -b 1400m ResultsDirs.tar.gz "ResultsDirs.tar.gz.part"`):

        cat ResultsDirs.tar.gz.parta* > ResultsDirs.tar.gz

-   -   \[optional: This step is only needed if you need to inspect the row result-files from the MGPM simulations.\] Decompress ResultsDirs.tar.gz:

            tar -zxvf ResultsDirs.tar.gz

-   Install the package from your local working copy (type the following command in your R interpreter). Note that the directory MGPMMammals should be in your current directory from which R is running:

``` r
install.packages("MGPMSimulations", repos=NULL, type="source")
```
