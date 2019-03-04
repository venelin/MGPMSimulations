
<!-- README.md is generated from README.Rmd. Please edit that file -->
MGPMSimulations
===============

The package MGPMSimulations contains all scripts and data that have been created or generated for the simulation-based study in the article "Automatic Generation of Evolutionary Hypotheses using Mixed Gaussian Phylogenetic Models". Due to size constraints, this R-package is not intended for a release on CRAN. Instead, the readers of the article are encouraged to install or clone the package from its code-repository (currently github).

Installation
------------

### Github

The latest version of the package can be installed using:

``` r
devtools::install_github("venelin/MGPMSimulations")
```

Cloning the repository
----------------------

Note that a traditional git clone of the package will fail to download all of the files containing the rww results from the simulations. These files have been stored using an extension of git for large binary files called git lfs. If you really wish to download these files (their total size is in the order of 10Gb), you should execute the following shell commands within the your cloned repository after doing a regular git clone:

``` shell
git lfs install
git pull
```
