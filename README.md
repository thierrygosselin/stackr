# stackr

[![Travis-CI Build Status](https://travis-ci.org/thierrygosselin/stackr.svg?branch=master)](https://travis-ci.org/thierrygosselin/stackr)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/thierrygosselin/stackr?branch=master&svg=true)](https://ci.appveyor.com/project/thierrygosselin/stackr)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/stackr)](http://cran.r-project.org/package=stackr)
[![DOI](https://zenodo.org/badge/14548/thierrygosselin/stackr.svg)](https://zenodo.org/badge/latestdoi/14548/thierrygosselin/stackr)

This is the development page of the **stackr**, an R package who's goal is
to make GBS/RAD data easy to analyse inside R.

## Use stackr to:

* **Input**: various file formats are supported (VCF, PLINK tped/tfam, genind, genlight, genepop, haplotype file produced by [STACKS](http://catchenlab.life.illinois.edu/stacks/) and data frame of genotypes)
* **Explore**, **manipulate** and **visualize**: 
Caracteristics and statistics for important variables of GBS/RADseq data are available: read depth (coverage) of alleles and 
genotypes, genotype likelihood, genotyped individuals and populations, missing data, minor allele frequency (local and global MAF),
observed heterozygosity (Het obs) and inbreeding coefficient (Fis), duplicate individual, etc.
* **Filters**: avoid bad data exploration and control the impact of filters on your downstream genetic analysis. 
Alleles, genotypes, markers, individuals and populations can be filtered and/or selected in several ways.
* **Imputations**: Map-independent imputation of missing genotype/alleles 
using Random Forest or the most frequent category
* **Output**: painless conversion process. Conversion functions are integrated with important filters, blacklist and whitelist to convert genomic dataset in a **tidy** format and/or into *VCF*, *PLINK* , *genepop*, *genind*, *genlight*, *fstat*, *gtypes*, *betadiv* and *dadi* formats. 
Easy integration with other software or R packages like [adegenet] (https://github.com/thibautjombart/adegenet), [strataG] (https://github.com/EricArcher/strataG.devel/tree/master/strataG.devel), [hierfstat] (https://github.com/jgx65/hierfstat), [pegas] (https://github.com/emmanuelparadis/pegas), [poppr] (https://github.com/grunwaldlab/poppr) and [assigner](https://github.com/thierrygosselin/assigner)
* `ggplot2`-based plotting to view distributions of summary statistics and create publication-ready figures

 
## Installation
To try out the dev version of **stackr**, copy/paste the code below:

```r
# Install or load the package **devtools**
if (!require("devtools")) install.packages("devtools") # to install
library(devtools) # to load
# Install **stackr**
# devtools::install_github("thierrygosselin/stackr") # to install without vignettes
devtools::install_github("thierrygosselin/stackr", build_vignettes = TRUE)  # to install WITH vignettes
library(stackr) # to load
```

## Prerequisite - Suggestions - Troubleshooting
  * **Parallel computing**: Follow the steps in this [vignette](https://github.com/thierrygosselin/stackr/blob/master/vignettes/vignette_imputations_parallel.Rmd) 
  to install an OpenMP enabled [randomForestSRC](http://www.ccs.miami.edu/~hishwaran/rfsrc.html)
 package to do imputations in parallel.
  * **Installation problem:** see this
  [vignette](https://github.com/thierrygosselin/stackr/blob/master/vignettes/vignette_installation_problems.Rmd)
  * **Windows users**:  to have *stackr* run in parallel use [parallelsugar](https://github.com/nathanvan/parallelsugar).
  Easy to install and use ([instructions](https://github.com/nathanvan/parallelsugar#installation)).
  * For a better experience in **stackr** and in R in general, I recommend using [RStudio](https://www.rstudio.com/products/rstudio/download/). 
  The R GUI is unstable with functions using parallel ([more info](https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/mclapply.html)). 
  Below, the combination of packages and how I install/load them :
  
  ```r
  if (!require("pacman")) install.packages("pacman")
  library("pacman")
  pacman::p_load(devtools, reshape2, ggplot2, stringr, stringi, plyr, dplyr, tidyr, readr, purrr, data.table, ape, adegenet, parallel, lazyeval, randomForestSRC)
  pacman::p_load(devtools, reshape2, ggplot2, stringr, stringi, plyr, dplyr, tidyr, readr, purrr, data.table, ape, adegenet, parallel, lazyeval, randomForestSRC)
  # install_github("thierrygosselin/stackr", build_vignettes = TRUE) # uncomment to install
  library("stackr")
  ```

## Vignettes and examples

From a browser:
* [installation problems](https://github.com/thierrygosselin/stackr/blob/master/vignettes/vignette_installation_problems.Rmd)
* [parallel computing during imputations](https://github.com/thierrygosselin/stackr/blob/master/vignettes/vignette_imputations_parallel.Rmd) 
* [vcf2dadi](https://github.com/thierrygosselin/stackr/blob/master/vignettes/vignette_vcf2dadi.Rmd)
* [haplo2genind](https://github.com/thierrygosselin/stackr/blob/master/vignettes/vignette_haplo2genind.Rmd)
* [Missing data visualization and analysis](https://github.com/thierrygosselin/stackr/blob/master/vignettes/vignette_missing_data_analysis.Rmd)
Inside R:
```r
browseVignettes("stackr") # To browse vignettes
vignette("vignette_vcf2dadi") # To open specific vignette
```

Vignettes are in development, check periodically for updates.


## Citation:
To get the citation, inside R:
```r
citation("stackr")
```

## New features
Version, new feature and bug history now lives in the [NEWS.md file] (https://github.com/thierrygosselin/stackr/blob/master/NEWS.md)

**v.0.2.9**
* bug fix in `tidy_genomic_data`
* bug fix between stackr -> devtools -> github -> travis, [this page helped] (http://itsalocke.com/using-travis-make-sure-use-github-pat/)


For previous news:
[NEWS.md file] (https://github.com/thierrygosselin/stackr/blob/master/NEWS.md)

## Roadmap of future developments:

* Until publication **stackr** will change rapidly (see contributions below for bug reports).
* Updated filters: more efficient, interactive and visualization included.
* Better integration with other GBS/RADseq approaches, beside [STACKS](http://catchenlab.life.illinois.edu/stacks/).
* Integrated converter function to input and output several file formats.
* Workflow tutorial that links functions and points to specific vignettes to further explore some problems.
* Integration of several functions with [STACKS](http://catchenlab.life.illinois.edu/stacks/) and [DArT](http://www.diversityarrays.com) database.
* Use Shiny and ggvis when subplots or facets becomes available...
* Suggestions ?


## Contributions:

This package has been developed in the open, and it wouldn’t be nearly as good without your contributions. There are a number of ways you can help me make this package even better:  
* If you don’t understand something, please let me know. 
* Your feedback on what is confusing or hard to understand is valuable. 
* If you spot a typo, feel free to edit the underlying page and send a pull request.

New to pull request on github ? The process is very easy:  
* Click the edit this page on the sidebar.
* Make the changes using github’s in-page editor and save.
* Submit a pull request and include a brief description of your changes. 
* “Fixing typos” is perfectly adequate.


## GBS workflow
The **stackr** package fits currently at the end of the GBS workflow. Below, a flow chart using [STACKS] (http://catchenlab.life.illinois.edu/stacks/) and other software. You can use the [STACKS] (http://catchenlab.life.illinois.edu/stacks/) workflow [used in the Bernatchez lab] (https://github.com/enormandeau/stacks_workflow). ![](vignettes/GBS_workflow.png)

## stackr workflow 
Currently under construction. Come back soon!

**Table 1: Quality control and filtering RAD/GBS data**

| Parameter | Libraries/Seq.Lanes | Allele | Genotype | Individual | Sampling sites | Populations | Globally |
|:----|:----:|:----:|:----:|:----:|:----:|:----:|:----:|
| Quality | x | | | | | | |
| Coverage | | x | x | | | | |
| Genotype Likelihood | | | x | | | | |
| Prop. Genotyped | | | | x | x | x | x |
| MAF | | | | | x | x | x |
| HET | | | | | | x | |
| FIS | | | | | | x | |
| SNP number/reads | | | | | | | x |


Step 1 as a quality insurance step. We need to modify the data to play with it efficiently in R. To have reliable summary statistics, you first need good coverage of your alleles to call your genotypes, good genotype likelihood, enough individuals in each sampling sites and enough putative populations with your markers... 

Step 2 is where the actual work is done to remove artifactual and uninformative markers based on summary statistics of your markers.


