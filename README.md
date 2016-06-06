# stackr

[![Travis-CI Build Status](https://travis-ci.org/thierrygosselin/stackr.svg?branch=master)](https://travis-ci.org/thierrygosselin/stackr)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/thierrygosselin/stackr?branch=master&svg=true)](https://ci.appveyor.com/project/thierrygosselin/stackr)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/stackr)](http://cran.r-project.org/package=stackr)
[![DOI](https://zenodo.org/badge/14548/thierrygosselin/stackr.svg)](https://zenodo.org/badge/latestdoi/14548/thierrygosselin/stackr)

The goal of **stackr** is to make GBS/RAD data produced by [STACKS] (http://catchenlab.life.illinois.edu/stacks/) and other pipelines easy to analyse in R.

This is the development page of the **stackr** package for the R software, optimized for *de novo* and population genetics.

## Use stackr to:

* Avoid bad data exploration. Use different filters and convert to the appropriate data file format.
* Filter genetic markers based on: alleles and 
genotype coverage (read depth), genotype likelihood, proportion/percentage/number 
of genotyped individuals and populations, minor allele frequency (local and global),
observed heterozygosity and inbreeding coefficient (Fis)
* Read and modify several output files: VCF, PLINK (tped/tfam), haplotype file produced by [STACKS] (http://catchenlab.life.illinois.edu/stacks/), genepop, genind and dataframes in wide or long/tidy formats
* Transform genomic file format into a tidy format useful to quickly visualise and filter summary statistics within R
* `ggplot2`-based plotting to view distributions of summary statistics and create publication-ready figures
* Convert data into *genepop*, *genind*, *fstat*, *gtypes*, *betadiv* and *dadi* files or objects for easy integration with other software or R packages like [adegenet] (https://github.com/thibautjombart/adegenet), [strataG] (https://github.com/EricArcher/strataG.devel/tree/master/strataG.devel), [hierfstat] (https://github.com/jgx65/hierfstat), [pegas] (https://github.com/emmanuelparadis/pegas) and [poppr] (https://github.com/grunwaldlab/poppr)
* Map-Independent Imputation of missing genotype or allele using Random Forest within several functions: *haplo2genepop*, *haplo2genind*, *haplo2hierfstat*, *haplo2gtypes*, *haplo2colony*, *vcf2genind*, *vcf2hierfstat*, *vcf2betadiv*, *vcf2dadi* and *vcf_imputation*. 

## Installation
You can try out the dev version of **stackr**. Follow the 3 steps below:

**Step 1:** Install or load the package **devtools**
```r
if (!require("devtools")) install.packages("devtools") # to install
library(devtools) # to load
```

**Step 2:** Install **stackr**
```r
devtools::install_github("thierrygosselin/stackr") # to install without vignettes
devtools::install_github("thierrygosselin/stackr", build_vignettes = TRUE)  # to install WITH vignettes
library(stackr) # to load
```

**Step 3 (optional): Parallel computing** Install an OpenMP enabled [randomForestSRC](http://www.ccs.miami.edu/~hishwaran/rfsrc.html) package to do imputation in parallel. Follow the steps in this [vignette](https://github.com/thierrygosselin/stackr/blob/master/vignettes/vignette_imputations_parallel.Rmd). You don't need to do this when updating **stackr**.


**Problems during installation:** see this [vignette](https://github.com/thierrygosselin/stackr/blob/master/vignettes/vignette_installation_problems.Rmd)

**Dependencies**

  * **Imports:** adegenet, data.table, ggplot2, lazyeval, parallel, purrr, randomForestSRC, readr, stringi, stringr, tidyr, utils, plyr, dplyr

  * **Suggests:** devtools, knitr, rmarkdown

  * **Remotes:** github::thierrygosselin/assigner

A quick way to install/load required packages and start using my packages (copy/paste the whole block):
```r
if (!require("pacman")) install.packages("pacman")
library("pacman")
pacman::p_load(devtools, reshape2, ggplot2, stringr, stringi, plyr, dplyr, tidyr, readr, purrr, data.table, ape, adegenet, parallel, lazyeval, randomForestSRC)
if (!require("stackr")){
  install_github("thierrygosselin/stackr", build_vignettes = TRUE)
  library("stackr")
}
if (!require("assigner")) {
  install_github("thierrygosselin/assigner", build_vignettes = TRUE)
  # if assigner was re-installed, uncomment and run the next line to install gsi_sim:
  #install_gsi_sim(fromSource = TRUE) 
  library("assigner")
}
```

## New features
Version, new feature and bug history now lives in the [NEWS.md file] (https://github.com/thierrygosselin/stackr/blob/master/NEWS.md)

**v.0.2.9**
* bug fix in `tidy_genomic_data`
* bug fix between stackr -> devtools -> github -> travis, [this page helped] (http://itsalocke.com/using-travis-make-sure-use-github-pat/)


For previous news:
[NEWS.md file] (https://github.com/thierrygosselin/stackr/blob/master/NEWS.md)

## Roadmap of future developments:

* Very soon: Joint Allele Frequency Spectrum from a *batch_x.sumstats.tsv* or a *batch_x.haplotypes.tsv* files
* Re-Integration with [strataG] (https://github.com/EricArcher/strataG.devel/tree/master/strataG.devel)
* Improved documentation and vignette
* Workflow tutorial
* More linkage map tools
* Use Shiny and ggvis when subplots or facets are available
* CRAN
* Interaction with [STACKS] (http://creskolab.uoregon.edu/stacks/) database (Web-interface)
* Reference genome tools
* Maybe try some integration with other GBS approaches: AftrRAD, pyRAD, dDocent
* ...suggestions ?

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


## Vignettes and examples

From a browser:
* [installation problems](https://github.com/thierrygosselin/stackr/blob/master/vignettes/vignette_installation_problems.Rmd)
* [parallel computing during imputations](https://github.com/thierrygosselin/stackr/blob/master/vignettes/vignette_imputations_parallel.Rmd) 
* [vcf2dadi](https://github.com/thierrygosselin/stackr/blob/master/vignettes/vignette_vcf2dadi.Rmd)
* [haplo2genind](https://github.com/thierrygosselin/stackr/blob/master/vignettes/vignette_haplo2genind.Rmd)

Inside R:
```r
browse(Vignettes("stackr")) # To browse vignettes
vignette("vignette_vcf2dadi") # To open specific vignette
```

Vignettes are in development, check periodically for updates.


## Citation:
To get the citation, inside R:
```r
citation("stackr")
```
