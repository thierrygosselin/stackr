# stackr

[![Travis-CI Build Status](https://travis-ci.org/thierrygosselin/stackr.svg?branch=master)](https://travis-ci.org/thierrygosselin/stackr)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/stackr)](http://cran.r-project.org/package=stackr)
[![DOI](https://zenodo.org/badge/14548/thierrygosselin/stackr.svg)](https://zenodo.org/badge/latestdoi/14548/thierrygosselin/stackr)

The goal of **stackr** is to make GBS/RAD data produced by [STACKS] (http://catchenlab.life.illinois.edu/stacks/) easy to analyse in R.

This is the development page of the **stackr** package for the R software, optimized for *de novo* and population genetics.

**Use stackr to:**
* Read and modify *batch_x.sumstats.tsv* and *batch_x.haplotypes.tsv* files
* Transform the VCF file, *batch_x.vcf*, into a tidy format to visualise and filter summary statistics within R
* Filters genetic markers based on: coverage (read depth, REF and ALT allele depth), genotype likelihood, number of individuals, number of populations, minor allele frequency (local and global), observed heterozygosity and inbreeding coefficient (Fis)
* `ggplot2`-based plotting to view distributions of summary statistics and create publication-ready figures
* Convert data into *genepop*, *genind*, *fstat*, *gtypes*, *betadiv* and *dadi* files or objects for easy integration with other software or R packages like [adegenet] (https://github.com/thibautjombart/adegenet), [strataG] (https://github.com/EricArcher/strataG.devel/tree/master/strataG.devel), [hierfstat] (https://github.com/jgx65/hierfstat), [pegas] (https://github.com/emmanuelparadis/pegas) and [poppr] (https://github.com/grunwaldlab/poppr)
* Impute GBS markers without a genetic map using Random Forest within the *haplo2genepop*, *haplo2genind*, *haplo2hierfstat*, *haplo2gtypes*, *haplo2colony*, *vcf2genind*, *vcf2hierfstat*, *vcf2betadiv*, *vcf2dadi* and *vcf_imputation* functions 

**Requirement:**
Because STACKS is always under development (more than 100 versions so far!), 
**stackr** will work best with Stacks version >= 1.29. 
Send me an e-mail if you desperately need to use prior versions.

## Installation
You can try out the dev version of **stackr**. Follow the 3 steps below:

Step 1 You will need the package *devtools*
```r
if (!require("devtools")) install.packages("devtools") # to install
library(devtools) # to load
```

Step 2 Install **stackr**:
```r
install_github("thierrygosselin/stackr") # to install
library(stackr) # to load
```

Step 3 For faster imputations, you need to install an OpenMP enabled **randomForestSRC package** [website](http://www.ccs.miami.edu/~hishwaran/rfsrc.html).

Option 1: From source (Linux & Mac OSX)

```r
# Terminal
cd ~/Downloads
curl -O https://cran.r-project.org/src/contrib/randomForestSRC_2.0.7.tar.gz
tar -zxvf randomForestSRC_2.0.7.tar.gz
cd randomForestSRC
autoconf
# Back in R:
install.packages(pkgs = "~/Downloads/randomForestSRC", repos = NULL, type = "source")
```
Option 2: Use a pre-compiled binary (Mac OSX & Windows) [instructions here] (http://www.ccs.miami.edu/~hishwaran/rfsrc.html) or quick copy/paste solution below:

```r
# Mac OSX
library("devtools")
install_url(url = "http://www.ccs.miami.edu/~hishwaran/rfsrc/randomForestSRC_2.0.7.tgz")
```

```r
# Windows
library("devtools")
install_url(url = "http://www.ccs.miami.edu/~hishwaran/rfsrc/randomForestSRC_2.0.7.zip")
```

**Problems during installation:**

Sometimes you'll get warnings while installing dependencies required for **stackr** or other R packages.
```r
Warning: cannot remove prior installation of package ‘stringi’
```

To solve this problem: 

Option 1. Delete the problematic packages manually and reinstall. On MAC computers, in the **Finder**, use the shortcut **cmd+shift+g**, or in the menu bar : **GO -> Go to Folder**, copy/paste the text below:
```r
/Library/Frameworks/R.framework/Resources/library
#Delete the problematic packages.
```

Option 2. If you know your way around the terminal and understand the consequences of using **sudo rm -R** command, here something faster to remove problematic packages:
```r
sudo rm -R /Library/Frameworks/R.framework/Resources/library/package_name
# Changing 'package_name' to the problematic package.
# Reinstall the package.
```
## Parallel computing in R and stackr

On Mac OSX using OpenMP greatly reduce the computation time for the imputations. Follow the instructions [here] (http://gbs-cloud-tutorial.readthedocs.org/en/latest/03_computer_setup.html#update-your-computer-s-compiler) to update your computer's compiler (5 min step). 

You need to tell R which compilers to use. Use TextWrangler or follow the lines below:
```r
cd ~
nano .R/Makevars
```

Enter the text below:
```r
CC=/usr/local/bin/gcc
CXX=/usr/local/bin/g++
FC=/usr/local/bin/gfortran
F77=/usr/local/bin/gfortran
PKG_LIBS = -fopenmp -lgomp
PKG_CFLAGS= -O3 -Wall -pipe -pedantic -std=gnu99 -fopenmp
CFLAGS= -O3 -Wall -pipe -pedantic -std=gnu99 -fopenmp
SHLIB_OPENMP_CFLAGS = -fopenmp
SHLIB_OPENMP_CXXFLAGS = -fopenmp
SHLIB_OPENMP_FCFLAGS = -fopenmp
SHLIB_OPENMP_FFLAGS = -fopenmp
```
Save and Exit with: crt-o, enter, crt-x. Preferably, re-install all packages depending on OpenMP or Rcpp:

```r
install.packages("Rcpp", type = "source")
install.packages("dplyr", type = "source")
```
## New features
Version, new feature and bug history now lives in the [NEWS.md file] (https://github.com/thierrygosselin/stackr/blob/master/NEWS.md)

**v.0.2.7**
* Added a `NEWS.md` file to track changes to the package.
* New function: `individuals2strata`. Several functions in **stackr** and 
[assigner] (https://github.com/thierrygosselin/assigner) requires a `strata`
argument, i.e. a data frame with the individuals and associated groupings. 
You can do it manually, however, if your individuals have a consistent naming scheme 
(e.g. SPECIES-POPULATION-MATURITY-YEAR-ID = CHI-QUE-ADU-2014-020), 
use this function to rapidly create a strata file.
* New function: `tidy_genomic_data`. 
Transform common genomic dataset format in a tidy data frame. Used internally in
**stackr** and [assigner] (https://github.com/thierrygosselin/assigner)
and might be of interest for users.
* New function: `read_long_tidy_wide`. Read genomic data frames in long/tidy and wide format.
Used internally in **stackr** and [assigner] (https://github.com/thierrygosselin/assigner)
and might be of interest for users.
* New function: `stackr_imputations_module`. 
Map-independent imputation of missing genotype using Random Forest
or the most frequent category. Impute genotypes or alleles. 
Used internally in **stackr** and [assigner] (https://github.com/thierrygosselin/assigner)
and might be of interest for users.
* New function: `find_duplicate_id`
Compute pairwise genome similarity to highligh potential duplicate individuals.

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
Look at step 1 as a quality insurance step. We need to modify the data to play with it efficiently in R. To have reliable summary statistics, you first need good coverage of your alleles to call your genotypes, good genotype likelihood, enough individuals in each sampling sites and enough putative populations with your markers... Step 2 is where the actual work is done to remove artifactual and uninformative markers based on summary statistics of your markers.
![](vignettes/stackr_workflow.png)

## Example 

**Using *haplo2genind* function to do a DAPC analysis of your data (5 steps).**

Step 1. Load the necessary librairies, here is an example of how to do this:
```r
library(adegenet)
library(stackr)
```


*Dependencies*: here the list of packages that **stackr** is depending on.
```r
if (!require("reshape2")) install.packages("reshape2")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("stringr")) install.packages("stringr")
if (!require("stringi")) install.packages("stringi")
if (!require("plyr")) install.packages("plyr")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("readr")) install.packages("readr")
if (!require("purrr")) install.packages("purrr")
if (!require("data.table")) install.packages("data.table")
if (!require("lazyeval")) install.packages("lazyeval")
if (!require("adegenet")) install.packages("adegenet")
if (!require("parallel")) install.packages("parallel")
if (!require("stringdist")) install.packages("stringdist")
if (!require("foreach")) install.packages("foreach")
if (!require("doParallel")) install.packages("doParallel")
```
If you don't have them, no worries, it's intalled automatically during **stackr** installation. If you have them, it's your job to update them, because i'm using the latest versions...

*When working in R I usually do this:*
```r
# Clean my desk
rm(list=ls())

# load the required libraries
library(reshape2)
library(ggplot2)
library(stringr)
library(stringi)
library(plyr)
library(dplyr) # load this package after plyr to work properly
library(tidyr)
library(readr)
library(randomForestSRC)
library(doParallel)
library(stackr)
library(purrr)
```


Step 2. Set your working directory (e.g. the path to your **stacks** output files and 
where you want the output to be saved):

```r
setwd("/Users/thierry/Dropbox/brook_charr_pop/01_stacks_populations")
```
Step 3. Missing genotypes: 


First remove individuals with more than 30% missing genotypes 
from the *batch_1.haplotypes.tsv* file. Explore this parameter with different values. 


You can also provide the function with a whitelist of loci to keep (after filtering).
We are interested in the the blacklisted id output ("blacklisted.id.30.txt"),
but the function also outputs many things, see the function documentation.
```r
blacklisted.id <- missing_genotypes(haplotypes.file = "batch_1.haplotypes.tsv", 
whitelist.loci = "new.whitelist.txt", pop.id.start = 5, pop.id.end = 7, 
missing.geno.threshold = 30)
```

Step 4. Use the *haplo2genind* function to convert the haplotype file created by 
**stacks** into a genind object ready to use in **adegenet**. 

I use the whitelist of loci created after filtering the data and filter out the individuals with more than 30% of missing genotypes (with the blacklisted individuals, created above). 

I also ask for imputation of the data using Random Forest.

```r
genind.sturgeon <- haplo2genind(haplotypes.file = "batch_1.haplotypes.tsv", whitelist.loci = "my.whitelist.txt", blacklist.id = "blacklisted.id.30.txt", pop.levels = c("LSL", "DRM", "JEN", "LAN", "GRA", "BUR", "GUL", "LLI", "ANG", "WEI", "FOX", "HAY", "GOD", "CHU"), pop.id.start = 5, pop.id.end = 7, imputations = "rf", imputations.group = "populations", num.tree = 100, split.number = 100, iteration.rf = 10, verbose = FALSE)
```

You can see that the object created is not yet a genind object because it contains 2 things: the imputed data and the data without imputation. To access both genind dataset:
```r
names(genind.sturgeon)
genind.sturgeon.noimputation <- genind.sturgeon$no.imputation
genind.sturgeon.imputed <- genind.sturgeon$imputed
```

Step 5. These 2 genind objects can be use directly in **adegenet**:
```r
dapc.optim.a.score <- optim.a.score(dapc(genind.sturgeon.imputed, n.da = 100, n.pca = 50))
dapc.optim.a.score$best
```
## All-in-one filter ![](vignettes/stackr_all-in-one_filters.png)
## Functions found in **stackr** ![](vignettes/stackr_functions.png)

Vignettes are in development, check periodically for updates.

## Citation:
To get the citation, inside R:
```r
citation("stackr")
```
