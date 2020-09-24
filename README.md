
<!-- badges: start -->

[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://tidyverse.org/lifecycle/#experimental)
[![Travis-CI Build
Status](https://travis-ci.org/thierrygosselin/stackr.svg?branch=master)](https://travis-ci.org/thierrygosselin/stackr)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/thierrygosselin/stackr?branch=master&svg=true)](https://ci.appveyor.com/project/thierrygosselin/stackr)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/stackr)](http://cran.r-project.org/package=stackr)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![DOI](https://zenodo.org/badge/14548/thierrygosselin/stackr.svg)](https://zenodo.org/badge/latestdoi/14548/thierrygosselin/stackr)
[![packageversion](https://img.shields.io/badge/Package%20version-2.1.0-orange.svg)](commits/master)
[![Last-changedate](https://img.shields.io/badge/last%20change-2020--09--24-brightgreen.svg)](/commits/master)
<!-- badges: end -->

# stackr: an R package to run stacks software pipeline

This is the development page of the **stackr**.

**What’s the difference with running stacks directly in the terminal?**

Besides running stacks within R, not much, tiny differences here and
there that speed up my RADseq workflow:

  - The philosophy of working by project with pre-organized folders.
  - Some important steps are **parallelized**.
  - You have more than 1 sequencing chip/lane ? This workflow will save
    you lots of time.
  - **Technical replicates**, inside or across chip/lanes are managed
    uniquely.
  - Noise reduction.
  - Data normalization.
  - **nightmares because of a crashed computer/cluster/server?** stackr
    manage stacks unique integer (previously called SQL IDs) throughout
    the pipeline. It’s integrated from the start, making it a breeze to
    just re-start your pipeline after a crash\!
  - **mismatch testing:** *de novo* mismatch threshold series is
    integrated inside `run_ustacks` and stackr will produce tables and
    figures automatically.
  - **catalog**: for bigger sampling size project, breaking down the
    catalog into several separate *cstacks* steps makes the pipeline
    more rigorous if your computer/cluster/server crash.
  - **logs** generated by stacks are read and transferred in
    human-readable tables/tibbles. Detecting problems is easier.
  - summary of different stacks modules: available automatically inside
    stackr pipeline, but also available for users who didn’t use stackr
    to run stacks.
  - For me all this = increased reproducibly.

**Who’s it for?**

  - It’s currently developed with my own projects in mind.
  - To help collaborators to get the most out of stacks.

It’s not for R or stacks beginners. stacks related issues should be
highlighted on [stacks google
group](https://groups.google.com/forum/?fromgroups#!forum/stacks-users).

## Installation

To try out the dev version of **stackr**, copy/paste the code below:

``` r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("thierrygosselin/stackr")
library(stackr)
```

## Citation:

To get the citation, inside R:

``` r
citation("stackr")
```

Web site with additional info:
<http://thierrygosselin.github.io/stackr/>

  - [Computer setup and
    troubleshooting](https://thierrygosselin.github.io/radiator/articles/rad_genomics_computer_setup.html)
  - [Vignettes](https://thierrygosselin.github.io/radiator/articles/index.html)

## Life cycle

stackr is maturing, but in order to make the package better, changes are
inevitable. Argument names are very stable and follows stacks
development closely.

  - Philosophy, major changes and deprecated functions/arguments are
    documented in life cycle section of functions.
  - The latest changes are documented in [changelog, versions, new
    features and bug
    history](http://thierrygosselin.github.io/stackr/news/index.html)
  - [issues](https://github.com/thierrygosselin/stackr/issues/new/choose)
    and
    [contributions](https://github.com/thierrygosselin/stackr/issues/new/choose)

## Stacks modules and RADseq typical workflow

**stackr** package provides wrapper functions to run
[STACKS](http://catchenlab.life.illinois.edu/stacks/)
*process\_radtags*, *ustacks*, *cstacks*, *sstacks*, *rxstacks* and
*populations* inside R.

Below, a flow chart showing the corresponding stacks modules and stackr
corresponding functions. ![](vignettes/stackr_workflow.png)
