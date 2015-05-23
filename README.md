# stackr
The goal of *stackr* is to make GBS/RAD data produced by STACKS easy to analyse in R.

This is the development page of the *stackr* package for the R software.

* Read, modify *batch_x.sumstats.tsv and batch_x.haplotypes.tsv* files.
* Create a tidy format of VCF file, *batch_x.vcf*, to visualise and filter summary statistics.
* Filters genetic markers based on: coverage (read depth, REF and ALT allele depth), genotype likelihood, the number of individuals, the number of populations, observed heterozygosity and inbreeding coefficient (Fis).
* View distribution of summary sttistics and create publication-ready figures
* Convert data into genind object for easy integration with adegenet, hierfstat, pegas.




### Installation*
You can try out the dev version of *stackr*. You will need the package *devtools*.

```r
install.packages("devtools")
install_github("thierrygosselin/stackr")
library("stackr")
```
