# stackr v.0.5.4 2017-04-06

* `detect_duplicate_genomes` : huge speed bump for pairwise genome similarity method.
Instead of hours the range is more in minutes.


# stackr v.0.5.3 2017-04-05

* `stackr_imputations_module` : better integration of VCF with haplotypes so
that nucleotide information is kept during imputations. 
* `filter_fis` : bug fix when no heterozygote were found. Thanks to Manuel Lamothe.


# stackr v.0.5.2 2017-03-27

* `stackr_imputations_module` : work on faster on-the-fly random forest and
extreme gradient tree boosting algorithm.


# stackr v.0.5.1 2017-03-21

Major work on `tidy_genomic_data`:
* `platypus` vcf files are correctly imported
* more efficient when working with vcf files
* better parallelization during parsing and cleaning

# stackr v.0.5.0

* better parallelization of `summary_haplotypes` function. With progress bar...
* bug fix with `summary_haplotypes` not properly summarizing info when no assembly artifacts were found



# stackr v.0.4.9

* bug fix with detection of mixed bi/multi allelic dataset. The bug was detected
in `tidy_genomic_data` and `genomic_converter` functions.


# stackr v.0.4.8

* safer use and better parsing of `strataG` object to work with tidy data and pass Travis CI.


# stackr v.0.4.7

* Better parsing of genepop file with 2 characters for allele coding
* Only 10% of markers are now used for increase speed of bi-allelic markers detection
* Work on imputation module that will be functional in the next version bump
* Better code for the progress bar (Linux and Mac only) that now shows an ETA along the progress.




# stackr v.0.4.6

* I'm pleased to announce that `stackr` parallel mode now works with **Windows**!
Nothing to install, just need to choose the number of CPU,
the rest is done automatically.
* `haplo2colony` is deprecated. Use the new function called `write_colony`!
* `write_colony`: works similarly to the deprecated function `haplo2colony`,
      * with the major advantage that it's no longer restricted to STACKS 
      haplotypes file. 
      * The function is using the `tidy_genomic_data` module
      to import files. So you can choose one of the 10 input file formats supported
      by `stackr`!
      * other benefits also include the possibility to efficiently test MAF, 
      snp.ld, haplotypes/snp approach, whitelist of markes, 
      blacklist of individuals, blacklist of genotypes, etc. with the buit-it
      arguments.
      * the function only **keeps markers in common** between populations/groups
      and **is removing monomorphic markers**.
      * **Note:** there are several *defaults* in the function and
      it's a complicated file format, so make sure to read the function
      documentation, please, and `COLONY` manual.


# stackr v.0.4.5

* temporary fix to `tidy_genomic_data` to read unconventional Tassel VCF
* new function `ibdg_fh` computes the FH measure that was previously
computed in `summary_haplotypes`.
It now works with biallelic and multiallelic data.
    The FH measure is based on the excess in the observed number of homozygous
genotypes within an individual relative to the mean number of homozygous
genotypes expected under random mating (see function for details). The `IBDg` in
the name is because the measure is a proxy of the realized proportion of the genome
that is identical by descent by reference to the current population 
under hypothetical random mating.
* `missing_visualization` now computes the FH measure and look for correlation
with average missingness per individual.
* `tidy_stacks_haplotypes_vcf` is now deprecated in favor of using `tidy_genomic_data`
that will import haplotypic vcf files.

# stackr v.0.4.4

* several updates to make functions faster.
* `stackr_imputations_module` no longer imputes globally after imputations
by populations. Instead, use `common.markers` or not to test impacts.
*  bug fix with `ref_alt_alleles` that was not working properly inside
the imputation module.
* `snp_ld` is not a separate module available for users. Check documentation.
* `missing_visualization` now show the proportion of variance with plot axis text.


# stackr v.0.4.3
* bug fix in `summary_haplotypes` stemming from a new `readr` version
* `artifacts` replace `paralogs` in `summary_haplotypes`

# stackr v.0.4.2
* `gtypes` object from [strataG] (https://github.com/EricArcher/strataG) package
can now be read/write in/out of **Stackr** using the `tidy_genomic_data` and 
`genomic_converter` functions.

# stackr v.0.4.1
* update `missing_visualization` function to include more PCoA plots

# stackr v.0.4.0
* couple of bug fix for detecting file formats

# stackr v.0.3.9
* several performance update
* couple of bug fix for detecting file formats

# stackr v.0.3.8
* fixed a bug in `filter_genotype_likelihood`, since the updated function to the 
interactive mode, some old code where still present in if/else sentences, breaking 
the code. Thanks to Jaromir Guzinski for the bug report.

# stackr v.0.3.7
* fixed a bug in `write_vcf`, the function was using REF/ALT coding in integer 
not character format. This function is used inside `vcf_imputation` and 
sometimes inside `genomic_converter`. Thanks to @jeansebastienmoore for 
highlighting the problem.

# stackr v.0.3.6
* fixed a bug in `vcf_imputation`, the function now calls `genomic_converter` 
with all the bells and whistles of that function (updated vcf import and imputations modules)


# stackr v.0.3.5
* updated tidy_genepop to read other flavors of the famous file format
* extracted a code block to create a new function called `tidy_fstat`

# stackr v.0.3.4
* updated documentation
* bug fix in `summary_haplotypes` introduced by the new version of `dplyr::distinct` (v.0.5.0)
* calculations of Pi is done in parallel inside `summary_haplotypes`

# stackr v.0.3.3
* `tidy_genomic_data`: added a check to throw an error when pop.levels != the pop.id in strata

# stackr v.0.3.2
* `genomic_converter` including all the `vcf2...` function can now use phase/unphase genotypes.
Some **pyRAD** vcf (e.g. v.3.0.64) have a mix of GT format with `/` and `|`. 
e.g. missing GT  = `./.` and genotyped individuals = `0|0`. 
I'm not sure it follows [VCF specification](http://samtools.github.io/hts-specs/VCFv4.2.pdf), 
but **stackr** can now read those vcf files.
* `vcf2dadi` is more user-friendly for scientist with in- and out-group metadata, using STACKS or not.


# stackr v.0.3.1
* Bug fix: combined use of `if (getRversion() >= "2.15.1") utils::globalVariables("variable")` 
and `@inheritParams` was not showing all the argument description.

# stackr v.0.3.0
* Update that makes my coding life easier.
* Several internal functions to convert from a tidy dataframe to:
`vcf`, `plink`, `genind`, `genlight`, `gtypes`, `hierfstat`, `genepop`, `structure` and `betadiv` 
are now separate modules available to users (look for `write_...` with the outputformat)
* New function `genomic_converter`: If you want the to convert from the supported 
input file formats to many output formats, at once, this is the function. 
With the new function `genomic_converter`, import and imputations are only 
done once, saving time if you were generating different output WITH imputations.
* Change: all the `vcf2...` functions (excep `vcf2dadi`) are now a shorcut of `genomic_converter`. 
 This is particularly interesting and faster if you were generating different
 output WITH imputations. This makes the functions `vcf2...` and `genomic_converter`
 easier to debug for me and more stable for users.
* Deprecated: the `haplo2...` functions are all deprecated and replaced by 
`genomic_converter`, **except haplo2colony** that requires so many arguments that it 
would be too complicated, for now, to integrate with `genomic_converter`.
* New feature: when arguments `pop.select`, `blacklist.id` and `imputation.method`
are used, the REF and ALT alleles are now re-computed to account for the filters 
and imputations.


# stackr v.0.2.9
* bug fix in `tidy_genomic_data`
* bug fix between stackr -> devtools -> github -> travis, [this page helped] (http://itsalocke.com/using-travis-make-sure-use-github-pat/)


# stackr v.0.2.8
* bug fix in `tidy_genomic_data` while using data.table::melt.data.table instead 
of tidyr::gather, and forgot to 
(i) add variable.factor = FALSE when melting the vcf and (ii) use as_data_frame
at the end of the melting to be able to continue working with dplyr verbs.


# stackr v.0.2.7
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


# stackr v.0.2.6
* dart2df_genind_plink: swiss army knife tool to prepare DArT output file (wide 
or binary format) for population genetics analysis. Import, filter and transform 
a DArT output file to different format: tidy data frame of genotypes, genind object 
and/or PLINK `tped/tfam` format. Map-independent imputation also available.


# stackr v.0.2.5
* vcf2plink: to easily convert a VCF file created in STACKS to a PLINK input 
file (tped/tfam format). This function comes with the commonly used arguments 
in **stackr**: map-independent imputation, whitelist, blacklist, common marker filtering, etc.

* data_pruning: to prune your dataset with whitelist, blacklist of individuals, 
erase genotypes, use common markers and other filtering (see function argument 
while waiting for the upcomming documentation).

# stackr v.0.2.4
* updated the vcf_imputation function for the commonly used arguments in **stackr**.

# stackr v.0.2.3
* vcf2dadi: to easily convert a VCF file created in STACKS to a dadi input file.
This function comes with the commonly used arguments in **stackr**: 
map-independent imputation, whitelist, blacklist, common marker filtering, etc.

# stackr v.0.2.2
* vcf2genepop: to easily convert a VCF file created in STACKS to a genepop input file.
This function comes with the commonly used arguments in **stackr**: 
map-dependent imputation, whitelist, blacklist, etc. For the haplotype version, see
haplo2genepop.

# stackr v.0.2.1
* 'read_stacks_vcf' can now use a whitelist or blacklist of loci that works with CHROM and/or SNP and/or LOCUS.
* 'filter_maf', 'filter_fis', 'filter_het' and 'filter_genotype_likelihood' now works by haplotypes or SNP.

# stackr v.0.2.0
Introducing several new functions: 
* vcf2betadiv: to easily convert a VCF file created in STACKS to a betadiv input file.
* vcf2genind: same as haplo2genind but works with SNP instead of haplotypes.
* vcf2hierfstat: same as haplo2hierfstat but works with SNP instead of haplotypes.

# stackr v.0.1.5
Introducing *haplo2gsi_sim* function.
* Conversion of STACKS haplotypes file into a gsi_sim data input file.
* Markers can be subsampled.
* Map-independent imputations using Random Forest or the most frequent allele are options also available for this function.
* [gsi_sim] (https://github.com/eriqande/gsi_sim) is a tool developed by Eric C. Anderson for doing and simulating genetic stock identification.

# stackr v.0.1.4
Introducing *haplo2fstat* function.
Conversion of STACKS haplotypes file into a hierfstat object and fstat file.
Access all the functions in the R package [hierfstat] (https://github.com/jgx65/hierfstat).

# stackr v.0.1.3
Map-independent imputations of a VCF file created by STACKS. 
Two options are available for imputations: using Random Forest or the most frequent allele.

Before imputations, the VCF file can be filtered with:

* a whitelist of loci (to keep only specific loci...)
* a blacklist of individuals (to remove individuals or entire populations...)
* also, a list of genotypes with bad coverage and/or genotype likelihood can be supplied to erase the genotypes before imputations (for more details look at the function: blacklist_erase_genotype).

# stackr v.0.1.2
**The *summary_haplotypes* function now outputs:**
* Putative paralogs, consensus, monomorphic and polymorphic loci
* The haplotype statistics for the observed and expected homozygosity and 
heterozygosity
* Wright’s inbreeding coefficient (Fis)
* Proxy measure of the realized proportion of the genome that is identical
by descent (IBDG). The FH measure is based on the excess in the observed number
of homozygous genotypes within an individual relative to the mean number of 
homozygous genotypes expected under random mating (Keller et al., 2011; 
Kardos et al., 2015).
* Nucleotide diversity (Pi), considering the consensus loci in the catalog 
(i.e. reads with no variation between population). It's Nei & Li (1979) 
function, adapted to the GBS reality.

Keller MC, Visscher PM, Goddard ME. 2011. Quantification of inbreeding due to 
distant ancestors and its detection using dense single nucleotide polymorphism
data. Genetics, 189, 237–249.

Kardos M, Luikart G, Allendorf FW. 2015. Measuring individual inbreeding in the 
age of genomics: marker-based measures are better than pedigrees. 
Heredity, 115, 63–72.

Nei M, Li WH. 1979. Mathematical model for studying genetic variation in terms
of restriction endonucleases. Proceedings of the National Academy of Sciences 
of the United States of America, 76, 5269–5273.

**The *haplo2colony* function**
* Converts the file to the required *COLONY* input file
* Can filter the haplotypes file with a whitelist of loci 
and a blacklist of individuals
* Can impute the data with Random Forest or the most frequent category
* Use the *print.all.colony.opt* to output all COLONY options to the file.
This however requires manual curation of the file to work directly with COLONY. 
