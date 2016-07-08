# Write a adegenet genind object from STACKS VCF file

#' @name vcf2genind
#' @title VCF to \code{adegenet} \code{\link[adegenet]{genind}} object with filters and data imputation

#' @description For full details of the function, please use 
#' \pkg{stackr} \code{\link[stackr]{genomic_converter}}. This function is a shorcut
#' to output only genind object.
#' @inheritParams genomic_converter 
#' @inheritParams tidy_genomic_data 
#' @inheritParams write_genepop
#' @inheritParams write_genind 
#' @inheritParams write_genlight 
#' @inheritParams write_structure
#' @inheritParams write_plink
#' @inheritParams write_vcf
#' @inheritParams write_gtypes
#' @inheritParams write_hierfstat
#' @inheritParams stackr_imputations_module 

#' @export
#' @rdname vcf2genind
#' @import dplyr
#' @import stringi
#' @importFrom data.table fread
#' @import parallel
#' @import dplyr
#' @import stringi

#' @references Jombart T (2008) adegenet: a R package for the multivariate
#' analysis of genetic markers. Bioinformatics, 24, 1403-1405.
#' @references Jombart T, Ahmed I (2011) adegenet 1.3-1: 
#' new tools for the analysis of genome-wide SNP data. 
#' Bioinformatics, 27, 3070-3071.

#' @seealso \code{adegenet} is available on CRAN \url{http://cran.r-project.org/web/packages/adegenet/} and github \url{https://github.com/thibautjombart/}
#' \code{\link[stackr]{genomic_converter}}

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

vcf2genind <- function(
  data,
  output,
  filename = NULL,
  blacklist.id = NULL,
  blacklist.genotype = NULL,
  whitelist.markers = NULL,
  monomorphic.out = TRUE,
  snp.ld = NULL,
  common.markers = TRUE,
  maf.thresholds = NULL,
  maf.pop.num.threshold = 1,
  maf.approach = "SNP",
  maf.operator = "OR",
  max.marker = NULL,
  strata = NULL,
  pop.levels = NULL,
  pop.labels = NULL,
  pop.select = NULL,
  imputation.method = NULL,
  impute = "genotype",
  imputations.group = "populations",
  num.tree = 100,
  iteration.rf = 10,
  split.number = 100,
  verbose = FALSE,
  parallel.core = detectCores()-1
) {
  
  res <- genomic_converter(
    data,
    output = "genind",
    filename = filename,
    blacklist.id = blacklist.id,
    blacklist.genotype = blacklist.genotype,
    whitelist.markers = whitelist.markers,
    monomorphic.out = monomorphic.out,
    snp.ld = snp.ld,
    common.markers = common.markers,
    maf.thresholds = maf.thresholds,
    maf.pop.num.threshold = maf.pop.num.threshold,
    maf.approach = maf.approach,
    maf.operator = maf.operator,
    max.marker = max.marker,
    strata = strata,
    pop.levels = pop.levels,
    pop.labels = pop.labels,
    pop.select = pop.select,
    imputation.method = imputation.method,
    impute = impute,
    imputations.group = imputations.group,
    num.tree = num.tree,
    iteration.rf = iteration.rf,
    split.number = split.number,
    verbose = verbose,
    parallel.core = parallel.core
  )
  return(res)
}

