# VCF to plink with data imputation using Random Forest

#' @name vcf2plink
#' @title VCF to plink with filters and data imputation

#' @description For full details of the function, please use 
#' \pkg{stackr} \code{\link[stackr]{genomic_converter}}. This function is a shorcut
#' to output only plink tped/tmap files.
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
#' @rdname vcf2plink
#' @import dplyr
#' @import stringi
#' @importFrom data.table fread

#' @references Danecek P, Auton A, Abecasis G et al. (2011)
#' The variant call format and VCFtools.
#' Bioinformatics, 27, 2156-2158.
#' @references Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, 
#' Bender D, et al. 
#' PLINK: a tool set for whole-genome association and population-based linkage 
#' analyses. 
#' American Journal of Human Genetics. 2007; 81: 559–575. doi:10.1086/519795


#' @seealso \code{\link[stackr]{genomic_converter}}


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


vcf2plink <- function(
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
    output = "plink",
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
