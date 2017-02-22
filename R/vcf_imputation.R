# VCF data imputation using Random Forest


#' @name vcf_imputation
#' @title VCF data imputation
#' @description For full details of the function, please use 
#' \pkg{stackr} \code{\link[stackr]{genomic_converter}}. This function is a shorcut
#' to output only genepop file.
#' @param data A VCF file.
#' @inheritParams genomic_converter 
#' @inheritParams tidy_genomic_data 
#' @inheritParams write_vcf
#' @inheritParams stackr_imputations_module 
#' @seealso \code{\link[stackr]{genomic_converter}}

#' @export
#' @rdname vcf_imputation
#' @import dplyr
#' @import stringi
#' @importFrom data.table fread
#' @importFrom purrr flatten_chr
#' @import parallel


#' @examples
#' \dontrun{
#' vcf_imputation(
#' data = "batch_1.vcf", 
#' whitelist.markers = "whitelist.loci.txt",
#' strata = "my.strata.tsv",
#' pop.levels = c("PAN", "COS"),
#' common.markers = TRUE, 
#' imputation.method = "max", 
#' impute = "allele", 
#' hierarchical.levels = "populations", 
#' num.tree = 100, 
#' iteration.rf = 10, 
#' split.number = 100, 
#' verbose = FALSE, 
#' parallel.core = 8
#' )
#' }



#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011) 
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences. 
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013) 
#' Stacks: an analysis tool set for population genomics. 
#' Molecular Ecology, 22, 3124-3140.
#' @references Ishwaran H. and Kogalur U.B. (2015). Random Forests for Survival,
#'  Regression and Classification (RF-SRC), R package version 1.6.1.
#' @references Ishwaran H. and Kogalur U.B. (2007). Random survival forests
#' for R. R News 7(2), 25-31.
#' @references Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S. (2008).
#' Random survival forests. Ann. Appl. Statist. 2(3), 841--860.
#' @seealso \code{randomForestSRC} is available on 
#' CRAN \url{http://cran.r-project.org/web/packages/randomForestSRC/} 
#' and github \url{https://github.com/ehrlinger/randomForestSRC}
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

vcf_imputation <- function(
  data,
  output = "vcf",
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
  hierarchical.levels = "populations",
  num.tree = 50,
  pred.mean.matching = 0,
  random.seed = NULL,
  verbose = FALSE,
  parallel.core = detectCores() - 1
) {
  res <- genomic_converter(
    data,
    output = "vcf",
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
    hierarchical.levels = hierarchical.levels,
    num.tree = num.tree,
    pred.mean.matching = pred.mean.matching,
    random.seed = random.seed,
    verbose = verbose,
    parallel.core = parallel.core
  )
  return(res)
}
