#' @name plink_ibm
#' @title Identity-by-missingness (IBM) analysis.
#' @description Diagnose your missing genotypes pattern with PLINK's 
#' identity-by-missingness (IBM) analysis.
#' @param vcf.file The VCF file created by STACKS.
#' @param whitelist (optional) A whitelist of markers to keep, 
#' for VCFTools 2 columns are required: 'CHROM' and 'POS'. 
#' Default: \code{whitelist = NULL}.
#' @param denovo (logical) Default: \code{denovo = TRUE}. For de novo analysis 
#' the vcf file created by stacks as a \code{un} in the \code{CHROM} column. 
#' We need to change that internally for the function to work in VCFTools and 
#' PLINK.
#' @param blacklist.id (optional) A blacklist with individual ID. NO 
#' column header. The blacklist is in the directory (e.g. "blacklist.id.txt").
#' Default: \code{blacklist.id = NULL}
#' @param strata File with column header 'INDIVIDUALS' and any other 
#' columns containing e.g. population id, sequencer id, lanes, etc. 
#' Any info that might impact missing data. 
#' Willl be used in \code{strata.select} argument.
#' @param strata.select Control the colour in ggplot2 plot. 
#' e.g. use \code{"POP_ID")} to see the strata column \code{"POP_ID"} 
#' or \code{"LANES"}, if you have a column in your strata file that starts 
#' with LANES (for sequencing lane...).
#' @return The function returns a list with the identity-by-missingness results 
#' ($ibm) and a MDS plot (MultiDimensional Scaling) of ibm results ($plot)
#' In the working directory:
#' the modified vcf file, if de novo was selected, plink tfiles and ibm results.
#' @examples
#' \dontrun{
#' With raw VCF without filtering the loci:
#' ibm <- plink_ibm(
#' vcf.file = "batch_1.vcf", 
#' denovo = TRUE, 
#' strata = "population.map.txt", 
#' strata.select = "POP_ID"
#' )
#' 
#' With filtered data:
#' ibm <- plink_ibm(
#' vcf.file = "batch_1.vcf", 
#' whitelist = "whitelist.vcf.txt", 
#' denovo = TRUE, 
#' strata = "population.map.txt", 
#' strata.select = "POP_ID"
#' )
#' 
#' If the object for the function is 'ibm' then:
#' to access the figure
#' ibm.plot <- ibm$plot
#' The ibm results of plink in a dataframe with your strata:
#' ibm.results <- ibm$ibm
#' }
#' @export
#' @rdname plink_ibm
#' @import reshape2
#' @import dplyr
#' @references Purcell S, Neale B, Todd-Brown K et al. (2007) 
#' PLINK: a tool set for whole-genome association and population-based linkage 
#' analyses. American Journal of Human Genetics, 81, 559–575.
#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011) 
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences. 
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013) 
#' Stacks: an analysis tool set for population genomics. 
#' Molecular Ecology, 22, 3124-3140.
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

plink_ibm <- function(vcf.file, whitelist = NULL, denovo = TRUE, blacklist.id = NULL, strata, strata.select){
  
  IID <- NULL
  FID <- NULL
  C1 <- NULL
  C2 <- NULL
  
  # Checking for missing and/or default arguments ******************************
  if (missing(vcf.file)) stop("VCF file required")
  if (missing(whitelist)) whitelist <- NULL # no Whitelist
  if (missing(denovo)) denovo <- TRUE
  if (missing(blacklist.id)) blacklist.id <- NULL
  if (missing(strata)) stop("strata file required")
  if (missing(strata.select)) stop("strata.select argument is required")
  
  
  
  # empty list for results
  res <- list()
  
  # de novo
  if(denovo == TRUE){
    # For de novo assembly, modify Stacks vcf file chrom column from "un" to "1"
    system(paste("perl -pe 's/^un/1/' ", vcf.file, " > ", "batch_1.modified.ibm.vcf", sep = ""))
    plink.input.file <- "batch_1.modified.ibm.vcf"
  } else {
    plink.input.file <- vcf.file
  }
  
  # no whitelist
  if (is.null(whitelist)) {
    message("No whitelist to filter the VCF")
    
    # no blakclist id
    if(is.null(blacklist.id)){
      message("No blacklisted id to apply to the VCF")
      
      # creates a PLINK tfile
      system(paste("vcftools --vcf ", plink.input.file, " --plink-tped", " --out plink.tfile.ibm"))
      
      # IBM analysis
      system(paste("plink --tfile ", "plink.tfile.ibm", "--cluster missing --mds-plot 4", "--out plink.ibm"))
      
      ibm.data <- "plink.ibm.mds"
      
    } else { # with blacklisted id
      message("Blacklisted id used to filter the VCF")
      
      # remove individuals
      system(paste("vcftools --vcf ", plink.input.file, " --remove ", blacklist.id, " --recode", " --out batch_1.modified.ibm.id.removed"))
      plink.input.file <- "batch_1.modified.ibm.id.removed.recode.vcf"
      
      # creates a PLINK tfile
      system(paste("vcftools --vcf ", plink.input.file, " --plink-tped", " --out plink.tfile.ibm.id.removed"))
      
      # IBM analysis
      system(paste("plink --tfile ", "plink.tfile.ibm.id.removed", "--cluster missing --mds-plot 4", "--out plink.ibm.id.removed"))
      
      ibm.data <- "plink.ibm.id.removed.mds"
    }
    
  } else { #with whitelist
    message("Filtering the VCF with the whitelist from your directory")
    if(is.null(blacklist.id)){ # no blacklist id
      message("No blacklisted id to apply to the VCF")
      
      system(paste("vcftools --vcf ", plink.input.file, " --positions ", whitelist, " --recode", "--out batch_1.modified.ibm.filtered"))
      plink.input.file <- "batch_1.modified.ibm.filtered.recode.vcf"
      
      # creates a PLINK tfile
      system(paste("vcftools --vcf ", plink.input.file, " --plink-tped", " --out plink.tfile.ibm.filtered"))
      
      # IBM analysis
      system(paste("plink --tfile ", "plink.tfile.ibm.filtered", "--cluster-missing --mds-plot 4", "--out plink.ibm.filtered"))
      ibm.data <- "plink.ibm.filtered.mds"
      
    } else { # with blacklisted id
      message("Blacklisted id used to filter the VCF")
      
      #whitelist
      system(paste("vcftools --vcf ", plink.input.file, " --positions ", whitelist, " --recode", "--out batch_1.modified.ibm.filtered"))
      plink.input.file <- "batch_1.modified.ibm.filtered.recode.vcf"
      
      # remove individuals
      system(paste("vcftools --vcf ", plink.input.file, " --remove ", blacklist.id, " --recode", " --out batch_1.modified.ibm.filtered.recode.id.removed"))
      plink.input.file <- "batch_1.modified.ibm.filtered.recode.id.removed.recode.vcf"
      
      # creates a PLINK tfile
      system(paste("vcftools --vcf ", plink.input.file, " --plink-tped", " --out plink.tfile.ibm.filtered.id.removed"))
      
      # IBM analysis
      system(paste("plink --tfile ", "plink.tfile.ibm.filtered.id.removed", "--cluster-missing --mds-plot 4", "--out plink.ibm.filtered.id.removed"))
      ibm.data <- "plink.ibm.filtered.id.removed.mds"
    }
  }
  
  # strata
  ibm <- read_table(ibm.data, col_names = T, col_types = "ccidddd") %>%
    select(-c(IID), INDIVIDUALS=FID) %>%
    left_join(
      read_tsv(strata, 
               col_names = T), 
      by = "INDIVIDUALS")
  
  # Include data frame in results output of the function
  res$ibm <- ibm
  
  # IBM Figures
  ibm.plot <- ggplot(ibm, aes(x = C1, y = C2), environment = environment())+
    geom_point(aes_string(colour = strata.select))+
    labs(title = "MultiDimensional Scaling (MDS)\n of Identity by Missing (IBM)")+
    labs(x = "PC1")+
    labs(y = "PC2")+
    theme(axis.title.x = element_text(size = 12, family = "Helvetica", face = "bold"),
          axis.title.y = element_text(size = 12, family = "Helvetica", face = "bold"),
          legend.title = element_text(size = 12, family = "Helvetica", face = "bold"), 
          legend.text = element_text(size = 12, family = "Helvetica", face = "bold"), 
          strip.text.x = element_text(size = 12, family = "Helvetica", face = "bold")
    )
  res$plot <- ibm.plot
  return(res)
}
