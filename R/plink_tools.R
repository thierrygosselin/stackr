#' @name make_tfam
#' @title Write a PLINK tfam file from a population map.
#' @description Make a tfam file for PLINK. Useful in ADMIXTURE for bootstrapping.
#' @param population.map The population map or individuals listed in one column. No headers.
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
#' @param pop.levels (optional) Character string with your populations ordered.
#' @param filename The name of the file written in the directory.
#' @export
#' @rdname make_tfam
#' @examples
#' \dontrun{
#' make_tfam(
#' population.map = "population_map_turle.txt", 
#' pop.id.start = 1, 
#' pop.id.end = 3, 
#' pop.levels = c("GAP", "ISA", "CRU"),
#' filename="turtle.tfam")
#' }
#' @references Purcell S, Neale B, Todd-Brown K et al. (2007) 
#' PLINK: a tool set for whole-genome association and population-based linkage 
#' analyses. American Journal of Human Genetics, 81, 559–575.


make_tfam <- function(population.map, pop.id.start, pop.id.end, pop.levels, filename) {
  
  if(missing(population.map)) stop("population.map argument required")
  if(missing(pop.id.start)) stop("pop.id.start argument required")
  if(missing(pop.id.end)) stop("pop.id.start argument required")
  if(missing(pop.levels)) pop.levels <- NULL
  if(missing(filename)) stop("filename argument required")
  
  
  X1 <- NULL
  INDIVIDUALS <- NULL
  POP_ID <- NULL
  
  if (is.vector(population.map) == "TRUE") {
    population.map <- read_tsv(population.map, col_names = F)
  } else {
    population.map <- population.map
  }
  
  if (is.null(pop.levels)) {
    tfam <- population.map %>%
      select(INDIVIDUALS=X1) %>%
      mutate(
        INDIVIDUALS = as.character(INDIVIDUALS),
        POP_ID = substr(INDIVIDUALS, pop.id.start, pop.id.end),
        POP_ID = factor(POP_ID)
      )
  } else {
    tfam <- population.map %>%
      select(INDIVIDUALS=X1) %>%
      mutate(
        INDIVIDUALS = as.character(INDIVIDUALS),
        POP_ID = substr(INDIVIDUALS, pop.id.start, pop.id.end),
        POP_ID = factor(POP_ID, levels = pop.levels, ordered =T)
      )
  }
  
  tfam <- tfam %>%
    select(POP_ID, INDIVIDUALS) %>%
    arrange(POP_ID, INDIVIDUALS) %>%
    mutate(
      COL1 = rep("0",n()),
      COL2 = rep("0",n()),
      COL3 = rep("0",n()),
      COL4 = rep("0",n())
    )
  
  write.table(tfam, filename, sep = "\t", row.names = F, col.names = F, quote = F)
  
  invisible(cat(sprintf(
    "tfam filename:
  %s\n
  Written in the directory:
  %s",
    filename, getwd()
  )))
}
