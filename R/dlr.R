#' @title Genotype likelihood ratio distance (Dlr)
#' @description The function computes Paetkau's et al. (1997) genotype likelihood
#' ratio distance (Dlr).
#' @param assignment The output assignment file (home likelihood or
#' likelihood ratio statistics) from GENODIVE.
#' @param l.skip The number of lines to skip before the individuals info
#' in GenoDive assignment results.
#' @param number.individuals The number of individuals analysed.
#' @param number.pop The number of populations analysed.
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
#' @param sites.levels An optional character string with your sites names in 
#' the same order as the pop levels below.
#' @param pop.labels The pop label to identify your sampling sites. Use the same
#' string as sites.levels if you are conducted the assignment test on sampling
#' sites that differ in naming with your populations.
#' @param pop.levels An optional character string with your populations ordered.
#' @import dplyr
#' @import readr
#' @import lazyeval
#' @importFrom stringr str_sub
#' @export 
#' @rdname dlr



dlr <- function (assignment, l.skip, number.individuals, number.pop, pop.id.start, pop.id.end, sites.levels, pop.labels, pop.levels) {
  
  Individuals <- NULL
  Current <- NULL
  Inferred <- NULL
  Lik_max <- NULL
  Lik_home <- NULL
  Lik_ratio <- NULL
  Populations <- NULL
  RATIO <- NULL
  DLR <- NULL
  DLR_RELATIVE <- NULL
  
  
  
  # create a new vector to assign the class of the column during the import
  # ccciii are default, integer are added based on the number of populations
  col.types <- stri_join("ccciii", stri_dup("i", times = number.pop), sep = "") 
  
  # import and modify the assignment file form GenoDive
  assignment <- read_delim(
    assignment,
    delim = "\t",
    skip = l.skip,
    n_max = number.individuals,
    col_names = c("Individuals", "Current", "Inferred", "Lik_max", "Lik_home", "Lik_ratio",
                  paste(pop.levels, sep = ",")),
    progress = interactive(),
    col_types = col.types) %>%
    mutate(
      Populations = str_sub(Individuals, pop.id.start, pop.id.end),
      Populations = factor(stri_replace_all_fixed(Populations, sites.levels, pop.labels, vectorize_all = F), levels = pop.levels, ordered =T),
      Populations = droplevels(Populations)
    ) %>%
    select(-c(Current, Inferred, Lik_max, Lik_home, Lik_ratio))
    
    # Dlr relative for one combination of pop
    dlr.relative <- function(pop1, pop2){
      
      dlr <- assignment %>%
        filter_(interp(~ Populations == as.name(pop1) | Populations == as.name(pop2))) %>%
        group_by(Individuals) %>%
        mutate_(
          RATIO1 = interp(~pop1 - pop2, pop1 = as.name(pop1), pop2 = as.name(pop2)),
          RATIO2 = interp(~pop2 - pop1, pop1 = as.name(pop1), pop2 = as.name(pop2)),
          RATIO = interp(~ifelse(Populations == pop1, c.RATIO1, c.RATIO2), Populations = quote(Populations), pop1 = as.name("pop1"), c.RATIO1 = quote(RATIO1), c.RATIO2 = quote(RATIO2))) %>% 
        group_by(Populations) %>%
        summarise(DLR_RELATIVE = (sum(RATIO)/length(RATIO)^2)) %>%
        ungroup %>%
        summarise(DLR_RELATIVE = sum(DLR_RELATIVE)/2)
      return(dlr)
    }
  
  # All combination of populations
  pop.pairwise <- combn(pop.levels, 2)
  pop.pairwise <- matrix(pop.pairwise,nrow=2)
  
  # Dlr for all pairwise populations
  dlr.all.pop <- as.numeric()
  for(i in 1:ncol(pop.pairwise)){
    dlr.all.pop[i] <- dlr.relative(pop1 = pop.pairwise[1,i], pop2 = pop.pairwise[2,i])
  }
  dlr.all.pop <- as.numeric(dlr.all.pop)
  
  # Table with Dlr
  names.pairwise <- combn(pop.levels, 2, paste, collapse = '-')
  
  dlr.table <- data_frame(PAIRWISE_POP = names.pairwise, DLR = dlr.all.pop) %>%
    mutate(DLR = round(as.numeric(DLR), 2))
  
  # Matrix
  dlr.dist <- dist(1:length(pop.levels))
  dlr.matrix <- dlr.all.pop
  attributes(dlr.matrix) <- attributes(dlr.dist)
  dlr.matrix <- as.matrix(dlr.matrix)
  colnames(dlr.matrix) <- rownames(dlr.matrix) <- pop.levels
  dlr.matrix <- as.dist(dlr.matrix)
  
  # Results
  dlr.results.list <- list()
  dlr.results.list$dlr.matrix <- dlr.matrix
  dlr.results.list$dlr.table <- dlr.table
  return(dlr.results.list)  
}

# Dlr absolute
# dlr.absolute <- assignment %>%
#   group_by(Individuals) %>%
#   mutate_(
#     RATIO1 = interp(~pop1 - pop2, pop1 = as.name(pop1), pop2 = as.name(pop2)),
#     RATIO2 = interp(~pop2 - pop1, pop1 = as.name(pop1), pop2 = as.name(pop2)),
#     RATIO = interp(~ifelse(Populations == pop1, c.RATIO1, c.RATIO2), Populations = quote(Populations), pop1 = as.name("pop1"), c.RATIO1 = quote(RATIO1), c.RATIO2 = quote(RATIO2))) %>% 
#     group_by(Populations) %>%
#     summarise(DLR_ABSOLUTE = sum(RATIO)/length(RATIO)) %>%
#     ungroup %>%
#     summarise(DLR_ABSOLUTE = sum(DLR_ABSOLUTE)/2)
