#' @title Prepare and summarise the assignment results from GENODIVE.
#' @description Import assignment results from GENODIVE. 
#' Current version needs results from both statistics, the home likelihood 
#' and likelihood ratio.
#' @param assignment.lmax The file with likelihood ratio assignment results
#'  from GENODIVE.
#' @param assignment.lhome The file with home likelihood assignment results
#'  from GENODIVE.
#' @param lmax.migrant.skip The number of lines to skip before the migrant
#'  info in the likelihood ratio assignment results. Usually = 11.
#' @param lmax.number.migrant The number of migrant detected
#'  in the likelihood ratio assignment results.
#' @param lmax.skip The number of lines to skip before the individuals info
#' in the likelihood ratio assignment results.
#' @param lhome.migrant.skip The number of lines to skip before 
#' the migrant info in the home likelihood assignment results. Usually = 11.
#' @param lhome.number.migrant The number of migrant detected
#' in the home likelihood assignment results.
#' @param lhome.skip The number of lines to skip before the individuals info 
#' in the home likelihood assignment results.
#' @param pop.levels An optional character string with your populations ordered.
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
#' @param number.individual The number of individuals analysed.
#' @param number.pop The number of populations analysed.
#' @import dplyr
#' @import readr
#' @importFrom stringr str_sub

assignment_genodive <- function(assignment.lmax, assignment.lhome, lmax.migrant.skip, lmax.number.migrant, lmax.skip, lhome.migrant.skip, lhome.number.migrant, lhome.skip, pop.levels, pop.id.start, pop.id.end, number.individual, number.pop) {
  
  Individual <- NULL
  Current <- NULL
  Inferred <- NULL
  Lik_max <- NULL
  Lik_home <- NULL
  Lik_ratio <- NULL
  Threshold <- NULL
  STATUS <- NULL
  IND <- NULL
  N_CORRECTED <- NULL
  DISCRIMINATE <- NULL
  STATISTICS <- NULL
  POUR <- NULL
 
  
  migrants_lmax <- read_delim(
    assignment.lmax,
    delim = "\t",
    skip = lmax.migrant.skip,
    n_max = lmax.number.migrant,
    col_names = c("Individual", "Current", "Inferred", "Lik_max", "Lik_home", "Lik_ratio", "Threshold"),
    progress = interactive(),
    col_types = "ccciiii")
  # View(migrants_lmax)
  
  migrants_lhome <- read_delim(
    assignment.lhome,
    delim = "\t",
    skip = lhome.migrant.skip,
    n_max = lhome.number.migrant,
    col_names = c("Individual", "Current", "Inferred", "Lik_max", "Lik_home", "Lik_ratio", "Threshold"),
    progress = interactive(),
    col_types = "ccciiii")
  # View(migrants_lhome)
  
  # Create a new vector with a unindified pop label
  pop.levels.lhome <- c(pop.levels, "?")
  
  # create a new vector to assign the class of the column during the import
  col.types <- stri_join("ccciii", stri_dup("i", times = number.pop), sep = "") # ccciii are default, integer are added based on the number of populations
  
  assignment_lmax <- read_delim(
    assignment.lmax,
    delim = "\t",
    skip = lmax.skip,
    n_max = number.individual,
    col_names = c("Individual", "Current", "Inferred", "Lik_max", "Lik_home", "Lik_ratio",
                  paste(pop.levels, sep = ",")),
    progress = interactive(),
    col_types = col.types) %>%
    mutate(
      Current = factor(str_sub(Current, pop.id.start, pop.id.end), levels = pop.levels, ordered =T),
      Inferred = factor(str_sub(Inferred, pop.id.start, pop.id.end), levels = pop.levels, ordered =T)
    ) %>%
    mutate(
      STATUS = ifelse(Individual %in% migrants_lmax$Individual, "migrant",
                      ifelse(Current != Inferred, "error", "resident")),
      STATUS = factor(STATUS, levels = c("resident", "migrant", "error"), ordered = T),
      Inferred = ifelse(STATUS == "migrant" & Current == Inferred, "?", as.character(Inferred)),
      Inferred = factor(Inferred, levels = pop.levels.lhome, ordered = T)
    ) %>%
    group_by(Current, Inferred, STATUS) %>%
    summarise(IND = as.numeric(n())) %>%
    group_by(Current) %>%
    mutate(
      N_CORRECTED = as.numeric(sum(IND) - sum(IND[STATUS == "migrant"])),
      POUR = as.numeric(round((IND/N_CORRECTED), 2))
    ) %>%
    ungroup() %>%
    mutate(
      DISCRIMINATE = POUR * 100,
      DISCRIMINATE = ifelse(as.character(Current) == as.character(Inferred), DISCRIMINATE, ""),
      STATISTICS = rep("Likelihood ratio", n())
    )
  
  
  # View(assignment_lmax)
  # names(assignment_lmax)
  # levels(assignment_lmax$STATUS)
  # class(assignment_lmax$Current)
  # class(assignment_lmax$Inferred)
  # levels(assignment_lmax$Current)
  # levels(assignment_lmax$Inferred)
  
  
  assignment_lhome <- read_delim(
    assignment.lhome,
    delim = "\t",
    skip = lhome.skip,
    n_max = number.individual,
    col_names = c("Individual", "Current", "Inferred", "Lik_max", "Lik_home", "Lik_ratio",
                  paste(pop.levels, sep = ",")),
    progress = interactive(),
    col_types = col.types) %>%
    mutate(
      Current = factor(str_sub(Current, pop.id.start, pop.id.end), levels = pop.levels.lhome, ordered =T),
      Inferred = factor(str_sub(Inferred, pop.id.start, pop.id.end), levels = pop.levels.lhome, ordered =T)
    ) %>%
    mutate(
      STATUS = ifelse(Individual %in% migrants_lhome$Individual, "migrant",
                      ifelse(Current != Inferred, "error", "resident")),
      STATUS = factor(STATUS, levels = c("resident", "migrant", "error"), ordered = T),
      Inferred = ifelse(STATUS == "migrant" & Current == Inferred, "?", as.character(Inferred)),
      Inferred = factor(Inferred, levels = pop.levels.lhome, ordered = T)
    ) %>%
    group_by(Current, Inferred, STATUS) %>%
    summarise(IND = as.numeric(n())) %>%
    group_by(Current) %>%
    mutate(
      N_CORRECTED = as.numeric(sum(IND) - sum(IND[STATUS == "migrant"])),
      POUR = as.numeric(round((IND/N_CORRECTED), 2))
    ) %>%
    ungroup() %>%
    mutate(
      DISCRIMINATE = POUR * 100,
      DISCRIMINATE = ifelse(as.character(Current) == as.character(Inferred), DISCRIMINATE, ""),
      STATISTICS = rep("Home likelihood", n())
    )
  # View(assignment_lhome)
  # names(assignment_lhome)
  
  # invert the vector for figure will be easier to inspect
  pop.levels.inverted <- rev(pop.levels)
  
  # Bind the 2 datasets
  assignment.summary <- assignment_lhome %>%
    rbind(assignment_lmax) %>%
    ungroup() %>%
    mutate(
      Current = factor(Current, levels = pop.levels.inverted, ordered = T),
      Inferred = factor(Inferred, levels = pop.levels.lhome, ordered = T),
      STATUS = factor(STATUS, levels = c("resident", "migrant", "error"), ordered = T),
      STATISTICS = factor(STATISTICS, levels = c("Home likelihood", "Likelihood ratio"), ordered = T)
    )
}


#' @title Cleveland dot plot figure of assignement results.
#' @description GGPLOT2 Cleveland dot plot figure of assignment results.
#' @param assignment.summary The assignment summary file created 
#' with assignment_genodive function.

figure_assignment <- function(assignment.summary) {
  
  Individual <- NULL
  Current <- NULL
  Inferred <- NULL
  STATUS <- NULL
  IND <- NULL
  DISCRIMINATE <- NULL
  POUR <- NULL

  ggplot(assignment.summary, aes(x = Inferred, y = Current, size = (POUR*100)))+
    #    geom_jitter(shape = 21, alpha = 0.5, aes(fill = STATUS), position = position_jitter(width = 0.05))+
    geom_point(shape = 21, alpha = 0.5, aes(fill = STATUS))+
    geom_text(aes(y = as.numeric(Current)-sqrt(POUR)/3, label = IND), vjust = 0.9, colour = "grey60", size = 4)+ # with 3 colors
    geom_text(aes(y = as.numeric(Current), label = DISCRIMINATE), vjust = 0.5, colour = "black", size = 4, face = "bold")+
    scale_fill_manual(name = "Assignment", values = c("darkgreen", "blue", "darkred"))+ # with 3 categories
    #  scale_fill_manual(name = "Assignment", values = c("darkgreen", "blue"))+ # with 2 categories
    scale_size_area(guide = FALSE, max_size = 20)+
    labs(x = "Inferred population")+
    labs(y = "Current population")+
    theme_bw()+
    theme(panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(), 
          panel.grid.major.y = element_line(colour = "grey60", linetype = "dashed"), 
          axis.title.x = element_text(size = 12, family = "Helvetica", face = "bold"), 
          axis.text.x = element_text(size = 12, family = "Helvetica", face = "bold"), 
          axis.title.y = element_text(size = 12, family = "Helvetica", face = "bold"), 
          axis.text.y = element_text(size = 12, family = "Helvetica", face = "bold"))+
    facet_wrap(~ STATISTICS, scales = "free_x")
} 




#' @title Classic stacked bar plot figure of assignement results.
#' @description GGPLOT2 stacked bar plot figure of assignment results.
#' @param assignment.summary The assignment summary file created 
#' with assignment_genodive function.
#' @param pop.levels An optional character string with your populations ordered.


figure_assignment_stacked_bar <- function(assignment.summary, pop.levels) {
  
  Current <- NULL
  Inferred <- NULL
  STATUS <- NULL
  IND <- NULL

  
  assignment_summary_stacked <- assignment.summary %>%
    mutate(Current = factor(Current, levels = pop.levels, ordered = T)) %>%
    arrange(STATUS)
  
  ggplot(assignment_summary_stacked, aes(x = Current, y = IND, fill = factor(STATUS)))+
    geom_bar(stat = "identity", position = "fill")+
    scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), expand = c(0, 0))+ #in order for the line not to expand beyond the graph!
    scale_fill_manual(name = "Assignment", 
                      breaks = c("resident", "migrant", "error"),
                      values = c("darkgreen", "blue", "darkred"))+ # with 3 categories
    labs(x = "Sampling sites")+
    labs(y = "Proportion")+
    facet_wrap(~ STATISTICS, nrow = 1, ncol = 2)+
    theme(
      panel.background = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.title.x = element_text(size = 14, family = "Helvetica", face = "bold"),
      axis.title.y = element_text(size = 14, family = "Helvetica", face = "bold"),
      axis.text.x = element_text(size = 12, family = "Helvetica", face = "bold"),
      axis.line = element_line(),
      axis.line.y = element_line(linetype = "solid", size = 0.5, colour = "grey50"),
      axis.line.x = element_blank(),
      legend.title = element_text(size = 12, family = "Helvetica", face = "bold"),
      legend.text = element_text(size = 12, family = "Helvetica", face = "bold"),
      strip.text.y = element_text(angle = 0, size = 12, family = "Helvetica", face = "bold"),
      strip.text.x = element_text(size = 12, family = "Helvetica", face = "bold")
    )
} # Figure function: assignment stacked bar function

