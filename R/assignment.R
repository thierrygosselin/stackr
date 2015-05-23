# ASSIGNMENT  
# rm(list=ls())
# ls()

# assignment function
assignment_genodive <- function(assignment.lmax, assignment.lhome, lmax.migrant.skip, lmax.number.migrant, lmax.skip, lhome.migrant.skip, lhome.number.migrant, lhome.skip, pop.levels, pop.id.start, pop.id.end, number.individual, number.pop) {

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

# Assignment Figures functions
figure_assignment <- function(assignment.summary) {
  
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

# Figure function : Basic assignment plot
figure_assignment_stacked_bar <- function(assignment.summary) {

  assignment_summary_stacked <- assignment.summary %>%
    mutate(Current = factor(Current, levels = SITES_LEVELS, ordered = T)) %>%
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

