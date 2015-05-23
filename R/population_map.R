## Population map

population_map <- function (data, pop.id.start, pop.id.end, pop.levels) {
  
  pop.map <- read_tsv(data, col_names = F) %>%
    select(INDIVIDUALS=X1) %>%
    mutate(
      INDIVIDUALS = as.character(INDIVIDUALS),
      POP_ID = str_sub(INDIVIDUALS, pop.id.start, pop.id.end),
      POP_ID = factor(POP_ID, levels = pop.levels, ordered =T)
    ) %>%
    arrange(INDIVIDUALS, POP_ID)
  
  invisible(cat(sprintf("This is the order of your POPULATION MAP file %s",
                        list(levels(pop.map$POP_ID))
                        )
                )
            )
  
pop.map
}
