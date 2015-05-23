# Choisir l'espace de travail
setwd("/Users/thierry/Dropbox/brook_charr_pop/02_admixture")

data <- "admixture.10pop.10.Q"
population.map <- "population.map.10pop.adu.txt"
pop.id.start <- 5
pop.id.end <- 7
pop.levels <- c("SKY", "LIM", "TWE" , "NIN", "CNC", "MOO", "SUN", "GOO", "WEI", "FRE")

figure_admixture <- function(data, population.map, pop.id.start, pop.id.end, pop.levels) {
  
  if (is.vector(population.map) == "TRUE") {
    population.map <- read_tsv(population.map, col_names = F) %>%
      select(INDIVIDUALS=X1) %>%
      mutate(
        INDIVIDUALS = as.character(INDIVIDUALS),
        POP_ID = str_sub(INDIVIDUALS, pop.id.start, pop.id.end),
        POP_ID = factor(POP_ID, levels = pop.levels, ordered =T)
        )
  } else {
    population.map <- population.map %>%
      select(INDIVIDUALS=X1) %>%
      mutate(
        INDIVIDUALS = as.character(INDIVIDUALS),
        POP_ID = str_sub(INDIVIDUALS, pop.id.start, pop.id.end),
        POP_ID = factor(POP_ID, levels = pop.levels, ordered =T)
        )
  }
  
  admixture <- population.map %>%
    bind_cols(
      read_delim(data,
                 delim = " ",
                 na = "NA",
                 progress = interactive(),
                 col_names = F
                 )
      )
  
  colnames(admixture) <- c("IND","POP", paste ("K", sep = "", seq(from = 1, to = ncol(admixture) - 2, by = 1)))

  admixture <- admixture  %>% 
    melt(id.vars = c("POP", "IND"), variable.name = "ANCESTRY", value.name = "VALUE")

  pop.name <- pop.levels 
  
figure.admixture <- ggplot(admixture, aes(x = IND, y = VALUE, fill = ANCESTRY)) + 
    geom_bar(stat = "identity", position = "fill") + 
#     scale_x_discrete(labels = pop.name, breaks = c()) +
#     scale_fill_manual(name = "Ancestry K", values = colour_palette_ancestry, limits = c("K1", "K2", "K3", "K4")) + 
    scale_y_continuous(expand = c(0, 0)) + labs(y = "Ancestry") + labs(x = "Individuals")
+ 
#     geom_vline(xintercept = 44, color = "black") + 
#     geom_vline(xintercept = 90, color = "black") + 
#     geom_vline(xintercept = 138, color = "black") + 
#     geom_vline(xintercept = 160, color = "black", linetype = "longdash") + 
#     geom_vline(xintercept = 206, color = "black", linetype = "longdash") + 
#     geom_vline(xintercept = 253, color = "black") + 
#     geom_vline(xintercept = 282, color = "black", linetype = "longdash")
figure.admixture

levels(admixture_long$ANCESTRY)

ggsave("admixture_neutral_4K_test.pdf",width=17,height=10,dpi=600,units="cm",useDingbats=F)
  
  
  
  
IND_ORDERED <- as.vector(id)
IND_ORDERED <- as.matrix(id)
IND_ORDERED <- as.list(id)
names(admixture_long)

class(admixture_long$ANCESTRY)
levels(admixture_long$ANCESTRY)
class(admixture$IND)
levels(admixture$IND)

####K=4 with neutral set
#utiliser admixture en format wide
BUR <- subset(admixture,subset=(POP=="BUR"))
BUR <- arrange(BUR,-K2,K3,K1,K4)
BUR <- subset(BUR,select=c("IND"))

GRA <- subset(admixture,subset=(POP=="GRA"))
GRA <- arrange(GRA,-K2,K3,K1,K4)
GRA <- subset(GRA,select=c("IND"))

GUL <- subset(admixture,subset=(POP=="GUL"))
GUL <- arrange(GUL,K3,K2,K4,K1)
GUL <- subset(GUL,select=c("IND"))

LLI <- subset(admixture,subset=(POP=="LLI"))
LLI <- arrange(LLI,K4,K1,K2,K3)
LLI <- subset(LLI,select=c("IND"))

ANG <- subset(admixture,subset=(POP=="ANG"))
ANG <- arrange(ANG,K4,K1,K2,K3)
ANG <- subset(ANG,select=c("IND"))

WEI <- subset(admixture,subset=(POP=="WEI"))
WEI <- arrange(WEI,-K4,K1,K2,K3)
WEI <- subset(WEI,select=c("IND"))

HAY <- subset(admixture,subset=(POP=="HAY"))
HAY <- arrange(HAY,K1,K4,K2,K3)
HAY <- subset(HAY,select=c("IND"))

GOD <- subset(admixture,subset=(POP=="GOD"))
GOD <- arrange(GOD,K1,K4,K2,K3)
GOD <- subset(GOD,select=c("IND"))

# CHU <- subset(admixture,subset=(POP=="CHU"))
# CHU <- arrange(CHU,K3,K1,K2,K4)
# CHU <- subset(CHU,select=c("IND"))

# ID_ORDERED <- rbind(BUR,GRA,GUL,LLI,ANG,WEI,HAY,GOD,CHU)
ID_ORDERED <- rbind(BUR,GRA,GUL,LLI,ANG,WEI,HAY,GOD)
ID_ORDERED <- as.vector(ID_ORDERED)
ID_ORDERED <- as.matrix(ID_ORDERED)
ID_ORDERED <- as.list(ID_ORDERED)


# revenir au dataframe admixture en format long
class(admixture_long$IND)
levels(admixture_long$IND)
admixture_long[["IND"]] <- factor(admixture_long[["IND"]],levels=ID_ORDERED,ordered=T)
levels(admixture_long$IND)

# colour_palette_barrier <- c("#318959","#26C1C5","#29F82F","#69B9FB","#8837DD","#083CF6","#FDA327","#FFFA38")
# colour_palette_barrier <- c("#318959","#26C1C5","#29F82F","#69B9FB","#8837DD","#083CF6","#FDA327","#FFFA38","RED")
colour_palette_ancestry <- c("#8837DD","#318959","#FDA327","#29F82F")



