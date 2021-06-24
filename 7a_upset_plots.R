# Upset plots 
# install.packages("UpSetR")
library(UpSetR)
library(tidyverse)
?upset

# Cow and sheep genes up and down ##### 

# Get data together - gene rows and 1s for where the genes are 
library(tidyverse)
head(sheep_DE_human_2WPI_NA)

y2WPI_sheep_up <- sheep_DE_human_2WPI_NA %>% filter(log2FoldChange_2WPI > 0) %>% dplyr::select(human_gene_stable_id) 
#y2WPI_sheep_up <- as.vector(y2WPI_sheep_up)
nrow(y2WPI_sheep_up) #666

y2WPI_sheep_down <- sheep_DE_human_2WPI_NA %>% filter(log2FoldChange_2WPI < 0) %>% dplyr::select(human_gene_stable_id) 
#y2WPI_sheep_down <- as.vector(y2WPI_sheep_down)
nrow(y2WPI_sheep_down) #416

head(sheep_DE_human_16WPI_NA)
y16WPI_sheep_up <- sheep_DE_human_16WPI_NA %>% filter(log2FoldChange_16WPI > 0) %>% dplyr::select(human_gene_stable_id) 
#y16WPI_sheep_up <- as.vector(y16WPI_sheep_up)
nrow(y16WPI_sheep_up) #1025

y16WPI_sheep_down <- sheep_DE_human_16WPI_NA %>% filter(log2FoldChange_16WPI < 0) %>% dplyr::select(human_gene_stable_id) 
#y16WPI_sheep_down <- as.vector(y16WPI_sheep_down)
nrow(y16WPI_sheep_down) #790

head(Andres_1WPI_NA)
y2WPI_cow_up <- Andres_1WPI_NA %>% filter(`1WPI_logFC` > 0) %>% dplyr::select(human_gene_stable_id) 
#y2WPI_cow_up <- as.vector(y2WPI_cow_up)
nrow(y2WPI_cow_up) #3

y2WPI_cow_down <- Andres_1WPI_NA %>% filter(`1WPI_logFC` < 0) %>% dplyr::select(human_gene_stable_id) 
#y2WPI_cow_down <- as.vector(y2WPI_cow_down)
nrow(y2WPI_cow_down) #14

head(Andres_14WPI_NA)
y16WPI_cow_up <- Andres_14WPI_NA %>% filter(`14WPI_logFC` > 0) %>% dplyr::select(human_gene_stable_id) 
#y16WPI_cow_up <- as.vector(y16WPI_cow_up)
nrow(y16WPI_cow_up) #741

y16WPI_cow_down <- Andres_14WPI_NA %>% filter(`14WPI_logFC` < 0) %>% dplyr::select(human_gene_stable_id) 
#y16WPI_cow_down <- as.vector(y16WPI_cow_down)
nrow(y16WPI_cow_down) #679

# Make 1s 
y2WPI_sheep_up$sheep_2WPI_up <- rep(1,nrow(y2WPI_sheep_up))
head(y2WPI_sheep_up)
y2WPI_sheep_down$sheep_2WPI_down <- rep(1,nrow(y2WPI_sheep_down))
y16WPI_sheep_up$sheep_16WPI_up <- rep(1,nrow(y16WPI_sheep_up))
y16WPI_sheep_down$sheep_16WPI_down <- rep(1,nrow(y16WPI_sheep_down)) 

y2WPI_cow_up$cow_2WPI_up <- rep(1,nrow(y2WPI_cow_up))
head(y2WPI_cow_up)
y2WPI_cow_down$cow_2WPI_down <- rep(1,nrow(y2WPI_cow_down))
y16WPI_cow_up$cow_16WPI_up <- rep(1,nrow(y16WPI_cow_up))
y16WPI_cow_down$cow_16WPI_down <- rep(1,nrow(y16WPI_cow_down))

# Join 
for_upset_species <- full_join(y2WPI_sheep_up,y2WPI_sheep_down, by = "human_gene_stable_id")
head(for_upset_species)
for_upset_species <- full_join(for_upset_species,y16WPI_sheep_up, by = "human_gene_stable_id")
for_upset_species <- full_join(for_upset_species,y16WPI_sheep_down, by = "human_gene_stable_id")
for_upset_species <- full_join(for_upset_species,y2WPI_cow_up, by = "human_gene_stable_id")
for_upset_species <- full_join(for_upset_species,y2WPI_cow_down, by = "human_gene_stable_id")
for_upset_species <- full_join(for_upset_species,y16WPI_cow_up, by = "human_gene_stable_id")
for_upset_species <- full_join(for_upset_species,y16WPI_cow_down, by = "human_gene_stable_id") 

colnames(for_upset_species)
# [1] "human_gene_stable_id" "sheep_2WPI_up"        "sheep_2WPI_down"      "sheep_16WPI_up"       "sheep_16WPI_down"     "cow_2WPI_up"         
# [7] "cow_2WPI_down"        "cow_16WPI_up"         "cow_16WPI_down" 

colnames(for_upset_species) <- c("human_gene_stable_id","sheep_2WPI_up","sheep_2WPI_down",
                                 "sheep_16WPI_up","sheep_16WPI_down","cattle_1WPI_up",         
                                 "cattle_1WPI_down","cattle_14WPI_up","cattle_14WPI_down")

nrow(for_upset_species) #3705
#for_upset_species <- for_upset_species %>% select(-X.x,-X.y,-X.x.x,-X.y.y)
View(for_upset_species)
for_upset_species[is.na(for_upset_species)] <- 0 

ncol(for_upset_species) #9

# Remove NAs and make zeros 

# ** upset plot ####
library(UpSetR)
plot_up_species <- upset(
  for_upset_species,
  # Careful: we have to manually specify how many subsets 
  # we want to display (we can also specify them by column name),
  # and also set nintersects to NA if we are interested in all (non-empty)
  # intersections. It is still a bit fuzzy when it has to display empty intersections
  # on so many sets, but there is an option to do that, if wanted.
  nsets = ncol(for_upset_species) - 1,
  nintersects = NA,
  # Display them from the most numerous intersection to the least
  order.by = "freq",
  line.size = 1.2,
  point.size = 3.5,
  text.scale = 2
)
plot_up_species 

# Highlight specific columns 
plot_up_species2 <-upset(
  for_upset_species,
  # Careful: we have to manually specify how many subsets 
  # we want to display (we can also specify them by column name),
  # and also set nintersects to NA if we are interested in all (non-empty)
  # intersections. It is still a bit fuzzy when it has to display empty intersections
  # on so many sets, but there is an option to do that, if wanted.
  nsets = ncol(for_upset) - 1,
  nintersects = NA,
  # Display them from the most numerous intersection to the least
  order.by = "freq",
  line.size = 1.2,
  point.size = 3.5,
  text.scale = 2, 
  queries = list(list(query = intersects, params = list(c("cattle_14WPI_down","sheep_16WPI_up")), color= "blue"),
            (list(query = intersects, params = list(c("cattle_14WPI_up","sheep_16WPI_down")), color= "purple")),
             (list(query = intersects, params = list(c("sheep_2WPI_up","cattle_14WPI_down")), color= "deeppink")),
            (list(query = intersects, params = list(c("sheep_2WPI_down","cattle_14WPI_up")), color= "magenta"))
            )
                      
  )
plot_up_species2

# Writing to file 
# pretty plot to file 
library(Cairo)
tiff(filename = "Upset_plot_species.tiff", compression = "lzw", width = 30, height = 22, 
     units = "cm", res = 600);
plot_up_species;
dev.off(); 

Cairo(file="Upset_plot_species.png", 
      type="png",
      units="in", 
      width=15, 
      height=12, 
      pointsize=12, 
      dpi=600)
plot_up_species
dev.off() 

tiff(filename = "Upset_plot_species2.tiff", compression = "lzw", width = 30, height = 22, 
     units = "cm", res = 600);
plot_up_species2;
dev.off(); 

# pretty plot to file 
Cairo(file="Upset_plot_species2.png", 
      type="png",
      units="in", 
      width=15, 
      height=12, 
      pointsize=12, 
      dpi=600)
plot_up_species2
dev.off()