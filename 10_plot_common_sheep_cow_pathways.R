# The common pathways in sheep and cattle were explored by hand and summarised in sheet 2 of the Excel file 
# That sheet 2 was originally to be Table 4 in the paper, but I decided to turn it into a figure 
# This code is for Figure 7 
library(tidyverse)
library(readxl)

# Upload data ####
IPA_for_plot <- read_excel("IPA_sheep_cow_canonical_sig.xlsx", sheet = 2)

View(IPA_for_plot)

IPA_for_plot <- IPA_for_plot[-c(1,5,13,22,23),]

IPA_for_plot$Pathway_overlap_type <- c(rep("Common acute pathways",3),rep("Acute phase in cattle, chronic phase in sheep",7),rep("Common chronic pathways",8))
colnames(IPA_for_plot)
#[1] "ingenuity_canonical_pathways" "log_p_value_2WPI"             "z_score_2WPI"                 "log_p_value_16WPI"            "z_score_16WPI"               
#[6] "log_p_value_1WPI"             "z_score_1WPI"                 "log_p_value_14WPI"            "z_score_14WPI"                "Pathway_overlap_type"   

# Big pivot ####
spec <- tribble(
  ~.name,            ~.value, ~WPI,      
  "log_p_value_2WPI","log_p_value", "2",           
  "log_p_value_16WPI","log_p_value","16",           
  "log_p_value_1WPI","log_p_value","1",            
  "log_p_value_14WPI","log_p_value","14",      
  "z_score_2WPI" ,         "z_score", "2",      
  "z_score_16WPI" ,           "z_score", "16",      
  "z_score_1WPI" ,            "z_score", "1", 
  "z_score_14WPI" ,            "z_score", "14",        
)

IPA_for_plot_long <- IPA_for_plot %>% pivot_longer_spec(spec = spec)

View(IPA_for_plot_long) 

IPA_for_plot_long$WPI <- factor(IPA_for_plot_long$WPI)
levels(IPA_for_plot_long$WPI) 

# New column changing WPI into species and stage 
IPA_for_plot_long$Species_Stage <- IPA_for_plot_long$WPI 
levels(IPA_for_plot_long$Species_Stage)
library(plyr)
IPA_for_plot_long$Species_Stage <- mapvalues(IPA_for_plot_long$Species_Stage, from = c("1","14","16","2"), 
                                       to = c("Cattle_acute", "Cattle_chronic", "Sheep_chronic", "Sheep_acute"))
detach(package:plyr)

levels(IPA_for_plot_long$Species_Stage) <- c("Sheep_acute", "Sheep_chronic", "Cattle_acute", "Cattle_chronic")
IPA_for_plot_long$Species_Stage <- factor(IPA_for_plot_long$Species_Stage, levels = c("Sheep_acute", "Sheep_chronic", "Cattle_acute", "Cattle_chronic"))

# Prep values for plot ####     
colnames(IPA_for_plot_long)

IPA_for_plot_long$z_score <- as.numeric(as.character(IPA_for_plot_long$z_score))
IPA_for_plot_long$log_p_value <- as.numeric(as.character(IPA_for_plot_long$log_p_value)) 

IPA_for_plot_long$Pathway_overlap_type <- factor(IPA_for_plot_long$Pathway_overlap_type) 

# PLOT ####

# This is Figure 7 
ggplot(data = IPA_for_plot_long) + 
  geom_col(mapping = aes(x = ingenuity_canonical_pathways, y = log_p_value, fill = z_score, color = Species_Stage), size = 1, position=position_dodge2(reverse = TRUE)) + 
  scale_fill_distiller(palette = "RdYlBu") +
  scale_color_brewer(palette = "Set1") +
  facet_grid(Pathway_overlap_type ~ ., scales = "free", space = "free") +
  coord_flip() 

?facet_grid

ggsave("IPA_common_pathways_fgrid.jpeg") 

# Other plot versions --- everything from here below are just different explorations for these plots, none have been used in paper 
ggplot(data = IPA_for_plot_long) + 
  geom_col(mapping = aes(x = ingenuity_canonical_pathways, y = log_p_value, fill = z_score, color = Species_Stage), size = 1, position=position_dodge2(reverse = TRUE)) + 
  scale_fill_distiller(palette = "RdYlBu") +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~Pathway_overlap_type, drop = TRUE) +
  coord_flip() 

ggsave("IPA_common_pathways.jpeg") 

ggplot(data = IPA_for_plot_long) + 
  geom_col(mapping = aes(x = ingenuity_canonical_pathways, y = log_p_value, fill = Species_Stage), position=position_dodge2(reverse = TRUE)) + 
  scale_fill_brewer(palette = "Set1") +
  facet_grid(Pathway_overlap_type ~ ., scales = "free", space = "free") +
  coord_flip() 

ggsave("IPA_common_pathways_fgrid_nozscore.jpeg") 

# ** Split by pathway overlap type ####
A <- ggplot(data = subset(IPA_for_plot_long, Pathway_overlap_type == "Common acute pathways")) + 
  geom_col(mapping = aes(x = ingenuity_canonical_pathways, y = log_p_value, fill = z_score, color = Species_Stage), size = 1,position=position_dodge2(reverse = TRUE)) + 
  scale_fill_distiller(palette = "RdYlBu", breaks = c(-1, 3)) +
  scale_color_brewer(palette = "Set1") + 
  ggtitle("Common acute pathways") + 
  coord_flip()  

A
ggsave("IPA_common_pathways_acute.jpeg")

B <- ggplot(data = subset(IPA_for_plot_long, Pathway_overlap_type == "Common chronic pathways")) + 
  geom_col(mapping = aes(x = ingenuity_canonical_pathways, y = log_p_value, fill = z_score, color = Species_Stage), size = 1,position=position_dodge2(reverse = TRUE)) + 
  scale_fill_distiller(palette = "RdYlBu", breaks = c(-1, 3)) +
  scale_color_brewer(palette = "Set1") + 
  ggtitle("Common chronic pathways") + 
  coord_flip() 

B
ggsave("IPA_common_pathways_chronic.jpeg")

C <- ggplot(data = subset(IPA_for_plot_long, Pathway_overlap_type == "Acute phase in cattle, chronic phase in sheep")) + 
  geom_col(mapping = aes(x = ingenuity_canonical_pathways, y = log_p_value, fill = z_score, color = Species_Stage), size = 1,position=position_dodge2(reverse = TRUE)) + 
  scale_fill_distiller(palette = "RdYlBu", breaks = c(-1, 3)) +
  scale_color_brewer(palette = "Set1") + 
  ggtitle("Acute phase in cattle, chronic phase in sheep") + 
  coord_flip() 

C
ggsave("IPA_common_pathways_mix.jpeg")

library(ggpubr)

ggarrange(A,B,C, nrow = 3, common.legend = TRUE, align = "v", heights = c(0.5,1,1)) 

ggsave("Arranged_common_IPA_pathways.jpeg")

# ** Split by pathway overlap type - and remove z score, confusing ####
A2 <- ggplot(data = subset(IPA_for_plot_long, Pathway_overlap_type == "Common acute pathways")) + 
  geom_col(mapping = aes(x = ingenuity_canonical_pathways, y = log_p_value, fill = Species_Stage), position=position_dodge2(reverse = TRUE)) + 
  scale_fill_brewer(palette = "Set1") + 
  ggtitle("Common acute pathways") + 
  coord_flip()  

A2
ggsave("IPA_common_pathways_acute.jpeg")

B2 <- ggplot(data = subset(IPA_for_plot_long, Pathway_overlap_type == "Common chronic pathways")) + 
  geom_col(mapping = aes(x = ingenuity_canonical_pathways, y = log_p_value, fill = Species_Stage), position=position_dodge2(reverse = TRUE)) + 
  scale_fill_brewer(palette = "Set1") + 
  ggtitle("Common chronic pathways") + 
  coord_flip() 

B2
ggsave("IPA_common_pathways_chronic.jpeg")

C2 <- ggplot(data = subset(IPA_for_plot_long, Pathway_overlap_type == "Acute phase in cattle, chronic phase in sheep")) + 
  geom_col(mapping = aes(x = ingenuity_canonical_pathways, y = log_p_value, fill = Species_Stage), position=position_dodge2(reverse = TRUE)) + 
  scale_fill_brewer(palette = "Set1") + 
  ggtitle("Acute phase in cattle, chronic phase in sheep") + 
  coord_flip() 

C2
ggsave("IPA_common_pathways_mix.jpeg")

library(ggpubr)

ggarrange(A2,B2,C2, nrow = 3, common.legend = TRUE, align = "v", heights = c(0.5,1,1)) 

ggsave("Arranged_common_IPA_pathways_nozscore.jpeg") 

ggarrange(A2,B2, nrow = 2, common.legend = TRUE, align = "v", heights = c(0.5,1)) 

ggsave("Arranged_common_IPA_pathways_nozscoreAB.jpeg")


