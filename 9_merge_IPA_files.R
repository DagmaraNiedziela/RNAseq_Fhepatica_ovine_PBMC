# For this file, I took Supplementary Table with IPA results from Garcia-Campos et al. 2019 as bovine pathways 
# Ovine IPA was performed using UCD Animal Genomics lab license 
# The files were then all merged here for Supplementary Table S7 

# Ovine pathways #### 
library(tidyverse)
library(readxl)
IPA_diff_2WPI <- read_excel("canonical 2 WPI.xls", sheet = 1)
IPA_diff_16WPI <- read_excel("canonical 16 WPI.xls", sheet = 1)
colnames(IPA_diff_2WPI)
# [1] "Ingenuity Canonical Pathways" "-log(p-value)"                "Ratio"                        "z-score"                      "Molecules"   
# [1] "ingenuity_canonical_pathways" "log_p_value"                  "ratio"                        "z_score"                      "molecules" 

library(janitor)
IPA_diff_2WPI <- IPA_diff_2WPI %>% clean_names()
IPA_diff_16WPI <- IPA_diff_16WPI %>% clean_names()
colnames(IPA_diff_16WPI) 

# This will take all pathways, including non-sognificant ones #### 

IPA_diff_join <- full_join(IPA_diff_2WPI, IPA_diff_16WPI, by = "ingenuity_canonical_pathways", suffix = c("_2WPI", "_16WPI"))
View(IPA_diff_join)
write.csv(IPA_diff_join, "IPA_canonical_diff_both_WPI.csv")

IPA_diff_join2 <- IPA_diff_join %>% select(ingenuity_canonical_pathways, log_p_value_2WPI, ratio_2WPI,z_score_2WPI, 
                                           log_p_value_16WPI, ratio_16WPI,z_score_16WPI)

# Bovine pathways #####
Andres_IPA <- read_excel("Andres_cattle_Suppl.xlsx", sheet = 7)
View(Andres_IPA)
colnames(Andres_IPA) <- c("ingenuity_canonical_pathways","log_p_value_1WPI","ratio_1WPI","z_score_1WPI",
                          "log_p_value_14-1WPI","ratio_14-1WPI","z_score_14-1WPI",
                          "log_p_value_14WPI","ratio_14WPI","z_score_14WPI")
Andres_IPA <- Andres_IPA[-1,]

IPA_sheep_cow_join <- full_join(IPA_diff_join2, Andres_IPA, by = "ingenuity_canonical_pathways", suffix = c("_sheep", "_cow"))
View(IPA_sheep_cow_join)
write.csv(IPA_sheep_cow_join, "IPA_sheep_cow_canonical.csv")

# Remove non-significant pathways - cut off p = 0.05 which is 0log P of 1.3 

IPA_diff_2WPI2 <- IPA_diff_2WPI %>% filter(log_p_value > 1.3)
IPA_diff_16WPI2 <- IPA_diff_16WPI %>% filter(log_p_value > 1.3)

IPA_diff_join2a <- full_join(IPA_diff_2WPI2, IPA_diff_16WPI2, by = "ingenuity_canonical_pathways", suffix = c("_2WPI", "_16WPI"))
View(IPA_diff_join2a)
write.csv(IPA_diff_join2a, "IPA_canonical_diff_both_WPI_sig.csv")

IPA_diff_join3 <- IPA_diff_join2a %>% select(ingenuity_canonical_pathways, log_p_value_2WPI, ratio_2WPI,z_score_2WPI, 
                                           log_p_value_16WPI, ratio_16WPI,z_score_16WPI)

Andres_IPA <- read_excel("Andres_cattle_Suppl.xlsx", sheet = 7)
View(Andres_IPA)
colnames(Andres_IPA) <- c("ingenuity_canonical_pathways","log_p_value_1WPI","ratio_1WPI","z_score_1WPI",
                          "log_p_value_14-1WPI","ratio_14-1WPI","z_score_14-1WPI",
                          "log_p_value_14WPI","ratio_14WPI","z_score_14WPI")
Andres_IPA <- Andres_IPA[-1,]

IPA_sheep_cow_join2 <- full_join(IPA_diff_join3, Andres_IPA, by = "ingenuity_canonical_pathways", suffix = c("_sheep", "_cow"))
View(IPA_sheep_cow_join2)
write.csv(IPA_sheep_cow_join2, "IPA_sheep_cow_canonical_sig.csv") 
# This is Table S7 