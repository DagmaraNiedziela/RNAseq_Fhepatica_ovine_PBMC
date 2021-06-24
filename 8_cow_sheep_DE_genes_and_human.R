# This is a comparison of DE genes with bovine PBMC response to Fasciola hepatica - paper by Garcia-Campos et al. 2019, Front Immunol 
# Longitudinal analysis was used for this comparison - "difference" DE gene data frames 
# Dagmara Niedziela

# Convert sheep gene names to human #### 
y2WPI_diff <- res4_LFC_diff_2WPI_sigDF %>% dplyr::select(Row.names) 

y16WPI_diff <- res4_LFC_diff_16WPI_sigDF %>% dplyr::select(Row.names) 

head(res4_LFC_diff_2WPI_sigDF)

# Reminder of ovine longitudinal DE genes file 
# ** Make a DE genes file 
sheep_diff_DE_genes <- full_join(res4_LFC_diff_2WPI_sigDF,res4_LFC_diff_16WPI_sigDF, by = "Row.names", suffix = c("_2WPI","_16WPI"))
View(sheep_diff_DE_genes)
sheep_diff_DE_genes <- sheep_diff_DE_genes %>% select(Row.names,starts_with("log2"),starts_with("padj"))
head(sheep_diff_DE_genes)
sheep_diff_DE_genes <- sheep_diff_DE_genes[,c(1,2,4,3,5)]
write.csv(sheep_diff_DE_genes,"sheep_diff_DE_genes.csv") 

# Annotate the sheep diff DE file - merge with biotype 
head(biotype2) # this is the one 
head(biotype3)
sheep_diff_DE_genes_annot <- merge(sheep_diff_DE_genes,biotype2,by.x=1,by.y=0, all=FALSE)
head(sheep_diff_DE_genes_annot)

library(dplyr)
sheep_diff_DE_genes_annot <- sheep_diff_DE_genes_annot %>% select(-strand,-gene_biotype)
write.csv(sheep_diff_DE_genes_annot,"sheep_diff_DE_genes_annot.csv")

# Convert to human - biomart 
sheep_human_conversion <- read.table(file="mart_export_sheep_ramb_to_human.txt", header=TRUE, sep = "\t")
head(sheep_human_conversion)
nrow(sheep_human_conversion) # 46768
View(sheep_human_conversion)
# There are some many2many orthologs, I will only keep the first one... 

library(janitor)
sheep_human_conversion <- sheep_human_conversion %>% clean_names()
colnames(sheep_human_conversion) 
[1] "gene_stable_id"                          "gene_stable_id_version"                 
[3] "transcript_stable_id"                    "transcript_stable_id_version"           
[5] "human_gene_stable_id"                    "human_gene_name"                        
[7] "human_homology_type"                     "human_orthology_confidence_0_low_1_high"

sheep_human_conversion$Row.names <- sheep_human_conversion$gene_stable_id

library(tidyverse)
sheep_human_conversion <- sheep_human_conversion%>% distinct(gene_stable_id, .keep_all = TRUE) 

nrow(sheep_human_conversion) # 26478

# join samples by the Ensembl ovine name 
sheep_DE_human <- left_join(sheep_diff_DE_genes, sheep_human_conversion, by = "Row.names")
head(sheep_DE_human)
nrow(sheep_DE_human) #2753 
View(sheep_DE_human) 

# 2 WPI table 
sheep_DE_human_2WPI <- sheep_DE_human %>% filter(log2FoldChange_2WPI != "<NA>") %>% select(human_gene_stable_id,log2FoldChange_2WPI,padj_2WPI)
sheep_DE_human %>% filter(log2FoldChange_2WPI != "<NA>")

sheep_DE_human_2WPI
nrow(sheep_DE_human_2WPI) #1148
sheep_DE_human_2WPI_NA <- sheep_DE_human_2WPI %>% filter(str_detect(human_gene_stable_id,'ENSG'))
sheep_DE_human_2WPI_NA 
nrow(sheep_DE_human_2WPI_NA) #1082

# 16 WPI table 
sheep_DE_human_16WPI <- sheep_DE_human %>% filter(log2FoldChange_16WPI != "<NA>") %>% select(human_gene_stable_id,log2FoldChange_16WPI,padj_16WPI)
sheep_DE_human %>% filter(log2FoldChange_16WPI != "<NA>")

sheep_DE_human_16WPI
nrow(sheep_DE_human_16WPI) #1927
sheep_DE_human_16WPI_NA <- sheep_DE_human_16WPI %>% filter(str_detect(human_gene_stable_id,'ENSG'))
sheep_DE_human_16WPI_NA 
nrow(sheep_DE_human_16WPI_NA) #1815

# Annotation DBI convert bovine entrez to Ensembl #### 
# Import Andres's data 
library(readxl)
cow_DE_genes <- read_excel("GarciaCampos_DE.xlsx", sheet = "Table S5")
# Because there are minuses instead of NAs, I previously deleted all - in Excel, so all negative LFC
# Need to redo Excel file by hand and recreate files 
View(cow_DE_genes)
nrow(cow_DE_genes) #2020 
table(is.na(cow_DE_genes$Ensembl_name))
# FALSE  TRUE 
#1844   176 
176/2020*100 # 8.712871 % genes that were not converted 

library(AnnotationDbi)
library(OrganismDbi)
library(org.Hs.eg.db)
library("org.Bt.eg.db") 
keytypes(org.Bt.eg.db)

cow_DE_genes
cow_DE_genes$EntrezID <- as.character(cow_DE_genes$EntrezID)
cow_DE_genes$Ensembl_name <- mapIds(org.Bt.eg.db,
                                       keys=cow_DE_genes$EntrezID,
                                       column="ENSEMBL",
                                       keytype="ENTREZID",
                                       multiVals="first")

# Convert bovine to human #####

cow_human_conversion <- read.table(file="mart_export_cow_to_human.txt", header=TRUE, sep = "\t")
head(cow_human_conversion)
nrow(cow_human_conversion) # 48827 
View(cow_human_conversion)
# There are some many2many orthologs, I will only keep the first one... 

library(janitor)
cow_human_conversion <- cow_human_conversion %>% clean_names()
colnames(cow_human_conversion) 
[1] "gene_stable_id"                          "gene_stable_id_version"                 
[3] "transcript_stable_id"                    "transcript_stable_id_version"           
[5] "human_gene_stable_id"                    "human_gene_name"                        
[7] "human_orthology_confidence_0_low_1_high" "human_homology_type" 

cow_human_conversion$Ensembl_name <- cow_human_conversion$gene_stable_id

library(tidyverse)
cow_human_conversion <- cow_human_conversion%>% distinct(gene_stable_id, .keep_all = TRUE) 

nrow(cow_human_conversion) # 27607 

# join samples by the bovine Ensembl gene name 
Andres_DE_human <- left_join(cow_DE_genes, cow_human_conversion, by = "Ensembl_name")
head(Andres_DE_human)
nrow(Andres_DE_human) #2020 
View(Andres_DE_human) 

# Convert 

# 1 WPI table 
Andres_1WPI <- Andres_DE_human %>% filter(!is.na(`1WPI_logFC`)) %>% select(human_gene_stable_id,Ensembl_name,gene_symbol,`1WPI_logFC`,`1WPI_FDR`)
Andres_DE_human %>% filter(!is.na(`1WPI_logFC`))

Andres_1WPI
nrow(Andres_1WPI) #21
Andres_1WPI_NA <- Andres_1WPI %>% filter(str_detect(human_gene_stable_id,'ENSG'))
Andres_1WPI_NA 
nrow(Andres_1WPI_NA) #17

# 1 WPI human names vector 
Andres_1WPI_NA$human_gene_stable_id

# 14 WPI table 
Andres_14WPI <- Andres_DE_human %>% filter(!is.na(`14WPI_logFC`)) %>% select(human_gene_stable_id,Ensembl_name,gene_symbol,`14WPI_logFC`,`14WPI_FDR`)
Andres_DE_human %>% filter(!is.na(`14WPI_logFC`))

Andres_14WPI
nrow(Andres_14WPI) #1624
Andres_14WPI_NA <- Andres_14WPI %>% filter(str_detect(human_gene_stable_id,'ENSG'))
Andres_14WPI_NA 
nrow(Andres_14WPI_NA) #1420 

# 14 WPI names vector 
Andres_14WPI_NA$human_gene_stable_id

# Join cow and sheep datasets - human name, LFC, P value & annotate with gene symbol and name #### 

# Human gene attributes file 
human_attr <- read.table(file="mart_export_human_gene_attr.txt", header=TRUE, sep = "\t", fill = TRUE)
head(human_attr)
nrow(human_attr) # 57865
View(human_attr)
# There are some many2many orthologs, I will only keep the first one... 

library(janitor)
human_attr <- human_attr %>% clean_names() 

# Remove duplicates 
human_attr <-  human_attr %>% distinct(gene_stable_id, .keep_all = TRUE) 
nrow(human_attr) #37540

human_attr$human_gene_stable_id <- human_attr$gene_stable_id

# Join cow and sheep 
diff_cow_sheep_DE <- full_join(Andres_DE_human,sheep_DE_human, by = "human_gene_stable_id", suffix = c("_cow","_sheep"))
View(diff_cow_sheep_DE)
colnames(diff_cow_sheep_DE)

diff_cow_sheep_DE <- diff_cow_sheep_DE %>% select(human_gene_stable_id,`1WPI_logFC`,`1WPI_FDR`,
                                                  `14WPI_logFC`,`14WPI_FDR`,log2FoldChange_2WPI,
                                                  padj_2WPI,log2FoldChange_16WPI,padj_16WPI,
                                                  gene_stable_id_cow,gene_stable_id_sheep)

# Annotate this file 
diff_cow_sheep_DE_annot <- left_join(diff_cow_sheep_DE,human_attr,by = "human_gene_stable_id")
View(diff_cow_sheep_DE_annot)

# Annotate with AnnotationDbi instead
library(AnnotationDbi)
library(OrganismDbi)
library(org.Hs.eg.db)
library("org.Bt.eg.db") 
columns(org.Hs.eg.db)

diff_cow_sheep_DE_annot <- diff_cow_sheep_DE
diff_cow_sheep_DE_annot$alias <- mapIds(org.Hs.eg.db,
                                    keys=diff_cow_sheep_DE_annot$human_gene_stable_id,
                                    column="ALIAS",
                                    keytype="ENSEMBL",
                                    multiVals="first")
diff_cow_sheep_DE_annot$symbol <- mapIds(org.Hs.eg.db,
                                        keys=diff_cow_sheep_DE_annot$human_gene_stable_id,
                                        column="SYMBOL",
                                        keytype="ENSEMBL",
                                        multiVals="first")
# Alias is different from that on the heatmap - better to merge with biotype.... 
diff_cow_sheep_DE_annot$gene_name <- mapIds(org.Hs.eg.db,
                                        keys=diff_cow_sheep_DE_annot$human_gene_stable_id,
                                        column="GENENAME",
                                        keytype="ENSEMBL",
                                        multiVals="first")

diff_cow_sheep_DE_annot <- apply(diff_cow_sheep_DE_annot,2,as.character) 
diff_cow_sheep_DE_annot <- as.data.frame(diff_cow_sheep_DE_annot)
View(diff_cow_sheep_DE_annot)

# ** Common acute genes table #### 
diff_cow_sheep_common_acute <- left_join(ycs_2_2,diff_cow_sheep_DE_annot, by = "human_gene_stable_id")
View(diff_cow_sheep_common_acute)
write.csv(diff_cow_sheep_common_acute, "diff_cow_sheep_common_acute.csv")

# ** Common chronic genes table ####
diff_cow_sheep_common_chronic <- left_join(ycs_16_16,diff_cow_sheep_DE_annot, by = "human_gene_stable_id")
View(diff_cow_sheep_common_chronic)
write.csv(diff_cow_sheep_common_chronic, "diff_cow_sheep_common_chronic.csv")