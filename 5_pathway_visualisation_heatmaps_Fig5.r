# Make LOGS tables of my significant genes ##### 
# LOGS is a data frame that has log2 fold changes only and ENSEMBL Gene IDs as row names 
# Use tidyverse 
library(DESeq2)
library(tidyverse)

# Infected vs Control analysis P < 0.1 
res_0WPI_sig01_DF <- as.data.frame(subset(res_0WPI_ann, padj < 0.1)) 
res_2WPI_sig01_DF <- as.data.frame(subset(res_2WPI_ann, padj < 0.1)) 
res_16WPI_sig01_DF <- as.data.frame(subset(res_16WPI_ann, padj < 0.1)) 

# Always take the converted human IDs 
head(res_2WPI_sig01_DF)
nrow(res_2WPI_sig01_DF) # 453 
colnames(res_2WPI_sig01_DF) 
#[1] "Row.names"      "baseMean"       "log2FoldChange" "lfcSE"          "pvalue"         "padj"          
#[7] "strand"         "gene_biotype"   "gene_name"      "Human_name" 

logs_2WPI_P01 <- res_2WPI_sig01_DF %>% select(Human_name, log2FoldChange) 
head(logs_2WPI_P01)
logs_0WPI_P01 <- res_0WPI_sig01_DF %>% select(Human_name, log2FoldChange) 
logs_16WPI_P01 <- res_16WPI_sig01_DF %>% select(Human_name, log2FoldChange) 

# ** Join all log2fold changes ####
logs_P01 <- full_join(logs_0WPI_P01, logs_2WPI_P01, by="Human_name", suffix = c("_0WPI", "_2WPI"))
logs_P01 <- full_join(logs_P01, logs_16WPI_P01, by = "Human_name", suffix = c("_2WPI", "_16WPI"))

head(logs_P01) 
nrow(logs_P01) #1096 
# Remove all the ones with no gene names 
logs_P01 <- logs_P01 %>% filter(!Human_name == "NULL")
nrow(logs_P01) # 433 

colnames(logs_P01) <- c("Human_name", "log2FoldChange_0WPI", "log2FoldChange_2WPI", "log2FoldChange_16WPI")
write.csv(logs_P01, "logs_P01.csv") # don't save this 

# ** HGNC logs ####
logs_2WPI_P01_hg <- res_2WPI_sig01_DF %>% select(gene_name, log2FoldChange) 
head(logs_2WPI_P01_hg)
logs_0WPI_P01_hg <- res_0WPI_sig01_DF %>% select(gene_name, log2FoldChange) 
logs_16WPI_P01_hg <- res_16WPI_sig01_DF %>% select(gene_name, log2FoldChange) 

# ** Join all log2fold changes ####
logs_P01_hg <- full_join(logs_0WPI_P01_hg, logs_2WPI_P01_hg, by="gene_name", suffix = c("_0WPI", "_2WPI"))
logs_P01_hg <- full_join(logs_P01_hg, logs_16WPI_P01_hg, by = "gene_name", suffix = c("_2WPI", "_16WPI"))

head(logs_P01_hg) 
View(logs_P01_hg)
nrow(logs_P01_hg) #1096 
# Remove all the ones with no gene names 
logs_P01_hg <- logs_P01_hg %>% filter(!gene_name == "<NA>")
nrow(logs_P01_hg) # 436 

#colnames(logs_P01) <- c("HGNC_name", "log2FoldChange_0WPI", "log2FoldChange_2WPI", "log2FoldChange_16WPI")
#write.csv(logs_P01, "logs_P01.csv") # don't save this 

# Create rownames - not done 
rownames(logs_P01) <- logs_P01$Human_name
logs_P01 <- logs_P01[,-1] 

# Make log fold changes numeric - not done, makes genes numeric, do it before heatmap when rownames are made 
logs_P01_hg[] <- lapply(logs_P01_hg, function(x) as.numeric(as.character(x)))

# Heatmaps from gprofiler KEGG and reactome results, and IPA #### 
head(gprofilerresult_2WPI_01_nf)
head(logs_P01)
logs_P01_2 <- logs_P01 %>% rownames_to_column() 
head(logs_P01_2)
colnames(logs_P01_2) <- c("Human_name", "log2FoldChange_0WPI", "log2FoldChange_2WPI", "log2FoldChange_16WPI") 

# ** TP53 regulation Reactome #### 

n_col <- gprofilerresult_2WPI_01_nf %>% filter(term.id == "REAC:R-HSA-3700989") %>% dplyr::select(overlap.size)
n_col[1,1]
n_col
x
tp53_REA <- gprofilerresult_2WPI_01_nf %>% filter(term.id == "REAC:R-HSA-3700989") %>% dplyr::select(term.id,intersection) %>% 
  separate(intersection, into = c(stringi::stri_rand_strings(n_col[1,1], 3)), sep = ",")

tp53_REA <- t(tp53_REA) 
tp53_REA <- as.data.frame(tp53_REA)
tp53_REA <- rownames_to_column(tp53_REA)
colnames(tp53_REA) <- c("random","Human_name")
tp53_REA <- tp53_REA[-1,]
tp53_REA
tp53_REA_logs <- left_join(tp53_REA, logs_P01_2, by = "Human_name")
head(tp53_REA_logs) 
tp53_REA_logs2 <- tp53_REA_logs %>% dplyr::select(-random) 

# Make a file with logs and gene descriptions 
library(AnnotationDbi)
library(OrganismDbi)
library(org.Hs.eg.db)

tp53_REA_logs$gene_name <- mapIds(org.Hs.eg.db,
                                      keys=tp53_REA_logs$Human_name,
                                      column="ALIAS",
                                      keytype="ENSEMBL",
                                      multiVals="first")
# Alias is different from that on the heatmap - better to merge with biotype.... 
tp53_REA_logs$gene_description <- mapIds(org.Hs.eg.db,
                                             keys=tp53_REA_logs$Human_name,
                                             column="GENENAME",
                                             keytype="ENSEMBL",
                                             multiVals="first")
tp53_REA_logs <- apply(tp53_REA_logs,2,as.character) 
tp53_REA_logs <- as.data.frame(tp53_REA_logs)
write.csv(tp53_REA_logs, "REA_tp53_logs_ann.csv")

# Continue to heatmap 

View(tp53_REA_logs2)

head(tp53_REA_logs2)
tp53_REA_logs2 <- tp53_REA_logs2 %>% arrange(rownames(tp53_REA_logs2))

typeof(tp53_REA_logs2)

tp53_REA_logs2 <- left_join(tp53_REA_logs2,biotype3, by = "Human_name")
tp53_REA_logs2 <- tp53_REA_logs2 %>% distinct(gene_name, .keep_all = TRUE)
?distinct
rownames(tp53_REA_logs2) <- tp53_REA_logs2$gene_name
nrow(tp53_REA_logs2) #58
head(tp53_REA_logs2) 
tp53_REA_logs2 <- tp53_REA_logs2 %>% dplyr::select(-Human_name,-strand, -gene_biotype, -gene_name)


library(pheatmap)

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red2"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(tp53_REA_logs2$log2FoldChange_2WPI), 0, length.out=ceiling(paletteLength/2 + 1)), 
              seq((max(tp53_REA_logs2$log2FoldChange_2WPI)/paletteLength), max(tp53_REA_logs2$log2FoldChange_2WPI), length.out=floor(paletteLength/2)))
# make breaks the same for all heatmaps to compare expression

REA_tp53_heatmap <- pheatmap(tp53_REA_logs2, na_col = "grey", 
                                 cluster_cols = FALSE, cluster_rows = FALSE, 
                                 cellwidth = 15, fontsize = 8, angle_col = 315,
                             color=myColor, breaks = myBreaks) 

REA_tp53_heatmap

tiff("tp53_REA_heatmap.tiff", height = 20, width = 24, units = 'cm', 
     compression = "lzw", res = 600) 
REA_tp53_heatmap
dev.off() 

library(Cairo)
Cairo(file="tp53_REA_heatmap.png", 
      type="png",
      units="cm", 
      width=24, 
      height=20, 
      pointsize=60, 
      dpi=600)
REA_tp53_heatmap
dev.off() 

# ** Toll all ####
# Reactome 
# MAP2K3,CASP8,MEF2A,MEF2C,PTPN4,MAP3K1,CREB1,DUSP4,TRAF3,CNPY3,PTPN11,MAPK11,TLR9
# KEGG 
# MAP2K3,CASP8,CD80,TRAF3,RAC1,IFNAR2,MAPK11,TLR9
# IPA 
# MAP2K3,MAP3K1,MAPK11,PPARA,TLR9
# 

Toll_gene_list <- c("CASP8","CD80","CNPY3","CREB1","DUSP4","IFNAR2","MAP2K3","MAP3K1","MAPK11","MEF2A2","MEF2C","PPARA","PTPN4","PTPN11","RAC1","TLR9","TRAF3")
length(Toll_gene_list) #17 

Toll_all <- Toll_gene_list 

Toll_all <- t(Toll_all) 
Toll_all <- as.data.frame(Toll_all)
Toll_all <- t(Toll_all)
colnames(Toll_all) <- c("gene_name")
Toll_all <- as.data.frame(Toll_all)
Toll_all

Toll_all_logs <- left_join(Toll_all, logs_P01_hg, by = "gene_name")
head(Toll_all_logs) 
Toll_all_logs2 <- Toll_all_logs 

# Make a file with logs and gene descriptions - best to merge with biotype! 
library(AnnotationDbi)
library(OrganismDbi)
library(org.Hs.eg.db)

Toll_all_logs$Ensembl_name <- mapIds(org.Hs.eg.db,
                                     keys=Toll_all_logs$gene_name,
                                     column="ENSEMBL",
                                     keytype="SYMBOL",
                                     multiVals="first")
# Alias is different from that on the heatmap - better to merge with biotype.... 
Toll_all_logs$gene_description <- mapIds(org.Hs.eg.db,
                                         keys=Toll_all_logs$gene_name,
                                         column="GENENAME",
                                         keytype="SYMBOL",
                                         multiVals="first")
Toll_all_logs <- apply(Toll_all_logs,2,as.character) 
Toll_all_logs <- as.data.frame(Toll_all_logs)
write.csv(Toll_all_logs, "REA_Toll_logs_ann.csv")
View(Toll_all_logs)

# Continue to heatmap 

View(Toll_all_logs2)

head(Toll_all_logs2)

typeof(Toll_all_logs2)

# Create rownames - not done 
rownames(Toll_all_logs2) <- Toll_all_logs2$gene_name
Toll_all_logs2 <- Toll_all_logs2[,-1] 

# Make log fold changes numeric - not done, makes genes numeric, do it before heatmap when rownames are made 
Toll_all_logs2[] <- lapply(Toll_all_logs2, function(x) as.numeric(as.character(x)))

library(pheatmap)

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red2"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(Toll_all_logs2$log2FoldChange_2WPI), 0, length.out=ceiling(paletteLength/2 + 1)), 
              seq((max(Toll_all_logs2$log2FoldChange_2WPI)/paletteLength), max(Toll_all_logs2$log2FoldChange_2WPI), length.out=floor(paletteLength/2)))


Toll_all_heatmap <- pheatmap(Toll_all_logs2, na_col = "grey", 
                             cluster_cols = FALSE, cluster_rows = FALSE, 
                             cellwidth = 15, fontsize = 8, angle_col = 315,
                             # breaks = seq(-range, range, length.out = 50))
                             color=myColor, breaks=myBreaks) 

tiff("Toll_all_heatmap.tiff", height = 20, width = 24, units = 'cm', 
     compression = "lzw", res = 600) 
Toll_all_heatmap
dev.off() 

library(Cairo)
Cairo(file="Toll_all_heatmap.png", 
      type="png",
      units="cm", 
      width=24, 
      height=20, 
      pointsize=60, 
      dpi=600)
Toll_all_heatmap
dev.off() 

# ** RIG all #### 

# KEGG 
# CASP8,MAP3K1,DDX58,DDX3X
# IPA 
# CASP8,DDX58,NFKBID,TRAF3
# Add TRAF3 and MAPK11 by hand 

RIG_gene_list <- c("CASP8", "DDX58","DDX3X","MAP3K1","MAPK11","NFKBID","TRAF3")
length(RIG_gene_list) #7

RIG_all <- RIG_gene_list 

RIG_all <- t(RIG_all) 
RIG_all <- as.data.frame(RIG_all)
RIG_all <- t(RIG_all)
colnames(RIG_all) <- c("gene_name")
RIG_all <- as.data.frame(RIG_all)
RIG_all

RIG_all_logs <- left_join(RIG_all, logs_P01_hg, by = "gene_name")
head(RIG_all_logs) 
RIG_all_logs2 <- RIG_all_logs 

# Make a file with logs and gene descriptions - best to merge with biotype! 
library(AnnotationDbi)
library(OrganismDbi)
library(org.Hs.eg.db)

RIG_all_logs$Ensembl_name <- mapIds(org.Hs.eg.db,
                                    keys=RIG_all_logs$gene_name,
                                    column="ENSEMBL",
                                    keytype="SYMBOL",
                                    multiVals="first")
# Alias is different from that on the heatmap - better to merge with biotype.... 
RIG_all_logs$gene_description <- mapIds(org.Hs.eg.db,
                                        keys=RIG_all_logs$gene_name,
                                        column="GENENAME",
                                        keytype="SYMBOL",
                                        multiVals="first")
RIG_all_logs <- apply(RIG_all_logs,2,as.character) 
RIG_all_logs <- as.data.frame(RIG_all_logs)
write.csv(RIG_all_logs, "RIG_all_logs_ann.csv")
View(RIG_all_logs)

# Continue to heatmap 

View(RIG_all_logs2)

head(RIG_all_logs2)

typeof(RIG_all_logs2)

# Create rownames - not done 
rownames(RIG_all_logs2) <- RIG_all_logs2$gene_name
RIG_all_logs2 <- RIG_all_logs2[,-1] 

# Make log fold changes numeric - not done, makes genes numeric, do it before heatmap when rownames are made 
RIG_all_logs2[] <- lapply(RIG_all_logs2, function(x) as.numeric(as.character(x)))

library(pheatmap)

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red2"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(RIG_all_logs2$log2FoldChange_2WPI), 0, length.out=ceiling(paletteLength/2 + 1)), 
              seq((max(RIG_all_logs2$log2FoldChange_2WPI)/paletteLength), max(RIG_all_logs2$log2FoldChange_2WPI), length.out=floor(paletteLength/2)))


RIG_all_heatmap <- pheatmap(RIG_all_logs2, na_col = "grey", 
                            cluster_cols = FALSE, cluster_rows = FALSE, 
                            cellwidth = 15, fontsize = 8, angle_col = 315,
                            # breaks = seq(-range, range, length.out = 50))
                            color=myColor, breaks=myBreaks) 

tiff("RIG_all_heatmap.tiff", height = 10, width = 24, units = 'cm', 
     compression = "lzw", res = 600) 
RIG_all_heatmap
dev.off() 

library(Cairo)
Cairo(file="RIG_all_heatmap.png", 
      type="png",
      units="cm", 
      width=24, 
      height=10, 
      pointsize=60, 
      dpi=600)
RIG_all_heatmap
dev.off() 

# ** PPAR IPA #### 
PPAR_gene_list <- c("CASP8","CASP9","DDX58","HSPA5","IFNAR2","IRF9","MAP2K3","MAPK11","NFKBID","TLR9","TRAF3")
length(PPAR_gene_list) #11 

PPAR_IPA <- PPAR_gene_list 

PPAR_IPA <- t(PPAR_IPA) 
PPAR_IPA <- as.data.frame(PPAR_IPA)
PPAR_IPA <- t(PPAR_IPA)
colnames(PPAR_IPA) <- c("gene_name")
PPAR_IPA <- as.data.frame(PPAR_IPA)
PPAR_IPA

PPAR_IPA_logs <- left_join(PPAR_IPA, logs_P01_hg, by = "gene_name")
head(PPAR_IPA_logs) 
PPAR_IPA_logs2 <- PPAR_IPA_logs 

# Make a file with logs and gene descriptions - best to merge with biotype! 
library(AnnotationDbi)
library(OrganismDbi)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

PPAR_IPA_logs$Ensembl_name <- mapIds(org.Hs.eg.db,
                                    keys=PPAR_IPA_logs$gene_name,
                                    column="ENSEMBL",
                                    keytype="SYMBOL",
                                    multiVals="first")
# Alias is different from that on the heatmap - better to merge with biotype.... 
PPAR_IPA_logs$gene_description <- mapIds(org.Hs.eg.db,
                                        keys=PPAR_IPA_logs$gene_name,
                                        column="GENENAME",
                                        keytype="SYMBOL",
                                        multiVals="first")
PPAR_IPA_logs <- apply(PPAR_IPA_logs,2,as.character) 
PPAR_IPA_logs <- as.data.frame(PPAR_IPA_logs)
write.csv(PPAR_IPA_logs, "PPAR_IPA_logs_ann.csv")
View(PPAR_IPA_logs)

# Continue to heatmap 

View(PPAR_IPA_logs2)

head(PPAR_IPA_logs2)

typeof(PPAR_IPA_logs2)

# Create rownames - not done 
rownames(PPAR_IPA_logs2) <- PPAR_IPA_logs2$gene_name
PPAR_IPA_logs2 <- PPAR_IPA_logs2[,-1] 

# Make log fold changes numeric - not done, makes genes numeric, do it before heatmap when rownames are made 
PPAR_IPA_logs2[] <- lapply(PPAR_IPA_logs2, function(x) as.numeric(as.character(x)))

library(pheatmap)

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red2"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(PPAR_IPA_logs2$log2FoldChange_2WPI), 0, length.out=ceiling(paletteLength/2 + 1)), 
              seq((max(PPAR_IPA_logs2$log2FoldChange_2WPI)/paletteLength), max(PPAR_IPA_logs2$log2FoldChange_2WPI), length.out=floor(paletteLength/2)))


PPAR_IPA_heatmap <- pheatmap(PPAR_IPA_logs2, na_col = "grey", 
                            cluster_cols = FALSE, cluster_rows = FALSE, 
                            cellwidth = 15, fontsize = 8, angle_col = 315,
                            # breaks = seq(-range, range, length.out = 50))
                            color=myColor, breaks=myBreaks) 

tiff("PPAR_IPA_heatmap.tiff", height = 10, width = 24, units = 'cm', 
     compression = "lzw", res = 600) 
PPAR_IPA_heatmap
dev.off() 

library(Cairo)
Cairo(file="PPAR_IPA_heatmap.png", 
      type="png",
      units="cm", 
      width=24, 
      height=10, 
      pointsize=60, 
      dpi=600)
PPAR_IPA_heatmap
dev.off() 

#### Combine TP53, RIG, Toll and PPAR #### 

library(ggpubr)
# 4 columns - not used 
h <- ggarrange(REA_tp53_heatmap[[4]], Toll_all_heatmap[[4]],RIG_all_heatmap[[4]],PPAR_IPA_heatmap[[4]],ncol=4, labels=c("A", "B", "C","D"))
h

# This for paper 
j <- ggarrange(REA_tp53_heatmap[[4]], Toll_all_heatmap[[4]],                                               # First column 
               ggarrange(RIG_all_heatmap[[4]], PPAR_IPA_heatmap[[4]], nrow = 2, labels = c("C", "D"),heights = c(1,1.3)), # Second column  
               ncol = 3, 
               labels = c("A","B")                                   # Labels of the first column
) 

j

library(Cairo)
Cairo(file="P53_Toll_Rig_PKR_ABCD.png", 
      type="png",
      units="cm", 
      width=24, 
      height=24, 
      pointsize=60, 
      dpi=600)
j
dev.off() 

tiff(filename = "Fig5_P53_Toll_Rig_PKR_ABCD.tiff", compression = "lzw", width = 24, height = 24, 
     units = "cm", res = 600);
j;
dev.off(); 