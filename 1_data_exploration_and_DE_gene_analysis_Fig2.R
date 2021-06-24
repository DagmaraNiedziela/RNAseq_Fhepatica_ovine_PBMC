# RNASeq DE gene analysis of ovine PBMC in response to Fasciola hepatica - file 1 ################
# Dagmara Niedziela 
# This file includes data exploration and DE gene analysis 1 - infected vs control 
library("DESeq2")
sessionInfo() # DESeq version 
# 1.18.1 

#Load matrix  ####
coldata <- read.csv("sample_metadata.csv")
countdata <- read.csv("gene_counts.csv", row.names = "Genename") #row.names clause essential to get equal numbers of rows vs columns in coldata and countdata
nrow(coldata) #48 check if the two are equal - otherwise the dds matrix won't load 
ncol(countdata) #48
head(countdata)

head(coldata)
coldata$WPI <- factor(coldata$WPI)
levels(coldata$WPI) 

coldata$condition <- paste(coldata$Group, coldata$WPI, sep = "_","WPI")  

max(countdata)
# 540457 
summary(countdata)
# Most counts really low? And mostly evenly distributed across samples 
#mean(countdata, na.rm = TRUE)

# Load my dds matrix with both strains - this is used for data exploration 
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata, 
                              design = ~ WPI + Group + WPI:Group) 

rowData(dds) <- anno 							  
head(dds) 
nrow(dds) 
#[1] 26478 

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
nrow(dds) 
#[1] 16721 

rowData(dds)
colnames(dds) 

# Get the gtf file and transform so that I can have rowData for DeSeq2 #### 
# This will all be used for annotations later 

# BiocManager::install("rtracklayer")
library(rtracklayer)
library(tidyverse) 

#Importing gtf and transform to data frame
rtracklayer::import("Ovis_aries_rambouillet.Oar_rambouillet_v1.0.101.gtf",
                    format = "gtf") %>%
  as.data.frame() -> gtf

#Check data frame
View(gtf) 
colnames(gtf)
nrow(gtf) # 991971

gtf_unique_genes <- gtf %>% distinct(gene_id, .keep_all = TRUE) 
nrow(gtf_unique_genes) # 26478
View(gtf_unique_genes)
head(gtf_unique_genes) 
27552-26098
140365-136163 
# width is gene length it seems 

gtf_gene_info <- gtf_unique_genes %>% select(gene_id, strand, gene_biotype, gene_name)
write.csv(gtf_gene_info, "gtf_gene_info.csv")

###### NORMALISE AND CLUSTER PLOTS ####

# ** Normalise data with a vst function ####
vsd <- vst(dds, blind=FALSE)

# ** sample distances heat map ####
sampleDists <- dist(t(assay(vsd)))
library("pheatmap")
library("RColorBrewer") 
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$Sheep_ID, vsd$WPI, vsd$Group, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap <- pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors) 

tiff("sample_distances_heatmap.tiff", height = 20, width = 25, units = 'cm', 
     compression = "lzw", res = 600) 
heatmap
dev.off() 

#PCA of samples 
plotPCA(vsd, intgroup=c("Group", "WPI")) # The simplest PCA 

PCAdata <- plotPCA(vsd, intgroup=c("Group", "WPI", "Sheep_ID"), returnData = TRUE) 

# MY PCA PLOT! Strain is shape, HPI is color ####
percentVar <- round(100 * attr(PCAdata, "percentVar")) 
library("ggplot2")
sp <- ggplot(PCAdata, aes(x = PC1, y = PC2, color = WPI, shape = Group)) +
  geom_point(size =3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) #+ 
  #coord_fixed() 
sp 

# ** My PCA formatted ####
myPCA <- sp + scale_color_brewer(palette="Set1", direction = -1) + 
  scale_fill_brewer(palette = "Set1", direction = -1) + 
  theme(text = element_text(size = 14, family = "Calibri")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
myPCA

ggsave("myPCA2.jpeg") 

#** plot PCA with each time point as separate grid ####
# and unfix the scales in each square 
PCAplot3 <- ggplot(data = PCAdata) + 
  geom_point(mapping = aes(x = PC1, y = PC2, color = Group)) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  facet_wrap(~ WPI, nrow = 1, scales = "free") +
  theme(text = element_text(size = 14, family = "Calibri")) 

myPCA_faceted <- PCAplot3 + scale_color_brewer(palette="Set1", direction = -1) 
myPCA_faceted
ggsave("PCA_plot_faceted_uncoord.jpeg")

##########################################
# ** FOR paper - PCA and faceted PCA combined - Figure 2
library(ggpubr)

plotPCA <- ggarrange(myPCA, myPCA_faceted, nrow =2, labels = c("A","B"),heights = c(1.4,1))
plotPCA 
ggsave("PCA_combined.jpeg", height = 9, width = 8, units = "in", dpi = 1200) 

# Data exploration DONE #####

# !!!Load matrix - protein coding genes only and grouping for individual animals #### 

# Remove non protein coding genes - with my gtf import 
biotype <- as.matrix(read.csv("gtf_gene_info.csv", row.names="gene_id")) # file attached, import it  
te <- merge(countdata, biotype, by=0, all=FALSE) #merge my counts file with the protein coding assignment file 
head(te) 
rownames(te) <- te[,1] # this is making rownames  
te<-te[,-1] #removing column of ENSEML that’s in twice 
ncol(te)
nrow(te) # 26478
newdata <- te[ which(te$gene_biotype=='protein_coding'), ] #remove genes that are not protein coding (miRNAs, pseudogenes, long ncRNAs etc)
head(newdata)
nrow(newdata) #20477 
ncol(newdata) #52
newdata <- newdata[,-c(49:52)] # remove all the gene descriptors 
# newdata contains protein coding genes only 
# countdata <- newdata #make newdata my countdata 

coldata_group <- coldata
head(coldata_group)
coldata_group$ind.n <- factor(rep(rep(1:8,each=1),6))
View(coldata_group)
# model.matrix(~ grp + grp:ind.n + grp:cnd, coldata) 

library("DESeq2") 

##### RUN DE GENE ANALYSIS - GROUP3 - group and protein coding genes only #####
# Please note - PCA has been checked for this dds object also and looks THE SAME 

dds_group3 <- DESeqDataSetFromMatrix(countData = newdata,
                                     colData = coldata_group, 
                                     design = ~ WPI + WPI:ind.n + WPI:Group) 
# some variables in design formula are characters, converting to factors 

head(dds_group3)
colnames(dds_group3) 
nrow(dds_group3) 
#[1] 20477

keep <- rowSums(counts(dds_group3)) >= 10
dds_group3 <- dds_group3[keep,]
nrow(dds_group3) 
#[1] 15363 

dds_group3 <- estimateSizeFactors(dds_group3) 

nrow(dds_group3) #15363 
dds_group3 <- DESeq(dds_group3) #runs analysis 
# 3 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest 
# Solution was found here https://support.bioconductor.org/p/65091/ 
dds_group3[which(mcols(dds_group3)$betaConv),] #rows 15360
dds_group3 <- dds_group3[which(mcols(dds_group3)$betaConv),]
dds_group3
res_group3 <- results(dds_group3) 
res_group3
#Wald test p-value: WPI16.GroupInfected 

resultsNames(dds_group3) 
#[1] "Intercept"           "WPI_2_vs_0"          "WPI_16_vs_0"         "WPI0.ind.n2"         "WPI2.ind.n2"        
#[6] "WPI16.ind.n2"        "WPI0.ind.n3"         "WPI2.ind.n3"         "WPI16.ind.n3"        "WPI0.ind.n4"        
#[11] "WPI2.ind.n4"         "WPI16.ind.n4"        "WPI0.ind.n5"         "WPI2.ind.n5"         "WPI16.ind.n5"       
#[16] "WPI0.ind.n6"         "WPI2.ind.n6"         "WPI16.ind.n6"        "WPI0.ind.n7"         "WPI2.ind.n7"        
#[21] "WPI16.ind.n7"        "WPI0.ind.n8"         "WPI2.ind.n8"         "WPI16.ind.n8"        "WPI0.GroupInfected" 
#[26] "WPI2.GroupInfected"  "WPI16.GroupInfected"

# Results tables per time point 
res_0WPI <- results(dds_group3, name="WPI0.GroupInfected")
res_2WPI <- results(dds_group3, name="WPI2.GroupInfected") 
res_16WPI <- results(dds_group3, name="WPI16.GroupInfected") 

# ** Results shrinkage #####
BiocManager::install("apeglm")
library(apeglm)
?lfcShrink
resLFC_0WPI <- lfcShrink(dds_group3, coef="WPI0.GroupInfected", type="apeglm")
resLFC_0WPI
resLFC_2WPI <- lfcShrink(dds_group3, coef="WPI2.GroupInfected", type="apeglm")
resLFC_16WPI <- lfcShrink(dds_group3, coef="WPI16.GroupInfected", type="apeglm")

mcols(res_group3, use.names = TRUE) #meaning of columns 
summary(res_group3) #gives info about outliers and low counts 

sum(!is.na(res_group3$pvalue)) #genes with a reported P value 
## 15360 

sum(res_0WPI$padj < 0.1, na.rm=TRUE) # 59
sum(res_2WPI$padj < 0.1, na.rm=TRUE) # 453
sum(res_16WPI$padj < 0.1, na.rm=TRUE) # 2 
# And alpha is set to 0.1 

sum(res_0WPI$padj < 0.05, na.rm=TRUE) # 24
sum(res_2WPI$padj < 0.05, na.rm=TRUE) # 229
sum(res_16WPI$padj < 0.05, na.rm=TRUE) # 2 

##### ANNOTATE GENES ##### 
# Inner join with gtf import
head(biotype) # use this, has row names 
head(gtf_gene_info)
head(res_0WPI)

biotype2 <- biotype[,-1]
head(biotype2)

# Use LFC shrunk values 
nrow(res_0WPI) #15360 
res_0WPI_ann <- merge(resLFC_0WPI, biotype2, by=0, all=FALSE)
head(res_0WPI_ann) # row name becomes a column, careful that might change things later
nrow(res_0WPI_ann) #15360 
res_2WPI_ann <- merge(resLFC_2WPI, biotype2, by=0, all=FALSE)
res_16WPI_ann <- merge(resLFC_16WPI, biotype2, by=0, all=FALSE)

# ** P < 0.1 Subset results for export #### 

res_0WPI_sig01 <- subset(resLFC_0WPI, padj < 0.1)
res_2WPI_sig01 <- subset(resLFC_2WPI, padj < 0.1)
res_16WPI_sig01 <- subset(resLFC_16WPI, padj < 0.1)

res_0WPI_sig01_DF <- as.data.frame(subset(res_0WPI_ann, padj < 0.1)) 
res_2WPI_sig01_DF <- as.data.frame(subset(res_2WPI_ann, padj < 0.1)) 
res_16WPI_sig01_DF <- as.data.frame(subset(res_16WPI_ann, padj < 0.1)) 
head(res_0WPI_sig01_DF)

######EXPORT TABLES #####
#all genes that are significant 
write.csv(res_0WPI_sigDF, file = "DE_genes_0WPI.csv") 
write.csv(res_2WPI_sigDF, file = "DE_genes_2WPI.csv") 
write.csv(res_16WPI_sigDF, file = "DE_genes_16WPI.csv") 

# P value < 0.1 
write.csv(res_0WPI_sig01_DF, file = "DE_genes_0WPI_P01.csv") 
write.csv(res_2WPI_sig01_DF, file = "DE_genes_2WPI_P01.csv") 
write.csv(res_16WPI_sig01_DF, file = "DE_genes_16WPI_P01.csv") 

#export tables with all genes, not just P < 0.05 
res_0WPI_DF <- as.data.frame(res_0WPI_ann)
write.csv(res_0WPI_DF, file = "All_results_0WPI.csv") 

res_2WPI_DF <- as.data.frame(res_2WPI_ann)
write.csv(res_2WPI_DF, file = "All_results_2WPI.csv")  

res_16WPI_DF <- as.data.frame(res_16WPI_ann)
write.csv(res_16WPI_DF, file = "All_results_16WPI.csv")  

# No LFC threshold datasets because I have very small log2FCs 

# DE gene numbers (Figure 3) are taken in "2_DE_gene_amounts_and_DE_gene_plot.R" ##### 
# Top 5 genes (Figure 4) are in "3_top5_genes_and_milk_genes.r" 
