# RNASeq DE gene analysis of ovine PBMC in response to Fasciola hepatica - file 1 ################
# Dagmara Niedziela 
# This file includes data exploration and DE gene analysis 2 - longitudinal 

#Load matrix ####

# Load with newdata - protein coding genes only, data exploration is always the same - as in file 1 (infected vs control)
head(coldata)

library("DESeq2")
dds4 <- DESeqDataSetFromMatrix(countData = newdata,
                              colData = coldata, 
                              design = ~ condition) 

nrow(dds4) 
#[1] 20477

keep <- rowSums(counts(dds4)) >= 10
dds4 <- dds4[keep,]
nrow(dds4) 
#[1] 15363 

head(dds4)
colnames(dds4) 

dds4 <- estimateSizeFactors(dds4) 

# DE gene analysis ####

dds4 <- DESeq(dds4) #runs analysis 
# -- replacing outliers and refitting for 34 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds) 

dds4
res4 <- results(dds4) 
res4
#Wald test p-value: condition Infected 2 WPI vs Control 0 WPI 

resultsNames(dds4) 
#[1] "Intercept"                                  "condition_Control_16_WPI_vs_Control_0_WPI" 
#[3] "condition_Control_2_WPI_vs_Control_0_WPI"   "condition_Infected_0_WPI_vs_Control_0_WPI" 
#[5] "condition_Infected_16_WPI_vs_Control_0_WPI" "condition_Infected_2_WPI_vs_Control_0_WPI" 

coldata$condition <- factor(coldata$condition)
levels(coldata$condition) # [1] "Control_0_WPI"   "Control_16_WPI"  "Control_2_WPI"   "Infected_0_WPI"  "Infected_16_WPI"
# [6] "Infected_2_WPI"

# Results tables per time point 
res4_0WPI <- results(dds4, contrast = c("condition", "Infected_0_WPI", "Control_0_WPI"))
sum(res4_0WPI$padj < 0.1, na.rm=TRUE) # 14 

head(subset(res4_0WPI, padj < 0.05))
#ENSOARG00020002731   novel gene 
#ENSOARG00020006258  glycoprotein Ib platelet subunit beta, GP1BB
#ENSOARG00020007641 diacylglycerol kinase delta, DGKD 
#ENSOARG00020008718 WRN helicase interacting protein 1, WRNIP1
#ENSOARG00020010668 novel gene, liver carboxylesterase-like

res4_2WPI <- results(dds4, contrast = c("condition", "Infected_2_WPI", "Control_2_WPI"))
sum(res4_2WPI$padj < 0.1, na.rm=TRUE) # 164  

res4_16WPI <- results(dds4, contrast = c("condition", "Infected_16_WPI", "Control_16_WPI"))
sum(res4_16WPI$padj < 0.1, na.rm=TRUE) # 8 
# The above are just checks to compare to the infected vs control analysis - ie what happens when you dont include animal as batch effect - we see less DE genes 

# DE genes vs time 0 in the Infected group 
res4_inf_2WPI <- results(dds4, contrast = c("condition", "Infected_2_WPI", "Infected_0_WPI"))
sum(res4_inf_2WPI$padj < 0.05, na.rm=TRUE) # 1522 
nrow(res4_inf_2WPI) #15363 

res4_inf_16WPI <- results(dds4, contrast = c("condition", "Infected_16_WPI", "Infected_0_WPI")) 
sum(res4_inf_16WPI$padj < 0.05, na.rm=TRUE) # 3949

# DE genes vs time 0 in the Control group 
res4_cont_2WPI <- results(dds4, contrast = c("condition", "Control_2_WPI", "Control_0_WPI"))
sum(res4_cont_2WPI$padj < 0.05, na.rm=TRUE) # 769
res4_cont_16WPI <- results(dds4, contrast = c("condition", "Control_16_WPI", "Control_0_WPI"))
sum(res4_cont_16WPI$padj < 0.05, na.rm=TRUE) # 2342 

# Take the infected vs T0 and control vs T0, shrink them, then intersect the 2WPI and 16 WPI lists/data frames 

# LFC shrinkage #### 
# BiocManager::install("ashr")
# res4_LFC_2WPI <- lfcShrink(dds4, contrast = c("condition", "Infected_2_WPI", "Control_2_WPI"), type="ashr") # This is just in case I want to use the non-animal batch analysis 

res4_LFC_inf_2WPI <- lfcShrink(dds4, contrast = c("condition", "Infected_2_WPI", "Infected_0_WPI"), type="ashr")
res4_LFC_inf_16WPI <- lfcShrink(dds4, contrast = c("condition", "Infected_16_WPI", "Infected_0_WPI"), type="ashr")
res4_LFC_cont_2WPI <- lfcShrink(dds4, contrast = c("condition", "Control_2_WPI", "Control_0_WPI"), type="ashr")
res4_LFC_cont_16WPI <- lfcShrink(dds4, contrast = c("condition", "Control_16_WPI", "Control_0_WPI"), type="ashr")

# Error in lfcShrink(dds4, contrast = c("condition", "Infected_2_WPI", "Control_2_WPI"),  : 
#                      type='apeglm' shrinkage only for use with 'coef' 
# Other types will work with contrast - I will use ashr 
# Complex models are suggested to make all the coefs - like dds_group3
# Or relevelling and rerunning the model 

# Check this link: https://support.bioconductor.org/p/98833/ 
# To find out how the change in DESeq2 now calculates p values for unshrunken matrix and 
# then shrinks log2foldchanges with lfcShrink 
# LFCs for contrast will be different to LFCs for coef 
# coef uses one column of data matrix (data versus reference), contrast compares two 
# Should you use shrinkage 
# https://www.biostars.org/p/340269/ 
# I think shrinkage should be used, because low count genes can be prone to sequencing error 

# Annotate LFC shrunk lists #### 

# Inner join with gtf import
head(biotype) # use this, has row names 
head(gtf_gene_info)
head(res4_LFC_inf_2WPI)

head(biotype2)

# Use LFC shrunk values 
res4_LFC_inf_2WPI_ann <- merge(res4_LFC_inf_2WPI, biotype2, by=0, all=FALSE)
head(res4_LFC_inf_2WPI_ann) # row name becomes a column, careful that might change things later
res4_LFC_inf_16WPI_ann <- merge(res4_LFC_inf_16WPI, biotype2, by=0, all=FALSE)
res4_LFC_cont_2WPI_ann <- merge(res4_LFC_cont_2WPI, biotype2, by=0, all=FALSE)
res4_LFC_cont_16WPI_ann <- merge(res4_LFC_cont_16WPI, biotype2, by=0, all=FALSE)

# Subset results for export ####
res4_LFC_inf_2WPI_sig <- subset(res4_LFC_inf_2WPI, padj < 0.05)
res4_LFC_inf_16WPI_sig <- subset(res4_LFC_inf_16WPI, padj < 0.05)
res4_LFC_cont_2WPI_sig <- subset(res4_LFC_cont_2WPI, padj < 0.05)
res4_LFC_cont_16WPI_sig <- subset(res4_LFC_cont_16WPI, padj < 0.05)

res4_LFC_inf_2WPI_sigDF <- as.data.frame(subset(res4_LFC_inf_2WPI_ann, padj < 0.05))
res4_LFC_inf_16WPI_sigDF <- as.data.frame(subset(res4_LFC_inf_16WPI_ann, padj < 0.05)) 
res4_LFC_cont_2WPI_sigDF <- as.data.frame(subset(res4_LFC_cont_2WPI_ann, padj < 0.05)) 
res4_LFC_cont_16WPI_sigDF <- as.data.frame(subset(res4_LFC_cont_16WPI_ann, padj < 0.05)) 
##head(res_group3_24h_LF_sigDF)

# DIFFERENCE - Gene list overlaps #### 
# Use the data frames 

library(dplyr)
nrow(res4_LFC_inf_2WPI_sigDF) #1506 
nrow(res4_LFC_cont_2WPI_sigDF) #769 
res4_LFC_diff_2WPI_sigDF <- anti_join(res4_LFC_inf_2WPI_sigDF,res4_LFC_cont_2WPI_sigDF, by = "Row.names")
length(setdiff(res4_LFC_inf_2WPI_sigDF$Row.names,res4_LFC_cont_2WPI_sigDF$Row.names)) #1148 
# setdiff works on vectors, anti_join better for data frames 
nrow(res4_LFC_diff_2WPI_sigDF) # 1148 
head(res4_LFC_diff_2WPI_sigDF)

nrow(res4_LFC_inf_16WPI_sigDF) #3949 
nrow(res4_LFC_cont_16WPI_sigDF) #2342 
length(setdiff(res4_LFC_inf_16WPI_sigDF$Row.names,res4_LFC_cont_16WPI_sigDF$Row.names)) #1927  
res4_LFC_diff_16WPI_sigDF <- anti_join(res4_LFC_inf_16WPI_sigDF,res4_LFC_cont_16WPI_sigDF, by = "Row.names")
nrow(res4_LFC_diff_16WPI_sigDF) # 1927  
head(res4_LFC_diff_16WPI_sigDF)

######EXPORT TABLES #####

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

# DE gene numbers are taken in "2_DE_gene_amounts_and_DE_gene_plot.R" ##### 
# Top 5 genes and milk protein genes are in "3_top5_genes_and_milk_genes.r" 
