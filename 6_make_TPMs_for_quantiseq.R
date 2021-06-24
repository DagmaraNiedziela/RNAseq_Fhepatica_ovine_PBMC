# RPKM versus TPM
# 
# RPKM and TPM are both normalized for library size and gene length.
#
# RPKM is not comparable across different samples.
#
# For more details, see: http://blog.nextgenetics.net/?e=51

head(newdata)# this is my counts matrix, with protein coding genes only 
nrow(newdata) #20477 
head(gtf_unique_genes) # This is from gtf file, with all start and ends of genes 
# width is gene length 
TPM_gene_info <- gtf_unique_genes %>% filter(gene_biotype == "protein_coding") %>% select(gene_id, width, gene_name)
head(TPM_gene_info)
genes <- TPM_gene_info %>% select(gene_id, width) 
colnames(genes) <- c("Gene", "Length")
head(genes)

# Function setups 
rpkm <- function(counts, lengths) {
  rate <- counts / lengths 
  rate / sum(counts) * 1e6
}

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

genes <- data.frame(
  Gene = c("A","B","C","D","E"),
  Length = c(100, 50, 25, 5, 1)
)

counts <- data.frame(
  S1 = c(80, 10,  6,  3,   1),
  S2 = c(20, 20, 10, 50, 400)
)

rpkms <- apply(counts, 2, function(x) rpkm(x, genes$Length))

# Make my tpms 
tpms <- apply(newdata, 2, function(x) tpm(x, genes$Length))
head(tpms)

# Input and output tables 
genes
#   Gene Length
# 1    A    100
# 2    B     50
# 3    C     25
# 4    D      5
# 5    E      1

counts
#   S1  S2
# 1 80  20
# 2 10  20
# 3  6  10
# 4  3  50
# 5  1 400

rpkms
#         S1    S2
# [1,]  8000 4e+02
# [2,]  2000 8e+02
# [3,]  2400 8e+02
# [4,]  6000 2e+04
# [5,] 10000 8e+05

tpms
#             S1         S2
# [1,] 281690.14    486.618
# [2,]  70422.54    973.236
# [3,]  84507.04    973.236
# [4,] 211267.61  24330.900
# [5,] 352112.68 973236.010

# Sample means should be equal.

colSums(rpkms)
#    S1     S2 
# 28400 822000
colSums(tpms)
#    S1    S2 
# 1e+06 1e+06

colMeans(rpkms)
#   S1     S2 
# 5680 164400
colMeans(tpms)
#    S1    S2 
# 2e+05 2e+05

# colsums and colmeans correct 

# Change the TPM matrix to one with gene symbols instead of Ensembl IDs 
wow <- merge(tpms, TPM_gene_info, by.x=0, by.y=1, all=FALSE)
head(wow)
View(wow)
library(tidyverse)
# Get rid of NA names which are way duplicated 
wow <- wow %>% filter(!gene_name == "<NA>") 
# Look at duplicates 
wow[duplicated(wow$gene_name), ] 
# 1945   0.00000000   0.00000000   0.00000000 0.000000e+00   0.04287995 0.000000e+00   426     RPS25 2.839765e-01
# 2304   0.03082046   0.02111936   0.04724112 4.259592e-02   0.01712301 4.885863e-02   426     RPS25 1.831740e+00 
# Example of a duplicated gene 
wow %>% filter(gene_name == "RPS25") 

# Create a new column with sums of rows of columns 2-49
wow <-  wow %>%  mutate(summed = rowSums(wow[, c(2:49)]))

# Remove duplicate gene names by picking the highest row sum per gene name 
tpm_for_save <- wow 
colnames(tpm_for_save)
head(tpm_for_save)

tpm_for_save <- tpm_for_save %>%
  arrange(gene_name, -summed) %>%
  filter(duplicated(gene_name) == FALSE) 

tpm_for_save %>% filter(gene_name == "RPS25") # just as a check before and after removing duplicates 

rownames(tpm_for_save) <-  tpm_for_save$gene_name 
nrow(tpm_for_save) 

tpm_for_save2 <- tpm_for_save[,c(2:49)]
head(tpm_for_save2)

# Save as a tab delimited file 
write.table(tpm_for_save2, file = "tpms_for_quantiseq.tab", sep = "\t", quote = FALSE)
