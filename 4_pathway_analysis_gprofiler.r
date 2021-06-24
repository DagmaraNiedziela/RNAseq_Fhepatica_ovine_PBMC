#### Pathway analysisi using gProfileR ##### 
# Dagmara Niedziela 

# Infected vs Control #####
# Data reminder 
res_0WPI_sig01_DF <- as.data.frame(subset(res_0WPI_ann, padj < 0.1)) 
res_2WPI_sig01_DF <- as.data.frame(subset(res_2WPI_ann, padj < 0.1)) 
res_16WPI_sig01_DF <- as.data.frame(subset(res_16WPI_ann, padj < 0.1)) 

#order results by adjusted p value 
res_0WPI_sig01_DF <- res_0WPI_sig01_DF[order (res_0WPI_sig01_DF$padj),] 
res_2WPI_sig01_DF <- res_2WPI_sig01_DF[order (res_2WPI_sig01_DF$padj),] 
res_16WPI_sig01_DF <- res_16WPI_sig01_DF[order (res_16WPI_sig01_DF$padj),]

head(res_2WPI_sigDF)

# Convert to human orthologs ####
# Use AnnotationDbi and add human Ensembl ID from gene symbols 
# This is not necessary as gProfileR can work with HGNC symbols as well 

library(AnnotationDbi)
library(OrganismDbi)
library(org.Hs.eg.db)

res_0WPI_sig01_DF$Human_name <- mapIds(org.Hs.eg.db,
                                    keys=res_0WPI_sig01_DF$gene_name,
                                    column="ENSEMBL",
                                    keytype="SYMBOL",
                                    multiVals="first") #column previously called Name

res_2WPI_sig01_DF$Human_name <- mapIds(org.Hs.eg.db,
                                    keys=res_2WPI_sig01_DF$gene_name,
                                    column="ENSEMBL",
                                    keytype="SYMBOL",
                                    multiVals="first") 

res_16WPI_sig01_DF$Human_name <- mapIds(org.Hs.eg.db,
                                     keys=res_16WPI_sig01_DF$gene_name,
                                     column="ENSEMBL",
                                     keytype="SYMBOL",
                                     multiVals="first") 

res_0WPI_sig01_DF <- apply(res_0WPI_sig01_DF,2,as.character)
res_2WPI_sig01_DF <- apply(res_2WPI_sig01_DF,2,as.character)
res_16WPI_sig01_DF <- apply(res_16WPI_sig01_DF,2,as.character)

res_2WPI_sig01_DF <- as.data.frame(res_2WPI_sig01_DF)
res_0WPI_sig01_DF <- as.data.frame(res_0WPI_sig01_DF)
res_16WPI_sig01_DF <- as.data.frame(res_16WPI_sig01_DF) 
head(res_2WPI_sig01_DF)

# Pathway analysis on human genes #####
library(gProfileR) 

names_2WPI_01 <- res_2WPI_sig01_DF$Human_name ## gets the gene name because that’s all we really need at this stage
names_0WPI_01 <- res_0WPI_sig01_DF$Human_name
names_16WPI_01 <- res_16WPI_sig01_DF$Human_name 

names_2WPI_01_HGNC <- res_2WPI_sig01_DF$gene_name
length(names_2WPI_01_HGNC) #453

# how many genes don't have a HGNC symbol 
sum(is.na(names_2WPI_01_HGNC)) #60
60/453*100 #13.24 

# Run pathway analysis 
head(names_2WPI_01)
typeof(names_2WPI_01)
class(names_2WPI_01)
names_2WPI_01 <- as.vector(names_2WPI_01)
?gprofiler
gprofilerresult_2WPI_01_nf <- gprofiler(names_2WPI_01, organism = "hsapiens", ordered_query = T, 
                                  significant = T, exclude_iea = F, underrep = F, evcodes = F, 
                                  region_query = F, max_p_value = 1, min_set_size = 0, 
                                  max_set_size = 0, min_isect_size = 0, 
                                  correction_method = "analytical", hier_filtering = "none", 
                                  domain_size = "annotated", custom_bg = "",numeric_ns = "", 
                                  png_fn = NULL, include_graph = F, src_filter = NULL) 

#hier_filtering = "none" ----> strong, apart from that default (ordered query true cause our genes are sorted by Padj value)
#gives output for gene ontology, kegg, reactome 

head(gprofilerresult_2WPI_01_nf)  

gprofilerresult_0WPI_01_nf <- gprofiler(names_0WPI_01, ordered_query = T, organism = "hsapiens", hier_filtering = "none")

head(gprofilerresult_0WPI_01_nf)  

write.csv(gprofilerresult_0WPI_01_nf , file="Gprofiler_0WPI_P01_nf.csv")
write.csv(gprofilerresult_2WPI_01_nf , file="Gprofiler_2WPI_P01_nf.csv")  

library(dplyr) 
gprofilerresult_2WPI_01_nf2 <- gprofilerresult_2WPI_01_nf %>% 
  filter(domain %in% c("BP", "keg", "rea")) %>% 
  arrange(domain, p.value) 

View(gprofilerresult_2WPI_01_nf2) 
write.csv(gprofilerresult_2WPI_01_nf2, "gprofilerresult_2WPI_P01_nf_edited.csv")

gprofilerresult_2WPI_01_nf3 <- gprofilerresult_2WPI_01_nf2 %>% 
  select(domain, term.name, p.value, overlap.size) 
colnames(gprofilerresult_2WPI_01_nf3) <- c("Pathway type", "Pathway", "P value", "Number of genes") 

write.csv(gprofilerresult_2WPI_01_nf3, "gprofilerresult_2WPI_P01_nf_fortables.csv") 
# This was used for Table 2 - with GO terms chosen based on relevance 

# Pathway at 2 wpi analysis using HGNC #### 

head(names_2WPI_01_HGNC)
names_2WPI_01 <- as.vector(names_2WPI_01)
library(gProfileR)
gprofilerresult_2WPI_01_nf_hg <- gprofiler(names_2WPI_01_HGNC, organism = "hsapiens", ordered_query = T, 
                                        significant = T, exclude_iea = F, underrep = F, evcodes = F, 
                                        region_query = F, max_p_value = 1, min_set_size = 0, 
                                        max_set_size = 0, min_isect_size = 0, 
                                        correction_method = "analytical", hier_filtering = "none", 
                                        domain_size = "annotated", custom_bg = "",numeric_ns = "", 
                                        png_fn = NULL, include_graph = F, src_filter = NULL) 

#hier_filtering = "none" ----> strong, apart from that default (ordered query true cause our genes are sorted by Padj value)
#gives output for gene ontology, kegg, reactome 

head(gprofilerresult_2WPI_01_nf_hg)  

gprofilerresult_2WPI_01_st_hg <- gprofiler(names_2WPI_01_HGNC, ordered_query = T, organism = "hsapiens", hier_filtering = "strong")

head(gprofilerresult_2WPI_01_st_hg)  

write.csv(gprofilerresult_2WPI_01_nf_hg , file="Gprofiler_2WPI_P01_nf_hg.csv")  
write.csv(gprofilerresult_2WPI_01_st_hg , file="Gprofiler_2WPI_P01_st_hg.csv")  

library(dplyr) 
gprofilerresult_2WPI_01_nf_hg2 <- gprofilerresult_2WPI_01_nf_hg %>% 
  filter(domain %in% c("BP", "keg", "rea")) %>% 
  arrange(domain, p.value) 

View(gprofilerresult_2WPI_01_nf_hg2) 
write.csv(gprofilerresult_2WPI_01_nf_hg2, "gprofilerresult_2WPI_P01_nf_hg_edited.csv")

gprofilerresult_2WPI_01_nf_hg3 <- gprofilerresult_2WPI_01_nf_hg2 %>% 
  select(domain, term.name, p.value, overlap.size) 
colnames(gprofilerresult_2WPI_01_nf_hg3) <- c("Pathway type", "Pathway", "P value", "Number of genes") 

write.csv(gprofilerresult_2WPI_01_nf_hg3, "gprofilerresult_2WPI_P01_nf_hg_fortables.csv") 