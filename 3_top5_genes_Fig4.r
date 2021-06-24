# Top 5 genes up and down genes from each time point ####
##### Infected vs Control PADJ < 0.1 #### 
head(res_0WPI_sigDF) 

library(tidyverse) 
library(RColorBrewer) 
library(ggrepel) 

res_0WPI_top5_up_01 <- res_0WPI_sig01_DF %>% filter(log2FoldChange > 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -pvalue) %>% 
  add_column(WPI = rep("0", 5), Direction = rep("Up", 5)) 
# AT 0 WPI there is only 1 gene up, so change the new vectors to 1  
res_0WPI_top5_down_01 <- res_0WPI_sig01_DF %>% filter(log2FoldChange < 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -pvalue) %>% 
  add_column(WPI = rep("0", 5), Direction = rep("Down", 5))

res_2WPI_top5_up_01 <- res_2WPI_sig01_DF %>% filter(log2FoldChange > 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -pvalue) %>% 
  add_column(WPI = rep("2", 5), Direction = rep("Up", 5)) 
# AT 0 WPI there is only 1 gene up, so change the new vectors to 1  
res_2WPI_top5_down_01 <- res_2WPI_sig01_DF %>% filter(log2FoldChange < 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -pvalue) %>% 
  add_column(WPI = rep("2", 7), Direction = rep("Down", 7)) 

#res_16WPI_top5_up_01 <- res_16WPI_sig01_DF %>% filter(log2FoldChange > 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -pvalue) %>% 
#  add_column(WPI = rep("16", 5), Direction = rep("Up", 5)) 
# 0 genes up - so no need for the above  
res_16WPI_top5_down_01 <- res_16WPI_sig01_DF %>% filter(log2FoldChange < 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -pvalue) %>% 
  add_column(WPI = rep("16", 2), Direction = rep("Down", 2))

res_top5_P01 <- rbind(res_0WPI_top5_up_01, res_0WPI_top5_down_01, res_2WPI_top5_up_01, res_2WPI_top5_down_01, res_16WPI_top5_down_01)
View(res_top5_P01) 
lapply(res_top5_P01, typeof) 
colnames(res_top5_P01) 
res_top5_P01$WPI <- as.numeric(as.character(res_top5_P01$WPI)) 
res_top5_P01$WPI <- as.factor(res_top5_P01$WPI) 

# Annotate top 5 genes by gene symbol and save #### 

library("AnnotationDbi")
library("OrganismDbi") 
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

# Normal genes 
res_top5_P01$Gene_description <- mapIds(org.Hs.eg.db,
                                    keys=res_top5_P01$gene_name,
                                    column="GENENAME",
                                    keytype="SYMBOL",
                                    multiVals="first")

res_top5_P01

res_top5_P01 <- apply(res_top5_P01,2,as.character)
write.csv(res_top5_P01, "top5_genes_pbmc_P01.csv")

# Change "NA" in Symbol/Alias with their Ensembl tag ####
View(res_top5_P01) 
library(dplyr)

res_top5_cp_01 <- res_top5_P01 %>%
  mutate(gene_name = coalesce(gene_name, Row.names))
res_top5_cp_01 
res_top5_cp_01$log2FoldChange <- as.numeric(as.character(res_top5_cp_01$log2FoldChange))

# plot - geom point with log2foldchange and Alias as a label ####
plotA1 <- ggplot(data = res_top5_cp_01) + 
  geom_point(mapping = aes(x = WPI, y = log2FoldChange, color = Direction, group = Direction)) + 
  scale_color_brewer(palette = "Set1", name="Direction", breaks = c("Up", "Down"), labels=c("Upregulated", "Downregulated"), direction = -1) + 
  xlab("Weeks post infection (wpi)") + ylab("Log2 fold change") + 
  geom_hline(mapping = aes(yintercept = 0), linetype = "dotted", size = 1) + 
  geom_text_repel(mapping = aes(x = WPI, y = log2FoldChange, group = Direction, label = gene_name), , fontface = "italic") 
ggsave("pbmc_top5_P01_2.jpeg")

######## Longitudinal #### 
head(res_0WPI_sigDF) 

library(tidyverse) 
library(RColorBrewer) 
library(ggrepel) 

# Infected 
res_2WPI_top5_up_inf <- res4_LFC_inf_2WPI_sigDF %>% filter(log2FoldChange > 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -pvalue) %>% 
  add_column(WPI = rep("2", 7), Direction = rep("Up", 7), Type = rep("Infected", 7)) 
# AT 0 WPI there is only 1 gene up, so change the new vectors to 1  
res_2WPI_top5_down_inf <- res4_LFC_inf_2WPI_sigDF %>% filter(log2FoldChange < 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -pvalue) %>% 
  add_column(WPI = rep("2", 5), Direction = rep("Down", 5), Type = rep("Infected", 5)) 

res_16WPI_top5_up_inf <- res4_LFC_inf_16WPI_sigDF %>% filter(log2FoldChange > 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -pvalue) %>% 
  add_column(WPI = rep("16", 5), Direction = rep("Up", 5), Type = rep("Infected", 5)) 
# 0 genes up - so no need for the above  
res_16WPI_top5_down_inf <- res4_LFC_inf_16WPI_sigDF %>% filter(log2FoldChange < 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -pvalue) %>% 
  add_column(WPI = rep("16", 5), Direction = rep("Down", 5), Type = rep("Infected", 5)) 

# Difference 
head(res4_LFC_diff_2WPI_sigDF)
res_2WPI_top5_up_diff <- res4_LFC_diff_2WPI_sigDF %>% filter(log2FoldChange > 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -pvalue) %>% 
  add_column(WPI = rep("2", 5), Direction = rep("Up", 5), Type = rep("Difference", 5)) 
# AT 0 WPI there is only 1 gene up, so change the new vectors to 1  
res_2WPI_top5_down_diff <- res4_LFC_diff_2WPI_sigDF %>% filter(log2FoldChange < 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -pvalue) %>% 
  add_column(WPI = rep("2", 5), Direction = rep("Down", 5), Type = rep("Difference", 5)) 

res_16WPI_top5_up_diff <- res4_LFC_diff_16WPI_sigDF %>% filter(log2FoldChange > 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -pvalue) %>% 
  add_column(WPI = rep("16", 5), Direction = rep("Up", 5), Type = rep("Difference", 5)) 
# 0 genes up - so no need for the above  
res_16WPI_top5_down_diff <- res4_LFC_diff_16WPI_sigDF %>% filter(log2FoldChange < 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -pvalue) %>% 
  add_column(WPI = rep("16", 5), Direction = rep("Down", 5), Type = rep("Difference", 5))

# Control
res_2WPI_top5_up_cont <- res4_LFC_cont_2WPI_sigDF %>% filter(log2FoldChange > 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -pvalue) %>% 
  add_column(WPI = rep("2", 6), Direction = rep("Up", 6), Type = rep("Control", 6)) 
# AT 0 WPI there is only 1 gene up, so change the new vectors to 1  
res_2WPI_top5_down_cont <- res4_LFC_cont_2WPI_sigDF %>% filter(log2FoldChange < 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -pvalue) %>% 
  add_column(WPI = rep("2", 5), Direction = rep("Down", 5), Type = rep("Control", 5)) 

res_16WPI_top5_up_cont <- res4_LFC_cont_16WPI_sigDF %>% filter(log2FoldChange > 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -pvalue) %>% 
  add_column(WPI = rep("16", 5), Direction = rep("Up", 5), Type = rep("Control", 5)) 
# 0 genes up - so no need for the above  
res_16WPI_top5_down_cont <- res4_LFC_cont_16WPI_sigDF %>% filter(log2FoldChange < 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -pvalue) %>% 
  add_column(WPI = rep("16", 5), Direction = rep("Down", 5), Type = rep("Control", 5)) 
  
res_top5_long <- rbind(res_2WPI_top5_up_inf, res_2WPI_top5_down_inf, res_16WPI_top5_up_inf, res_16WPI_top5_down_inf, 
                       res_2WPI_top5_up_diff, res_2WPI_top5_down_diff, res_16WPI_top5_up_diff, res_16WPI_top5_down_diff, 
                       res_2WPI_top5_up_cont, res_2WPI_top5_down_cont, res_16WPI_top5_up_cont, res_16WPI_top5_down_cont)
View(res_top5_long) 
lapply(res_top5_long, typeof) 
res_top5_long$WPI <- as.numeric(as.character(res_top5_long$WPI)) 
res_top5_long$WPI <- as.factor(res_top5_long$WPI) 

# Annotate top 5 genes by gene symbol and save ####
# Not necessary here 
library("AnnotationDbi")
library("OrganismDbi") 
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

# Normal genes 
res_top5_long$Gene_description <- mapIds(org.Hs.eg.db,
                                    keys=res_top5_long$gene_name,
                                    column="GENENAME",
                                    keytype="SYMBOL",
                                    multiVals="first")

res_top5_long

res_top5_long <- apply(res_top5_long,2,as.character)
write.csv(res_top5_long, "top5_genes_pbmc_long.csv")

# Change "NA" in Symbol with their Ensembl tag ####
View(res_top5_long)
View(res_top5) 

library(dplyr)

res_top5_cp_long <- res_top5_long %>%
  mutate(gene_name = coalesce(gene_name, Row.names))
res_top5_cp_long 
res_top5_cp_long$Type <- factor(res_top5_cp_long$Type, levels = c("Infected", "Control", "Difference"))

# plot - geom point with log2foldchange and Alias as a label ####
plotB1 <- ggplot(data = res_top5_cp_long) + 
  geom_point(mapping = aes(x = WPI, y = log2FoldChange, color = Direction, group = Direction)) + 
  scale_color_brewer(palette = "Set1", name="Direction", breaks = c("Up", "Down"), labels=c("Upregulated", "Downregulated"), direction = -1) + 
  xlab("Weeks post infection (wpi)") + ylab("Log2 fold change") + 
  geom_hline(mapping = aes(yintercept = 0), linetype = "dotted", size = 1) + 
  geom_text_repel(mapping = aes(x = WPI, y = log2FoldChange, group = Direction, label = gene_name), , fontface = "italic") + 
  facet_wrap(~res_top5_cp_long$Type)
ggsave("pbmc_top5_long_faceted.jpeg") 

##### Combine plots for paper #### 

library(ggpubr) 

ggarrange(plotA1, plotB1, nrow = 2, labels = c("A", "B"), common.legend = TRUE)
ggsave("Top5_genes_plotAB.jpeg") 
#Saving 14.8 x 10.3 in image 

# Split plotB1 into 3 plots 
plotB1_inf <- ggplot(data = subset(res_top5_cp_long, Type == "Infected")) + 
  geom_point(mapping = aes(x = WPI, y = log2FoldChange, color = Direction, group = Direction)) + 
  scale_color_brewer(palette = "Set1", name="Direction", breaks = c("Up", "Down"), labels=c("Upregulated", "Downregulated"), direction = -1) + 
  xlab("Weeks post infection (wpi)") + ylab("Log2 fold change") + 
  geom_hline(mapping = aes(yintercept = 0), linetype = "dotted", size = 1) + 
  geom_text_repel(mapping = aes(x = WPI, y = log2FoldChange, group = Direction, label = gene_name), fontface = "italic") + 
  ggtitle("Infected")
plotB1_inf 

plotB1_cont <- ggplot(data = subset(res_top5_cp_long, Type == "Control")) + 
  geom_point(mapping = aes(x = WPI, y = log2FoldChange, color = Direction, group = Direction)) + 
  scale_color_brewer(palette = "Set1", name="Direction", breaks = c("Up", "Down"), labels=c("Upregulated", "Downregulated"), direction = -1) + 
  xlab("Weeks post infection (wpi)") + ylab("Log2 fold change") + 
  geom_hline(mapping = aes(yintercept = 0), linetype = "dotted", size = 1) + 
  geom_text_repel(mapping = aes(x = WPI, y = log2FoldChange, group = Direction, label = gene_name), fontface = "italic") + 
  ggtitle("Control")
plotB1_cont

plotB1_diff <- ggplot(data = subset(res_top5_cp_long, Type == "Difference")) + 
  geom_point(mapping = aes(x = WPI, y = log2FoldChange, color = Direction, group = Direction)) + 
  scale_color_brewer(palette = "Set1", name="Direction", breaks = c("Up", "Down"), labels=c("Upregulated", "Downregulated"), direction = -1) + 
  xlab("Weeks post infection (wpi)") + ylab("Log2 fold change") + 
  geom_hline(mapping = aes(yintercept = 0), linetype = "dotted", size = 1) + 
  geom_text_repel(mapping = aes(x = WPI, y = log2FoldChange, group = Direction, label = gene_name), fontface = "italic") + 
  ggtitle("Difference")
plotB1_diff 

plotB1_diff2 <- ggplot(data = subset(res_top5_cp_long, Type == "Difference")) + 
  geom_point(mapping = aes(x = WPI, y = log2FoldChange, color = Direction, group = Direction)) + 
  scale_color_brewer(palette = "Set1", name="Direction", breaks = c("Up", "Down"), labels=c("Upregulated", "Downregulated"), direction = -1) + 
  xlab("Weeks post infection (wpi)") + ylab("Log2 fold change") + 
  geom_hline(mapping = aes(yintercept = 0), linetype = "dotted", size = 1) + 
  geom_text_repel(mapping = aes(x = WPI, y = log2FoldChange, group = Direction, label = gene_name), fontface = "italic") 
plotB1_diff2

# ** Final combined plot for paper - Figure 4 #### 
ggarrange(plotA1, plotB1_diff2, ncol = 1, nrow = 2, labels = c("A","B"), common.legend = TRUE)

ggsave("top5_finalAB.jpeg", width = 20, height = 24, units = "cm", dpi = 1200)
