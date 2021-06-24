# DE gene analysis file 2 - this will make plots of DE gene counts for both analysis 1 - infected vs control 
# and analtysis 2 - longitudinal 
# Aim of this file - to create a data table and then plot DE genes 

####### Infected vs Control ##### 

# P < 0.1, Get DE numbers 

# Reminder of extracted values 
res_0WPI <- results(dds_group3, name="WPI0.GroupInfected")
res_2WPI <- results(dds_group3, name="WPI2.GroupInfected") 
res_16WPI <- results(dds_group3, name="WPI16.GroupInfected") 

sum(res_0WPI$padj < 0.1, na.rm=TRUE) # 59
sum(res_2WPI$padj < 0.1, na.rm=TRUE) # 453
sum(res_16WPI$padj < 0.1, na.rm=TRUE) # 2 

res_0WPI_sig01 <- subset(resLFC_0WPI, padj < 0.1)
res_2WPI_sig01 <- subset(resLFC_2WPI, padj < 0.1)
res_16WPI_sig01 <- subset(resLFC_16WPI, padj < 0.1)

# Get DE numbers! 
#All DE genes 

a <- nrow(subset(res_0WPI_sig01, log2FoldChange > 0))
b <- nrow(subset(res_0WPI_sig01, log2FoldChange < 0))
c <- nrow(subset(res_2WPI_sig01, log2FoldChange > 0))
d <- nrow(subset(res_2WPI_sig01, log2FoldChange < 0))
e <- nrow(subset(res_16WPI_sig01, log2FoldChange > 0))
f <- nrow(subset(res_16WPI_sig01, log2FoldChange < 0))
Count <- c(a,b,c,d,e,f)
Count 
All_labels <- c("0_up", "0_down", "2_up", "2_down", "16_up", "16_down") 
WPI <- c(0,0,2,2,16,16)
Direction <- c("Up", "Down")
Direction <- rep(Direction, 3)
LFnumbers_pbmc_P01 <- cbind(All_labels, Count, WPI, Direction)
LFnumbers_pbmc_P01

# Plot DE gene amounts #####

write.csv(LFnumbers_pbmc_P01, "LFnumbers_pbmc_P01.csv") 

LFnumbers2 <- read.csv("LFnumbers_pbmc_P01.csv")
  LFnumbers2$WPI <- as.factor(LFnumbers2$WPI)

# To change plot order of labels,
# change the order of varible levels with factor()
LFnumbers2$Direction <- factor(LFnumbers2$Direction, levels = c("Up", "Down")) 

library(ggplot2) 
library(RColorBrewer) 

LFnumbers2
ggplot(data = LFnumbers2) + geom_col(mapping = aes(x = WPI, y = Count, fill = WPI)) + 
  facet_wrap(~ LFnumbers2$Direction) + 
  ggtitle("Significant DE genes (FDR < 0.05)") + 
  scale_fill_brewer(palette = "Set1", direction = -1)
ggsave("DE_genes_facet.jpeg")

# Up and downregulated gene subsets 
LFnumbers2_Up <- subset(LFnumbers2, Direction == "Up") 
LFnumbers2_Down <- subset(LFnumbers2, Direction == "Down") 

# ** My pet graph - this into presentation / paper ####  
LFnumbers2_Up
plotA <- ggplot(data = LFnumbers2_Up) + 
  geom_col(data = LFnumbers2_Up, 
           mapping = aes(x = WPI, y = Count, fill = Direction, group = 1), position = "dodge") +
  geom_text(data = LFnumbers2_Up, 
            mapping = aes(x = WPI, y = Count, label = Count, group = 1), vjust = -0.6) +
  geom_col(data = LFnumbers2_Down, 
           mapping = aes(x = WPI, y = -Count, fill = Direction), position = "dodge") + 
  geom_text(data = LFnumbers2_Down, 
            mapping = aes(x = WPI, y = -Count, label = Count), vjust = 1.2) + 
  scale_fill_brewer(palette = "Set1", direction = -1,
                    name="Direction",
                    breaks=c("Up", "Down"),
                    labels=c("Upregulated", "Downregulated")) + 
  scale_y_continuous(expand = c(.1, .1)) +
  xlab("Weeks post infection (wpi)") + 
  ylab("Number of genes") + 
  theme(text = element_text(size = 14, family = "Calibri")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

plotA
ggsave("DE_genes_pbmc_P01.jpeg") 

########### LONGITUDINAL #####

# Get DE gene numbers 
# Reminder of extracted values 

res4_inf_2WPI <- results(dds4, contrast = c("condition", "Infected_2_WPI", "Infected_0_WPI"), alpha = 0.05)
sum(res4_inf_2WPI$padj < 0.05, na.rm=TRUE) # 1522 
res4_inf_16WPI <- results(dds4, contrast = c("condition", "Infected_16_WPI", "Infected_0_WPI")) 
sum(res4_inf_16WPI$padj < 0.05, na.rm=TRUE) # 3949

res4_cont_2WPI <- results(dds4, contrast = c("condition", "Control_2_WPI", "Control_0_WPI"))
sum(res4_cont_2WPI$padj < 0.05, na.rm=TRUE) # 769
res4_cont_16WPI <- results(dds4, contrast = c("condition", "Control_16_WPI", "Control_0_WPI"))
sum(res4_cont_16WPI$padj < 0.05, na.rm=TRUE) # 2342 

res4_LFC_inf_2WPI_sig <- subset(res4_LFC_inf_2WPI, padj < 0.05)
res4_LFC_inf_16WPI_sig <- subset(res4_LFC_inf_16WPI, padj < 0.05)
res4_LFC_cont_2WPI_sig <- subset(res4_LFC_cont_2WPI, padj < 0.05)
res4_LFC_cont_16WPI_sig <- subset(res4_LFC_cont_16WPI, padj < 0.05)

res4_LFC_inf_2WPI_sigDF <- as.data.frame(subset(res4_LFC_inf_2WPI_ann, padj < 0.05))
res4_LFC_inf_16WPI_sigDF <- as.data.frame(subset(res4_LFC_inf_16WPI_ann, padj < 0.05)) 
res4_LFC_cont_2WPI_sigDF <- as.data.frame(subset(res4_LFC_cont_2WPI_ann, padj < 0.05)) 
res4_LFC_cont_16WPI_sigDF <- as.data.frame(subset(res4_LFC_cont_16WPI_ann, padj < 0.05))

res4_LFC_diff_2WPI_sigDF <- anti_join(res4_LFC_inf_2WPI_sigDF,res4_LFC_cont_2WPI_sigDF, by = "Row.names")
nrow(res4_LFC_diff_2WPI_sigDF) # 1148 

res4_LFC_diff_16WPI_sigDF <- anti_join(res4_LFC_inf_16WPI_sigDF,res4_LFC_cont_16WPI_sigDF, by = "Row.names")
nrow(res4_LFC_diff_16WPI_sigDF) # 1927  
head(res4_LFC_diff_16WPI_sigDF)

# Get DE numbers! 
#All DE genes 

a3 <- nrow(subset(res4_LFC_inf_2WPI_sig, log2FoldChange > 0))
b3 <- nrow(subset(res4_LFC_inf_2WPI_sig, log2FoldChange < 0))
c3 <- nrow(subset(res4_LFC_inf_16WPI_sig, log2FoldChange > 0))
d3 <- nrow(subset(res4_LFC_inf_16WPI_sig, log2FoldChange < 0))
e3 <- nrow(subset(res4_LFC_cont_2WPI_sig, log2FoldChange > 0))
f3 <- nrow(subset(res4_LFC_cont_2WPI_sig, log2FoldChange < 0))
g3 <- nrow(subset(res4_LFC_cont_16WPI_sig, log2FoldChange > 0))
h3 <- nrow(subset(res4_LFC_cont_16WPI_sig, log2FoldChange < 0))
i3 <- nrow(subset(res4_LFC_diff_2WPI_sigDF, log2FoldChange > 0))
j3 <- nrow(subset(res4_LFC_diff_2WPI_sigDF, log2FoldChange < 0))
k3 <- nrow(subset(res4_LFC_diff_16WPI_sigDF, log2FoldChange > 0))
l3 <- nrow(subset(res4_LFC_diff_16WPI_sigDF, log2FoldChange < 0))
a3

Count3 <- c(a3,b3,c3,d3,e3,f3,g3,h3,i3,j3,k3,l3)
Count3
length(Count3)
All_labels3 <- c("2_inf_up", "2_inf_down", "16_inf_up", "16_inf_down", 
                "2_cont_up", "2_cont_down", "16_cont_up", "16_cont_down", 
                "2_diff_up", "2_diff_down", "16_diff_up", "16_diff_down") 
Group <- c((rep("Infected",4)),(rep("Control",4)),(rep("Difference",4)))
Group
WPI3 <- (rep(c(2,2,16,16),3))
WPI3
Direction3 <- c("Up", "Down")
Direction3 <- rep(Direction3, 6)
LFnumbers_pbmc3 <- cbind(All_labels3, Count3, WPI3, Direction3, Group)
LFnumbers_pbmc3

# Plot DE gene amounts - condition #####

write.csv(LFnumbers_pbmc3, "LFnumbers_pbmc_condition.csv") 
LFnumbers3 <- read.csv("LFnumbers_pbmc_condition.csv")
 
LFnumbers3$WPI3 <- as.factor(LFnumbers3$WPI3)

# To change plot order of facet wrap,
# change the order of varible levels with factor()
LFnumbers3$Direction3 <- factor(LFnumbers3$Direction3, levels = c("Up", "Down")) 
LFnumbers3$Group <- factor(LFnumbers3$Group, levels = c("Infected", "Control", "Difference"))

library(ggplot2) 
library(RColorBrewer) 

LFnumbers3
ggplot(data = LFnumbers3) + geom_col(mapping = aes(x = WPI3, y = Count3, fill = Direction3, group = Direction3), position = "dodge") + 
  facet_wrap(~ LFnumbers3$Group) + 
  ggtitle("Significant DE genes (FDR < 0.05)") + 
  scale_fill_brewer(palette = "Set1", direction = -1)
ggsave("DE_genes_facet_condition.jpeg")

# Up and downregulated gene subsets 
LFnumbers3_Up <- subset(LFnumbers3, Direction == "Up") 
LFnumbers3_Down <- subset(LFnumbers3, Direction == "Down") 

# ** My pet graph - this into presentation / paper ####  
LFnumbers3_Up
plotB <- ggplot(data = LFnumbers3_Up) + 
  geom_col(data = LFnumbers3_Up, 
           mapping = aes(x = WPI3, y = Count3, fill = Direction3, group = 1), position = "dodge") +
  geom_text(data = LFnumbers3_Up, 
            mapping = aes(x = WPI3, y = Count3, label = Count3, group = 1), vjust = -0.6) +
  geom_col(data = LFnumbers3_Down, 
           mapping = aes(x = WPI3, y = -Count3, fill = Direction3), position = "dodge") + 
  geom_text(data = LFnumbers3_Down, 
            mapping = aes(x = WPI3, y = -Count3, label = Count3), vjust = 1.2) + 
  facet_wrap(~ LFnumbers3_Up$Group) + 
  scale_fill_brewer(palette = "Set1", direction = -1,
                    name="Direction",
                    breaks=c("Up", "Down"),
                    labels=c("Upregulated", "Downregulated")) + 
  scale_y_continuous(expand = c(.1, .1)) +
  xlab("Weeks post infection (wpi)") + 
  ylab("Number of genes") + 
  theme(text = element_text(size = 14, family = "Calibri")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

plotB
ggsave("DE_genes_pbmc_condition.jpeg") 

##### Combination of DE gene plots for paper - Figure 3 #### 
# It will be combined with Venn diagrams from file 8_biomarkers_and_Venn_diagrams.R 

library(ggpubr)

plotDE_AB <- ggarrange(plotA, plotB, ncol = 1, nrow = 2, common.legend = TRUE,
                       labels = c("A","C"))

plotDE_AB
ggsave("DE_Genes_AB.jpeg", width = 15, height = 24, units = "cm", dpi = 1200)

# Combine with Venn diagrams ####
library(ggplot2)
library(gridExtra)

# Give the Venns labels A B etc 
myplot1 <- arrangeGrob(grobTree(venn.plot_P01), top = textGrob("B", x = unit(0, "npc")
                                               , y   = unit(-0.7, "npc"), just=c("left","top"),
                                               gp=gpar(col="black", fontsize=16, fontfamily="Arial",fontface="bold")))

myplot2 <- arrangeGrob(grobTree(venn_2WPI_DR), top = textGrob("D", x = unit(0, "npc")
                                               , y = unit(-0.2, "npc"), just=c("left","top"),
                                               gp=gpar(col="black", fontsize=16, fontfamily="Arial", fontface="bold")))
# position of the label letter is in y

grid.arrange(plotDE_AB,arrangeGrob(myplot1, myplot2),
             ncol=2) 
# The individual plots before combining have been corrected for whitespace
# ggplot extend within y axis, venn margin 

# Writing to file for paper
?tiff
tiff(filename = "Fig3_Combined_DE_Venn.tiff", compression = "lzw", width = 10, height = 10, 
     units = "in", res = 300);
grid.arrange(plotDE_AB,arrangeGrob(myplot1, myplot2),
             ncol=2);
dev.off(); 

# pretty png file 
library(Cairo)
Cairo(file="Combined_DE_Venn.png", 
      type="png",
      units="in", 
      width=10, 
      height=10, 
      pointsize=12, 
      dpi=600)
grid.arrange(plotDE_AB,arrangeGrob(myplot1, myplot2),
             ncol=2);
dev.off() 

