# Import cell results ####

dagmara_cell_results <- read.table("quanTIseq_cell_fractions.txt" , sep = "\t")
View(dagmara_cell_results) 
colnames(dagmara_cell_results) <- dagmara_cell_results[1,]
dagmara_cell_results <- dagmara_cell_results[-1,]
head(dagmara_cell_results)
# Sample is a row, column names are cell types 

# Merge with coldata? 
head(coldata_group) 
library(tidyverse)
cell_results <- full_join(dagmara_cell_results, coldata_group, by = "Sample")
View(cell_results)
head(cell_results)

# Create a file with long format 
?pivot_longer
cell_results_long <- pivot_longer(cell_results, cols = c("B.cells", "Macrophages.M1","Macrophages.M2","Monocytes", "Neutrophils","NK.cells",
                                                         "T.cells.CD4","T.cells.CD8","Tregs","Dendritic.cells","Other"), 
                                  names_to = "Cell_type", values_to = "Cell_proportion")
head(cell_results_long)

# ** Plot individual cell populations #### 
ggplot(data = cell_results_long) + 
  geom_col(data = cell_results_long, 
           mapping = aes(x = WPI, y = Cell_proportion, fill = Cell_type, color = Group), position = "dodge") + 
  facet_wrap(~cell_results_long$Sheep_ID) + 
  scale_fill_manual(values = my_cell_colours, 
                    name="Cell type") + 
  xlab("Weeks post infection (wpi)") + 
  ylab("Proportion of cells") + 
  theme(text = element_text(size = 14, family = "Calibri")) 
ggsave("Cell_composition_individual.jpeg")

library(ggplot2)
library(RColorBrewer)
brewer.pal.info 
# Pick a color palette 
display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE, 
                   colorblindFriendly=FALSE) 
display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE, 
                   colorblindFriendly=TRUE)

# Because here we don't want to see a time course, I set HPI to factor 
cell_results_long$WPI <- as.factor(cell_results_long$WPI)  
levels(cell_results_long$WPI) 

cell_results_long$Cell_type <- as.factor(cell_results_long$Cell_type)  
levels(cell_results_long$Cell_type) 
# If I need a time course later, I can use + scale_x_continuous(breaks = c(0,24,48,72,168)) 

# There is 11 cell types - make Other gray, have a manual scale 
# From http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=12 - this is the Paired palette: 
# '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928' 
# #ff7f00 is the orange for Other 
# grey #808080 
# '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#808080','#cab2d6','#6a3d9a','#ffff99','#b15928' 

my_cell_colours <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#808080','#cab2d6','#6a3d9a','#ffff99','#b15928') 

# Make cell averages ####
typeof(cell_results_long$Cell_proportion)
cell_results_long$Cell_proportion <- as.numeric(as.character(cell_results_long$Cell_proportion))
cell_averages <- cell_results_long %>% group_by(WPI, Group, Cell_type) %>% 
  summarize(Proportion_mean = mean(Cell_proportion, na.rm = TRUE)) 
View(cell_averages) 
write.csv(cell_averages, "cell_averages.csv")

# ** Average proportions ####
ggplot(data = cell_averages) + 
  geom_col(data = cell_averages, 
           mapping = aes(x = WPI, y = Proportion_mean, fill = Cell_type)) + 
  facet_wrap(~cell_averages$Group) + 
  scale_fill_manual(values = my_cell_colours, 
                    name="Cell type") + 
  xlab("Weeks post infection (wpi)") + 
  ylab("Proportion of cells") + 
  theme(text = element_text(size = 14, family = "Calibri")) 
?scale_color_manual

ggsave("Cell_proportions_group.jpeg") 

# ** Supplementary Figure 1 ####
ggplot(data = cell_averages) + 
  geom_col(data = cell_averages, 
           mapping = aes(x = WPI, y = Proportion_mean, fill = Cell_type)) + 
  facet_wrap(~cell_averages$Group) + 
  scale_fill_brewer(palette="Spectral", 
                    name="Cell type") + 
  xlab("Weeks post infection (wpi)") + 
  ylab("Proportion of cells") + 
  theme(text = element_text(size = 14, family = "Calibri")) 
# Set3 is another one with more than 11 colors, but it is pastel 

ggsave("Cell_proportions_group2.jpeg")   

# Check what is the proportion of Other group #####
cell_averages %>% 
  filter(Cell_type == "Other")

cell_averages %>% 
  filter(Cell_type == "Other") %>% 
  group_by(Group) %>% 
  summarize(mean = mean(Proportion_mean * 100), sd = sd(Proportion_mean *100))

###### Statistical analysis - T regs ####
head(cell_results_long)
levels(cell_results_long$Cell_type)
Tregs <- cell_results_long %>% filter(Cell_type=="Tregs")
Tregs 

table(Tregs$Group, Tregs$WPI)
#          0 2 16
# Control  8 8  8
# Infected 8 8  8 

# Balanced design 

ggboxplot(Tregs, x = "WPI", y = "Cell_proportion", color = "Group",
          palette = c("#00AFBB", "#E7B800")) 

ggsave("ANOVA-boxplot.jpeg") 

ggline(Tregs, x = "WPI", y = "Cell_proportion", color = "Group",
       add = c("mean_se", "dotplot"),
          palette = c("#00AFBB", "#E7B800")) 

ggsave("ANOVA-lineplot.jpeg") 

# Two way ANOVA 
res.aov2 <- aov(Cell_proportion ~ Group + WPI, data = Tregs)
summary(res.aov2)
# Df   Sum Sq  Mean Sq F value   Pr(>F)    
# Group        1 0.000705 0.000705    3.21   0.0801 .  
# WPI          2 0.010691 0.005346   24.33 7.65e-08 ***
#   Residuals   44 0.009666 0.000220                     
#---
 # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Two way ANOVA with interaction 
res.aov3 <- aov(Cell_proportion ~ Group * WPI, data = Tregs)
res.aov3 <- aov(Cell_proportion ~ Group + WPI + Group:WPI, data = Tregs)
summary(res.aov3)

# Df   Sum Sq  Mean Sq F value   Pr(>F)    
# Group        1 0.000705 0.000705   3.231   0.0794 .  
# WPI          2 0.010691 0.005346  24.497 8.89e-08 ***
#   Group:WPI    2 0.000502 0.000251   1.149   0.3266    
# Residuals   42 0.009165 0.000218 

# Non-significant interaction 

# Multiple comparisons 
TukeyHSD(res.aov3, which = "WPI") 

# Tukey multiple comparisons of means
# 95% family-wise confidence level

# Fit: aov(formula = Cell_proportion ~ Group + WPI + Group:WPI, data = Tregs)

# $WPI
# diff         lwr         upr     p adj
# 2-0   0.03569472  0.02300625 0.048383186 0.0000001
# 16-0  0.02468146  0.01199300 0.037369930 0.0000754
# 16-2 -0.01101326 -0.02370172 0.001675209 0.1002122

TukeyHSD(res.aov3, which = "Group:WPI") 
TukeyHSD(res.aov3, which = "Group")
?TukeyHSD

# Residuals 
plot(res.aov3, 1)

# Normality 
plot(res.aov3, 2)

# Extract the residuals
aov_residuals <- residuals(object = res.aov3)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals ) 

# W = 0.97112, p-value = 0.2799 

# ** No time 0, make it a covariate #### 
View(Tregs)
Tregs2 <- Tregs %>% mutate(Time0 = rep(Tregs$Cell_proportion[1:16],3)) %>%
  filter(!WPI == 0)
View(Tregs2)  

# Two way ANOVA 
res.aov4 <- aov(Cell_proportion ~ Group + WPI, data = Tregs2)
summary(res.aov4)
#             Df   Sum Sq   Mean Sq F value Pr(>F)  
# Group        1 0.001056 0.0010560   3.783 0.0615 .
# WPI          1 0.000970 0.0009703   3.476 0.0724 .
# Residuals   29 0.008096 0.0002792                      
#---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Two way ANOVA with interaction 
res.aov5 <- aov(Cell_proportion ~ Group * WPI + Time0, data = Tregs2)
res.aov5 <- aov(Cell_proportion ~ Group + WPI + Time0 + Group:WPI, data = Tregs2)
summary(res.aov5)

# Df   Sum Sq   Mean Sq F value Pr(>F)  
# Group        1 0.001056 0.0010560   4.514 0.0429 *
#   WPI          1 0.000970 0.0009703   4.148 0.0516 .
# Time0        1 0.001629 0.0016287   6.962 0.0137 *
#  Group:WPI    1 0.000151 0.0001508   0.645 0.4291  
# Residuals   27 0.006316 0.0002339 

# Non-significant interaction 

# Multiple comparisons 
TukeyHSD(res.aov5, which = "WPI") 

TukeyHSD(res.aov5, which = "Group:WPI") 

TukeyHSD(res.aov5, which = "Group")
?TukeyHSD

# Residuals 
plot(res.aov5, 1)

# Normality 
plot(res.aov5, 2)

# Extract the residuals
aov_residuals <- residuals(object = res.aov5)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals ) 

# W = 0.9804, p-value = 0.8112