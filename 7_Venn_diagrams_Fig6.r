#####Venn diagrams ##### 
library(VennDiagram) 
library(tidyverse)
library(Cairo)
library(gridExtra)
library(grid)

###### Infected vs Control P01 #### 
library(VennDiagram) 

# Get gene lists #####

# Reminder of what the data is: 
res_0WPI_sig01 <- subset(resLFC_0WPI, padj < 0.1)
res_2WPI_sig01 <- subset(resLFC_2WPI, padj < 0.1)
res_16WPI_sig01 <- subset(resLFC_16WPI, padj < 0.1)

res_0WPI_sig01_DF <- as.data.frame(subset(res_0WPI_ann, padj < 0.1)) 
res_2WPI_sig01_DF <- as.data.frame(subset(res_2WPI_ann, padj < 0.1)) 
res_16WPI_sig01_DF <- as.data.frame(subset(res_16WPI_ann, padj < 0.1)) 
# check one to see if it's ordered, otherwise order by Padj (but this should have been done in pathway analysis file) 
head(res_0WPI_sig01_DF) 

# extract row names as a column 
library(tidyverse)

y0WPI_sig01 <- res_0WPI_sig01_DF %>% dplyr::select(Row.names) 
y0WPI_sig01 <- as.vector(y0WPI_sig01)
head(y0WPI_sig01)
length(y0WPI_sig01) #1
nrow(y0WPI_sig01) #59

y2WPI_sig01 <- res_2WPI_sig01_DF %>% dplyr::select(Row.names) 
y2WPI_sig01 <- as.vector(y2WPI_sig01)
nrow(y2WPI_sig01) #453

y16WPI_sig01 <- res_16WPI_sig01_DF %>% dplyr::select(Row.names) 
y16WPI_sig01 <- as.vector(y16WPI_sig01) 
head(y16WPI_sig01)
nrow(y16WPI_sig01) #2 

### Create overlap vectors ##### 

y0_2_P01 <- dplyr::intersect(y0WPI_sig01,y2WPI_sig01) 
nrow(y0_2_P01) #7
#at the end, row bind all overlaps with time 0, annotate and save 
y0_16_P01 <- dplyr::intersect(y0WPI_sig01,y16WPI_sig01) 
nrow(y0_16_P01) #1

y2_16_P01 <- dplyr::intersect(y2WPI_sig01, y16WPI_sig01)
nrow(y2_16_P01) #1 

overlaps_with_0WPI_sig01 <- rbind(y0_2_P01, y0_16_P01) 
overlaps_with_0WPI_sig01 
head(biotype2)
biotype3 <- biotype2 %>% rownames_to_column()
class(biotype2)
overlaps_with_0WPI_P01_DF <- merge(overlaps_with_0WPI_sig01, biotype2, by.x=1, by.y=0, all=FALSE)
?merge
head(overlaps_with_0WPI_P01_DF)
write.csv(overlaps_with_0WPI_P01_DF, "overlaps_with_0WPI_P01.csv") 

#for more lists than 2
# Reduce(intersect, list(a,b,c)) 

y0_2_16_sig01 <- Reduce(dplyr::intersect, list(y0WPI_sig01,y2WPI_sig01,y16WPI_sig01))

###Draw Venn diagrams - group3 ######

#Venn diagram with 3 circles 
# This is Figure 3B 

grid.newpage()
venn.plot_P01 <- draw.triple.venn(area1 = nrow(y0WPI_sig01), area2 = nrow(y2WPI_sig01), area3 = nrow(y16WPI_sig01), 
                                     n12 = nrow(y0_2_P01), n23 = nrow(y2_16_P01), n13 = nrow(y0_16_P01), n123 = nrow(y0_2_16_sig01), 
                                     category = c("0 WPI", "2 WPI", "16 WPI"), lty = "solid", 
                                     fill= c("#a1d7f4","#2aa7ea","#0c15cc"),
                                     euler.d=FALSE,
                                     scaled=FALSE,
                                     sub=substitute( paste(bolditalic('M. bovis'))),
                                     sub.fontfamily = "Arial",
                                     sub.cex=1.2,
                                     sub.pos=c(0.5,1),
                                     lwd=1,
                                     alpha           = rep(0.50,3),
                                     label.col       = "#003333",
                                     cex             = 1.5,
                                     fontfamily      = "Arial",            
                                     cat.pos=c(-5,10,180),
                                     cat.dist =c(0.05,0.05,.05),
                                     cat.col         = "black",
                                     cat.cex         = 1.5,
                                     cat.fontfamily  = "Arial",
                                     cat.fontface    = 2,
                                     rotation.degree = 360,
                                     margin          = 0.1,
                                     height          = 10,
                                     width           = 4,
                                     units           = 'cm',
                                     compression     = 'lzw',
                                     resolution      = 1200
)
grid.draw(venn.plot_P01) 
grid.newpage()

# 2 Add Venn plot title 
require(gridExtra)
table <- gTree(children=venn.plot_P01)
venn.plot_P01 <- grid.arrange(table, top=textGrob("Ovine PBMC DE genes, FDR < 0.1", gp=gpar(fontsize=20)))

# Writing to file 
?tiff
tiff(filename = "Triple_Venn_P01.tiff", compression = "lzw", width = 12, height = 12, 
     units = "cm", res = 300);
grid.draw(venn.plot_P01);
dev.off(); 

# 3 pretty file 
library(Cairo)
Cairo(file="Triple_Venn_P01.png", 
      type="png",
      units="in", 
      width=10, 
      height=10, 
      pointsize=12, 
      dpi=300)
grid.draw(venn.plot_P01)
dev.off() 

###### Longitudinal #### 

library(VennDiagram) 

# Double Venn - infected vs control and longitudinal difference 
# This is Figure 3D 

# Reminder of what the data is: 

res4_LFC_diff_2WPI_sigDF <- anti_join(res4_LFC_inf_2WPI_sigDF,res4_LFC_cont_2WPI_sigDF, by = "Row.names") 
nrow(res4_LFC_diff_2WPI_sigDF) #1148
res4_LFC_diff_16WPI_sigDF <- anti_join(res4_LFC_inf_16WPI_sigDF,res4_LFC_cont_16WPI_sigDF, by = "Row.names")
nrow(res4_LFC_diff_16WPI_sigDF) #1927 

# extract row names as a column 
library(tidyverse)

y2WPI_diff <- res4_LFC_diff_2WPI_sigDF %>% dplyr::select(Row.names) 
y2WPI_diff <- as.vector(y2WPI_diff)
nrow(y2WPI_diff) #1148

y16WPI_diff <- res4_LFC_diff_16WPI_sigDF %>% dplyr::select(Row.names) 
y16WPI_diff <- as.vector(y16WPI_diff) 
head(y16WPI_diff)
nrow(y16WPI_diff) #1927  

y2WPI_sig01 <- res_2WPI_sig01_DF %>% dplyr::select(Row.names) 
y2WPI_sig01 <- as.vector(y2WPI_sig01)
nrow(y2WPI_sig01) #453

y16WPI_sig01 <- res_16WPI_sig01_DF %>% dplyr::select(Row.names) 
y16WPI_sig01 <- as.vector(y16WPI_sig01) 
head(y16WPI_sig01)
nrow(y16WPI_sig01) #2 

### Create overlap vectors ##### 

y2WPI_DR <- dplyr::intersect(y2WPI_diff, y2WPI_sig01) 
nrow(y2WPI_DR) #33
#at the end, row bind all overlaps with time 0, annotate and save 
y16WPI_DR <- dplyr::intersect(y16WPI_diff, y16WPI_sig01) 
nrow(y16WPI_DR) #0


### Double Venns - Analysis 1 vs Analysis 2 at 2 WPI ######
grid.newpage()
venn_2WPI_DR <- draw.pairwise.venn(nrow(y2WPI_diff), nrow(y2WPI_sig01), nrow(y2WPI_DR), category = c("Longitudinal", "Infected_vs_Control"),
                                   lty = rep("solid",2), fill = c("#ef3456", "#2aa7ea"), alpha = rep(0.5, 2),  euler.d=FALSE,
                                scaled=FALSE,
                                lwd=1,
                                sub=substitute( paste(bolditalic('24 HPI'))),
                                sub.fontfamily = "Arial",
                                sub.cex=1.2,
                                sub.pos=c(0.5,1),
                                label.col       = "#003333",
                                cex             = 1.75,
                                fontfamily      = "Arial",
                                cat.pos=c(0,0.25),
                                cat.col         = "black",
                                cat.cex         = 1.5,
                                cat.fontfamily  = "Arial",
                                cat.fontface    = 2,
                                rotation.degree = 0,
                                height          = 10,
                                width           = 4,
                                units           = 'cm',
                                compression     = 'lzw',
                                resolution      = 1200, 
                                cat.just=list(c(0.9,0) , c(0.2,0)),
                                margin = 0.2,
                                main="2 WPI",main.cex=2,
                                 main.fontfamily  = "Arial",
                                main.fontface    = 2) 

venn_2WPI_DR <- grid.arrange(gTree(children=venn_2WPI_DR), top=textGrob("2 WPI", gp=gpar(fontsize=25)))
grid.draw(venn_2WPI_DR)

# Writing to file 
?tiff
tiff(filename = "Double_Venn_2WPI_diff_vs_reg.tiff", compression = "lzw", width = 12, height = 12, 
     units = "cm", res = 300);
grid.draw(venn_2WPI_DR);
dev.off(); 

# 3 pretty file 
library(Cairo)
Cairo(file="Double_Venn_2WPI_diff_vs_reg.png", 
      type="png",
      units="in", 
      width=10, 
      height=10, 
      pointsize=12, 
      dpi=300)
grid.draw(venn_2WPI_DR)
dev.off() 

###### Difference versus cattle ########################################################## 

# extract row names as a column - human names! 
library(tidyverse)
head(sheep_DE_human_2WPI_NA)

y2WPI_sheep <- sheep_DE_human_2WPI_NA %>% dplyr::select(human_gene_stable_id) 
y2WPI_sheep <- as.vector(y2WPI_sheep)
nrow(y2WPI_sheep) #1082

y16WPI_sheep <- sheep_DE_human_16WPI_NA %>% dplyr::select(human_gene_stable_id) 
y16WPI_sheep <- as.vector(y16WPI_sheep)
nrow(y16WPI_sheep) #1815 

y2WPI_cow <- Andres_1WPI_NA %>% dplyr::select(human_gene_stable_id) 
y2WPI_cow <- as.vector(y2WPI_cow)
nrow(y2WPI_cow) #17

y16WPI_cow <- Andres_14WPI_NA %>% dplyr::select(human_gene_stable_id) 
y16WPI_cow <- as.vector(y16WPI_cow)
nrow(y16WPI_cow) #1420 

### Quad Venn - Cow and sheep ######
#Make Venns for: MOK023 separately, MOK124 separately, and for both at each time point? 
#Definitely export overlap lists? 

nrow(y2WPI_sheep) #17
nrow(y16WPI_sheep) #1420 
nrow(y2WPI_cow) #17
nrow(y16WPI_cow) #1420 

ycow_2_16 <- intersect(y2WPI_cow, y16WPI_cow) # this intersect requires the tidyverse!!!! 
ysheep_2_16 <- intersect(y2WPI_sheep, y16WPI_sheep)
ycs_2_2 <- intersect(y2WPI_cow, y2WPI_sheep)
ycs_16_16 <- intersect(y16WPI_cow, y16WPI_sheep)
yc2_s16 <- intersect(y2WPI_cow,  y16WPI_sheep)
ys2_c16 <- intersect(y2WPI_sheep, y16WPI_cow)

nrow(ycs_2_2) #3
nrow(ycs_16_16) #232

y_cowsheep_1234 <- Reduce(intersect, list(y2WPI_cow, y16WPI_cow, y2WPI_sheep, y16WPI_sheep)) 
nrow(y_cowsheep_1234) # 1 

y_cowsheep_123 <- Reduce(intersect, list(y2WPI_cow, y16WPI_cow, y2WPI_sheep))

y_cowsheep_234 <- Reduce(intersect, list(y16WPI_cow, y2WPI_sheep, y16WPI_sheep))

y_cowsheep_124 <- Reduce(intersect, list(y2WPI_cow, y16WPI_cow, y16WPI_sheep))

y_cowsheep_134 <- Reduce(intersect, list(y2WPI_cow, y2WPI_sheep, y16WPI_sheep))


### **Draw Venn diagram - cowsheep ######

grid.newpage()
venn.plot_cow_sheep <- draw.quad.venn(area1 = nrow(y2WPI_cow), area2 = nrow(y16WPI_cow), area3 = nrow(y2WPI_sheep), area4 = nrow(y16WPI_sheep), n12 = nrow(ycow_2_16), n13 = nrow(ycs_2_2), n14 = nrow(yc2_s16), n23 = nrow(ys2_c16), n24 = nrow(ycs_16_16), 
                              n34 = nrow(ysheep_2_16), n123 = nrow(y_cowsheep_123), n124 = nrow(y_cowsheep_124), n134 = nrow(y_cowsheep_134), n234 = nrow(y_cowsheep_234), n1234 = nrow(y_cowsheep_1234), category = c("Cattle acute", "Cattle chronic", "Sheep acute", "Sheep chronic"), lty = "solid", 
                              fill            = c("#c6e1ef",
                                                  "#2aa7ea",
                                                  "#ddccff",
                                                  "#9966ff"),
                              euler.d=FALSE,
                              scaled=FALSE,
                              sub=substitute( paste(bolditalic('M. bovis'))),
                              sub.fontfamily = "Arial",
                              sub.cex=1.2,
                              sub.pos=c(0.5,1),
                              lwd=1,
                              alpha           = rep(0.50,4),
                              label.col       = "#003333",
                              cex             = 1.5,
                              fontfamily      = "Arial",            
                              cat.pos=c(-5,10,10,0),
                              cat.col         = "black",
                              cat.cex         = 1.5,
                              cat.fontfamily  = "Arial",
                              cat.fontface    = 2,
                              rotation.degree = 360,
                              margin          = 0,
                              height          = 10,
                              width           = 4,
                              units           = 'cm',
                              compression     = 'lzw',
                              resolution      = 1200
)
grid.draw(venn.plot_cow_sheep) 
grid.newpage()

# Writing to file 
tiff(filename = "Quad_Venn_cow_sheep.tiff", compression = "lzw");
grid.draw(venn.plot_cow_sheep);
dev.off(); 

# 3 pretty file 
Cairo(file="Quad_Venn_cow_sheep.png", 
      type="png",
      units="in", 
      width=10, 
      height=10, 
      pointsize=12, 
      dpi=72)
grid.draw(venn.plot_cow_sheep)
dev.off() 

##### Combine the Upset plot and Venn #### 
# This is Figure 6 

library(gridExtra)

# Give the Venns labels A B etc 
myplot3 <- arrangeGrob(grobTree(venn.plot_cow_sheep), top = textGrob("A", x = unit(0, "npc")
                                                                     , y = unit(-0.2, "npc"), just=c("left","top"),
                                                                     gp=gpar(col="black", fontsize=16, fontfamily="Arial", fontface="bold")))


# Call the upset plot from file 7a 
plot_up_species2
grid.edit('arrange',name='arrange2')
myplot7 = grid.grab()
myplot7

# Upset with a B label 
myplot8 <- arrangeGrob(grobTree(myplot7), top = textGrob("B", x = unit(0, "npc")
                                                              , y = unit(-0.2, "npc"), just=c("left","top"),
                                                              gp=gpar(col="black", fontsize=16, fontfamily="Arial", fontface="bold")))

myplot8

grid.arrange(arrangeGrob(myplot3, layout_matrix = matrix(c(NA,1,NA),1), widths = c(0.35,1,0.35), ncol = 3), 
             myplot8,
             nrow=2, heights = c(0.6,1)) #this one

# Writing to file 
?tiff
tiff(filename = "Fig6_Combined_Venn_Upset.tiff", compression = "lzw", width = 15, height = 23, 
     units = "in", res = 600);
grid.arrange(arrangeGrob(myplot3, layout_matrix = matrix(c(NA,1,NA),1), widths = c(0.35,1,0.35), ncol = 3), 
             myplot8,
             nrow=2, heights = c(0.6,1)) #this one;
dev.off(); 

# png file 
library(Cairo)
Cairo(file="Combined_Venn_Upset.png", 
      type="png",
      units="in", 
      width=15, 
      height=23, 
      pointsize=12, 
      dpi=600)
grid.arrange(arrangeGrob(myplot3, layout_matrix = matrix(c(NA,1,NA),1), widths = c(0.35,1,0.35), ncol = 3), 
             myplot8,
             nrow=2, heights = c(0.6,1)) #this one;
dev.off() 