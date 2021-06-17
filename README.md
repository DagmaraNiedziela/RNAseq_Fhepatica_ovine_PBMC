# RNAseq_Fhepatica_ovine_PBMC
Analysis of an RNAseq dataset from sheep infected with liver fluke Fasciola hepatica. 
Sixteen male Merino sheep were randomly assigned to infected or control groups (n = 8 per group) and orally infected with 120 F. hepatica metacercariae. Transcriptomic data was generated from peripheral blood mononuclear cells (PBMC) at 0, 2 and 16 weeks post-infection (wpi), and analysed for differentially expressed (DE) genes between infected and control animals at each time point (analysis 1), and for each group relative to time 0 (analysis 2). Analysis 2 was then compared to a similar study performed previously on bovine PBMC. 
Analyses performed in UNIX/bash: quality control - FastQC, trimming - fastp, alignment to ovine genome - STAR, cell composition - QuanTIseq. 
Analyses in R: DE gene analysis - DESeq2. Pathway analysis - gProfileR. Figures: ggplot2, pheatmap, VennDiagram, UpSetR. 
