library(readxl) 

# FASTQC results ####
merged_fastqc_summaries <- read_excel("merged_fastqc_summaries_ovine_PBMC.xlsx", sheet = 1)
View(merged_fastqc_summaries)
colnames(merged_fastqc_summaries)

library(tidyr)
pivoted_fastqc_summaries <- pivot_wider(merged_fastqc_summaries, id_cols = Sample, 
                                        names_from = Statistic, values_from = Result)
View(pivoted_fastqc_summaries)

write.csv(pivoted_fastqc_summaries, "Pivoted_fastqc_summaries.csv") 

# FASTP results #### 
library(tidyverse) 
library(dplyr)

merged_fastp_summaries <- read_excel("merged_fastp_summaries_ovine_PBMC.xlsx", sheet = 2)
View(merged_fastp_summaries) 
colnames(merged_fastp_summaries[,1:48]) # This gives me column names to replace from higher up file if needed

# Remove blank cells 
sapply(merged_fastp_summaries,as.character)

merged_fastp_summaries2 <- as.data.frame(merged_fastp_summaries)
merged_fastp_summaries2 <- lapply(1:nrow(merged_fastp_summaries2), function(x) merged_fastp_summaries2[x,][is.na(merged_fastp_summaries2[x,])==F])
merged_fastp_summaries2 <- as.data.frame(merged_fastp_summaries2)
View(merged_fastp_summaries2)
ncol(merged_fastp_summaries2) # 43
nrow(merged_fastp_summaries2) # 48 
# It looks like the lapply transposed my data and sample names are gone 

# Transpose - not done - this only if I lapply again, df gets transposed again
# merged_fastp_summaries_trans <- t(merged_fastp_summaries2) 
# rownames(merged_fastp_summaries_trans) <- NULL
# View(merged_fastp_summaries_trans)

# Add sample name column to the front 
samples <- colnames(merged_fastp_summaries[,1:48])
samples 
# Clean the untidy folder name 
samples <- gsub("/home/workspace/dniedziela/ovine_PBMC/fastp_results/","",samples)
samples <- gsub(".json","",samples)
merged_fastp_summaries3 <- cbind(samples,merged_fastp_summaries2) 
View(merged_fastp_summaries3)

# Add column names 
nrow(merged_fastp_summaries3) #48 
colnames(merged_fastp_summaries3) <- c("Sample_name", "total_reads_before_filtering", 
                                       "total_reads_after_filtering","total_reads_read1_before",
                                       "total_reads_read1_after","total_reads_read2_before",
	"total_reads_read2_after","total_bases_before_filtering","total_bases_after_filtering",
	"total_bases_read1_before","total_bases_read1_after","total_bases_read2_before","total_bases_read2_after",
	"q20_rate_before","q20_rate_after","q30_rate_before","q30_rate_after",
	"gc_content_before","gc_content_after","read1_mean_length_before","read1_mean_length_after",
	"read2_mean_length_before","read2_mean_length_after","passed_filter_reads","corrected_reads",
	"corrected_bases","low_quality_reads","too_many_N_reads","too_short_reads","too_long_reads",
	"q20_bases_overall_before","q20_bases_overall_after","q20_bases_read1_before","q20_bases_read1_after",
	"q20_bases_read2_before","q20_bases_read2_after","q30_bases_overall_before","q30_bases_overall_after",
	"q30_bases_read1_before","q30_bases_read1_after","q30_bases_read2_before","q30_bases_read2_after",
	"adapter_trimmed_reads","adapter_trimmed_bases")
nrow()

# Remove text from all columns - gsub in a transmute_all
# > gsub("cheap", "sheep's", "A wolf in cheap clothing")
# [1] "A wolf in sheep's clothing" --- pattern, replacement, text to modify 
# pattern *reads:, *bases:, *rate:, *content:, *length:

ncol(merged_fastp_summaries3) #44
substitute <- function(x)(gsub(".*reads:","", x))
substitute2 <- function(x)(gsub(".*bases:","", x))
substitute3 <- function(x)(gsub(".*rate:","", x))
substitute4 <- function(x)(gsub(".*content:","", x))
substitute5 <- function(x)(gsub(".*length:","", x))
substitute6 <- function(x)(gsub(",","", x))

merged_fastp_summaries4 <- merged_fastp_summaries3 %>% 
  transmute_all(substitute) %>% 
  transmute_all(substitute2) %>%
  transmute_all(substitute3) %>%
  transmute_all(substitute4) %>%
  transmute_all(substitute5) %>%
  transmute_all(substitute6) 

View(merged_fastp_summaries4)
merged_fastp_summaries4 <- as_tibble(merged_fastp_summaries4)
merged_fastp_summaries4 <- merged_fastp_summaries4 %>% 
  mutate_at(colnames(merged_fastp_summaries4[,2:44]),function(x)as.numeric(as.character(x)))

merged_fastp_summaries4
colnames(merged_fastp_summaries4)
ncol(merged_fastp_summaries4) #44

# Calculate q20 and q30 rate before and after for read 1 and read 2 --- so q20 bases/total bases*100, use mutate 
merged_fastp_summaries_final <- merged_fastp_summaries4 %>% 
  mutate(q20_rate_read1_before=q20_bases_read1_before/total_bases_read1_before*100,
         q20_rate_read1_after=q20_bases_read1_after/total_bases_read1_after*100,
         q30_rate_read1_before=q30_bases_read1_before/total_bases_read1_before*100,
         q30_rate_read1_after=q30_bases_read1_after/total_bases_read1_after*100,
         q20_rate_read2_before=q20_bases_read2_before/total_bases_read2_before*100,
         q20_rate_read2_after=q20_bases_read2_after/total_bases_read2_after*100,
         q30_rate_read2_before=q30_bases_read2_before/total_bases_read2_before*100,
         q30_rate_read2_after=q30_bases_read2_after/total_bases_read2_after*100,
         percent_trimmed=100-(total_bases_after_filtering/total_bases_before_filtering*100),
         percent_trimmed_read1=100-(total_bases_read1_after/total_bases_read1_before*100), 
         percent_trimmed_read2=100-(total_bases_read2_after/total_bases_read2_before*100))

View(merged_fastp_summaries_final)
ncol(merged_fastp_summaries_final) #46
colnames(merged_fastp_summaries_final)
# Should be: q20_rate_read1_before=q20_bases_read1_before/total_bases_read1_before*100
# Will have to grep "total bases" from the json files 
# Also need total bases before and after to get % trimmed value 

# Save a csv 
write.csv(merged_fastp_summaries_final,"merged_fastp_summaries_final.csv")

# FASTQC results trimmed #### 
# Need to write column names in Excel 
merged_fastqc_summaries_trimmed <- read_excel("merged_fastqc_summaries_ovine_PBMC.xlsx", sheet = 3)
View(merged_fastqc_summaries_trimmed)
colnames(merged_fastqc_summaries_trimmed)

library(tidyr)
pivoted_fastqc_summaries_trimmed <- pivot_wider(merged_fastqc_summaries_trimmed, id_cols = Sample, 
                                        names_from = Statistic, values_from = Result)
View(pivoted_fastqc_summaries_trimmed)

write.csv(pivoted_fastqc_summaries_trimmed, "Pivoted_fastqc_summaries_trimmed.csv") 

# END #### 