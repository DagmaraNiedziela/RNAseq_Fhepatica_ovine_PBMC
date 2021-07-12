###### 30 10 2020 Start of ovine PBMC RNAseq pipeline ### 
# Dagmara Niedziela 

# Data is in separate folders, i.e. there are 3 overall data folders, followed by a separate folder for each sample, A1, A2 etc, and an MD5.txt file for each 
# I want to copy all fastq.gz files into one folder, and combine all MD5.txt files into one file, followed by an MD5 sum check after file copying. 

# FILE REARRANGEMENTS ########################## 
################################################

# Script to copy fastq files 
for file in `find /home/workspace/ccorreia/ovinePBMC_RNAseq/fastq/ \
-name *.fq.gz`; \
do echo "cp $file /home/workspace/dniedziela/ovine_PBMC/fastq_files/" >> copy_fast_files.sh; \
done 

. copy_fast_files.sh 

# Append all MD5 files into one 
for file in `find /home/workspace/ccorreia/ovinePBMC_RNAseq/fastq/ \
-name MD5.txt`; \
do 
echo "cat $file >> /home/workspace/dniedziela/ovine_PBMC/fastq_files/MD5_all.txt" >> md5_cat.sh; \ 
done 

wc -l md5_cat.sh # 48 lines 
. md5_cat.sh 

wc -l /home/workspace/dniedziela/ovine_PBMC/fastq_files/MD5_all.txt # 96 lines - all correct 

# Md5 sum check of the copied files 
md5sum -c /home/workspace/dniedziela/ovine_PBMC/fastq_files/MD5_all.txt >> /home/workspace/dniedziela/ovine_PBMC/fastq_files/md5check_ovine_PMBC_all.txt 
# All ok, file has 96 lines 

# FASTQC ################################
######################################### FastQC version 0.11.8

# Fastqc takes multiple arguments and can multithread so a script is not necessary 
nohup fastqc -o /home/workspace/dniedziela/ovine_PBMC/fastqc -t 20 /home/workspace/dniedziela/ovine_PBMC/fastq_files/*.fq.gz > /home/workspace/dniedziela/ovine_PBMC/fastqc/fastqc.nohup & 
htop 

#SUMMARY File for FASTQC ################
#########################################

#Location of code: /home/workspace/dniedziela/ovine_PBMC/fastqc
mkdir fastqc_summary 

for file in `ls *.zip`; \
do outfile=`basename $file | perl -p -e 's/\.zip//'`; \
echo -e "unzip /home/workspace/dniedziela/ovine_PBMC/fastqc/$file; \
cp /home/workspace/dniedziela/ovine_PBMC/fastqc/$outfile/summary.txt /home/workspace/dniedziela/ovine_PBMC/fastqc/fastqc_summary/${outfile}_summary.txt" >> /home/workspace/dniedziela/ovine_PBMC/scripts/fastqc_summary.sh; \
done 

sh /home/workspace/dniedziela/ovine_PBMC/scripts/fastqc_summary.sh #run shell script 
#after I do it I can use cat to combine all summaries into one file 
cd fastqc_summary 
ls 
echo * #to make sure it will list all files 
cat * > merged_fastqc_summaries_ovine_PBMC.txt 

# The summary file is processed in R - code in 03_for_bash_manipulate_results.R 

# FASTP ##################
########################### fastp version 0.19.7 

# Illumina Universal Adapter: AGATCGGAAGAG - it was found in the FastQC results 
# Check if I only need the universal adapter or the TruSeq adapter 
# Main quality issues - base quality bad at the tail of read 2, per tile sequence quality bad in some samples, adapter content bad in some samples 

# TruSeq primer
fastp --in1 /home/workspace/dniedziela/ovine_PBMC/fastq_files/A1_1.fq.gz --in2 /home/workspace/dniedziela/ovine_PBMC/fastq_files/A1_2.fq.gz --out1=/home/workspace/dniedziela/ovine_PBMC/fastp_results/A1_R1.clean_tur.fq.gz --out2=/home/workspace/dniedziela/ovine_PBMC/fastp_results/A1_R2.clean_tru.fq.gz -j /home/workspace/dniedziela/ovine_PBMC/fastp_results/A1_tru.json -h /home/workspace/dniedziela/ovine_PBMC/fastp_results/A1_tru.html -p -c --detect_adapter_for_pe --length_required 30 --qualified_quality_phred 20 --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT 
# Result html file has an insane number of overrepresented seqs for read1 and cuts off

# Universal primer only
fastp --in1 /home/workspace/dniedziela/ovine_PBMC/fastq_files/A1_1.fq.gz --in2 /home/workspace/dniedziela/ovine_PBMC/fastq_files/A1_2.fq.gz --out1=/home/workspace/dniedziela/ovine_PBMC/fastp_results/A1_R1.clean_uni.fq.gz --out2=/home/workspace/dniedziela/ovine_PBMC/fastp_results/A1_R2.clean_uni.fq.gz -j /home/workspace/dniedziela/ovine_PBMC/fastp_results/A1_uni.json -h /home/workspace/dniedziela/ovine_PBMC/fastp_results/A1_uni.html -p -c --detect_adapter_for_pe --length_required 30 --qualified_quality_phred 20 --adapter_sequence=AGATCGGAAGAG --adapter_sequence_r2=AGATCGGAAGAG 
# This option trimmed slightly more adapters, but trimmed sequenced are the exact same - so this is ok probably, overall adapters trimmed in approx. 0.2% reads  
# Overrepresented sequence analysis turned up sequences with 8 or 11 of those in the dataset - 0.001% or less of the raw reads, so I am turning off this function 

# options used -j -h -p -c 
# -c base correction for PE data enabled 
# -p enable overrepresented sequence analysis 
#important -g - can specify for polyG correction - it is a default for NovaSeq sequencing, because no signal in those systems means G 
#no default correction options were disabled 
# -j and -h specify reporting 
# -j, --json the json format report file name (string [=fastp.json])
# -h, --html the html format report file name (string [=fastp.html]) 
# minimum length left 30 bases 

# Make a script for all files - inside fastq_files folder
for i in $(ls *_1.fq.*);do  sn=`echo $i | cut -d '_' -f 1`; r2=`echo $i | sed 's/_1./_2./'`; echo  -e "fastp --in1 /home/workspace/dniedziela/ovine_PBMC/fastq_files/$i --in2 /home/workspace/dniedziela/ovine_PBMC/fastq_files/$r2 --out1=/home/workspace/dniedziela/ovine_PBMC/fastp_results/${sn}_R1.clean.fq.gz --out2=/home/workspace/dniedziela/ovine_PBMC/fastp_results/${sn}_R2.clean.fq.gz -j /home/workspace/dniedziela/ovine_PBMC/fastp_results/${sn}.json -h /home/workspace/dniedziela/ovine_PBMC/fastp_results/${sn}.html -c --detect_adapter_for_pe --length_required 30 --qualified_quality_phred 20 --adapter_sequence=AGATCGGAAGAG --adapter_sequence_r2=AGATCGGAAGAG" >> /home/workspace/dniedziela/ovine_PBMC/scripts/fastp.sh;done 

less /home/workspace/dniedziela/ovine_PBMC/scripts/fastp.sh # looks good 

nohup sh /home/workspace/dniedziela/ovine_PBMC/scripts/fastp.sh > /home/workspace/dniedziela/ovine_PBMC/fastp_results/fastp.nohup & 

# FASTP SUMMARY #################### 
####################################

for i in $(ls /home/workspace/dniedziela/ovine_PBMC/fastp_results/*.json);do echo  -e "echo ${i} >> ${i}_json.txt; grep "total_reads" $i >> ${i}_json.txt;grep "total_bases" $i >> ${i}_json.txt;grep "q20_rate" $i >> ${i}_json.txt;grep "q30_rate" $i >> ${i}_json.txt;grep "gc_content" $i  >> ${i}_json.txt;grep "read1_mean_length" $i >> ${i}_json.txt ; 
grep "read2_mean_length" $i >> ${i}_json.txt;
grep "passed_filter_reads" $i >> ${i}_json.txt;
grep "corrected_reads" $i >> ${i}_json.txt;
grep "corrected_bases" $i >> ${i}_json.txt;
grep "low_quality_reads" $i >> ${i}_json.txt;
grep "too_many_N_reads" $i >> ${i}_json.txt;
grep "too_short_reads" $i >> ${i}_json.txt;
grep "too_long_reads" $i >> ${i}_json.txt;
grep "q20_bases" $i >> ${i}_json.txt;
grep "q30_bases" $i >> ${i}_json.txt;
grep "adapter_trimmed_reads" $i >> ${i}_json.txt;
grep "adapter_trimmed_bases" $i >> ${i}_json.txt;">> /home/workspace/dniedziela/ovine_PBMC/scripts/json_summary2.sh;done 
# Added this grep "total_bases" $i >> ${i}_json.txt; to the script code above, not ran --- done 

less /home/workspace/dniedziela/ovine_PBMC/scripts/json_summary.sh
sh /home/workspace/dniedziela/ovine_PBMC/scripts/json_summary2.sh

mkdir json_summary2 ;
mv *_json.txt ./json_summary2 ;
cd json_summary2 
paste * > columned_json_summary2.txt #will merge files column by column 

# File was cleaned up in R - removed blank cells and text from where numbers are, and new values calculated - code in 03_for_bash_manipulate_results.R  

# FASTQC on trimmed sequences ############### 
#############################################

nohup fastqc -o /home/workspace/dniedziela/ovine_PBMC/fastqc_trimmed -t 20 /home/workspace/dniedziela/ovine_PBMC/fastp_results/*.clean.fq.gz > /home/workspace/dniedziela/ovine_PBMC/fastqc_trimmed/fastqc_trimmed.nohup & 

# FastQC summary again ## 

#Location of code: /home/workspace/dniedziela/ovine_PBMC/fastqc_trimmed
mkdir fastqc_summary_trimmed

for file in `ls *.zip`; \
do outfile=`basename $file | perl -p -e 's/\.zip//'`; \
echo -e "unzip /home/workspace/dniedziela/ovine_PBMC/fastqc_trimmed/$file; \
cp /home/workspace/dniedziela/ovine_PBMC/fastqc_trimmed/$outfile/summary.txt /home/workspace/dniedziela/ovine_PBMC/fastqc_trimmed/fastqc_summary_trimmed/${outfile}_summary.txt" >> /home/workspace/dniedziela/ovine_PBMC/scripts/fastqc_summary_trimmed.sh; \
done 

for file in `ls *.zip`; \
do outfile=`basename $file | perl -p -e 's/\.zip//'`; \
echo -e "cp /home/workspace/dniedziela/ovine_PBMC/fastqc_trimmed/$outfile/summary.txt /home/workspace/dniedziela/ovine_PBMC/fastqc_trimmed/fastqc_summary_trimmed/${outfile}_summary.txt" >> /home/workspace/dniedziela/ovine_PBMC/scripts/fastqc_summary_trimmed_copy.sh; \
done 

less /home/workspace/dniedziela/ovine_PBMC/scripts/fastqc_summary_trimmed_copy.sh
sh /home/workspace/dniedziela/ovine_PBMC/scripts/fastqc_summary_trimmed.sh #run shell script 
#after I do it I can use cat to combine all summaries into one file 
cd fastqc_summary_trimmed 
cat * > merged_fastqc_summaries_trimmed_ovine_PBMC.txt 

# STAR ##########################################
################################################# STAR version 2.7.3a 

# Download Ovis aries reference genome, Oar rambouillet version 1.0 from Ensembl
mkdir /home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_Ensembl/source_file

wget ftp://ftp.ensembl.org/pub/release-101/fasta/ovis_aries_rambouillet/dna/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.dna.toplevel.fa.gz 
gunzip Ovis_aries_rambouillet.Oar_rambouillet_v1.0.dna.toplevel.fa.gz

# Download annotation file 
mkdir /home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_Ensembl/annotation_file 

wget ftp://ftp.ensembl.org/pub/release-101/gtf/ovis_aries_rambouillet/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.101.gtf.gz
gunzip Ovis_aries_rambouillet.Oar_rambouillet_v1.0.101.gtf.gz

# Generate index with STAR 
# --sjdbOverhang you give "max read length-1" for your sequencing data 
mkdir /home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_Ensembl/STAR-2.7.3a_index

nohup STAR --runThreadN 40 --runMode genomeGenerate \
--genomeDir /home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_Ensembl/STAR-2.7.3a_index \
--genomeFastaFiles \
/home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_Ensembl/source_file/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.dna.toplevel.fa \
--sjdbGTFfile /home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_Ensembl/annotation_file/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.101.gtf \
--sjdbOverhang 149 --outFileNamePrefix \
/home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_Ensembl/STAR-2.7.3a_index &
# There was a sign which was not recognised as a new line for some reason, after 149 

# Generate script for all files to run STAR 
# One file 

STAR --readFilesCommand zcat --readFilesIn /home/workspace/dniedziela/ovine_PBMC/fastp_results/A1_R1.clean.fq.gz /home/workspace/dniedziela/ovine_PBMC/fastp_results/A1_R2.clean.fq.gz --runMode alignReads --genomeLoad LoadAndRemove --genomeDir /home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_Ensembl/STAR-2.7.3a_index --runThreadN 40 --outFileNamePrefix /home/workspace/dniedziela/ovine_PBMC/STAR_results/A1_ --outSAMmode Full --outReadsUnmapped Fastx --quantMode GeneCounts --limitBAMsortRAM 20000000000.0 --outSAMtype BAM SortedByCoordinate 

# BAMoutput.cpp:27:BAMoutput: exiting because of *OUTPUT FILE* error: could not create output file /home/workspace/dniedziela/ovine_PBMC/STAR_results/A1__STARtmp//BAMsort/19/26
# SOLUTION: check that the path exists and you have write permission for this file. Also check ulimit -n and increase it to allow more open files.
ulimit -n # 1024 
# Writing permissions are there for each file in the directories 
ulimit -n 10000 
# This worked, script ran, all files generated 

# Create a script for all files 
for i in $(ls *_R1.clean.fq.*);do  sn=`echo $i | cut -d '_' -f 1`; r2=`echo $i | sed 's/_R1./_R2./'`; echo  -e "STAR --readFilesCommand zcat --readFilesIn /home/workspace/dniedziela/ovine_PBMC/fastp_results/$i /home/workspace/dniedziela/ovine_PBMC/fastp_results/$r2 --runMode alignReads --genomeLoad LoadAndRemove --genomeDir /home/workspace/genomes/ovisaries/Oar_rambouillet_v1.0_Ensembl/STAR-2.7.3a_index --runThreadN 40 --outFileNamePrefix /home/workspace/dniedziela/ovine_PBMC/STAR_results/${sn}_ --outSAMmode Full --outReadsUnmapped Fastx --quantMode GeneCounts --limitBAMsortRAM 20000000000.0 --outSAMtype BAM SortedByCoordinate" >> /home/workspace/dniedziela/ovine_PBMC/scripts/STAR.sh;done 

less /home/workspace/dniedziela/ovine_PBMC/scripts/STAR.sh 
nohup sh /home/workspace/dniedziela/ovine_PBMC/scripts/STAR.sh & # Wooohooo it's running!!!! 


# Generate gene counts files from STAR output #######################
#####################################################################

grep "Uniquely mapped reads %" *_Log.final.out >> uniquely_mapped_ovine_PBMC.txt #this gives me a file with all sample names and uniquely mapped read % for all! 

# tidied up
mkdir BAM_files
mv *.bam ./BAM_files/
mkdir log_files
mv *.final.out ./log_files/
mkdir other_information
mv *.out ./other_information/
mkdir splice_junctions
mv *SJ.out.tab ./splice_junctions/ 
mv *mate* ./other_information/
mv nohup.out ../
# This only leaves gene counts files, uniquely mapped reads summary and the nohup file in the main directory 

mkdir STAR_ReadsPerGene 
cp ./*_ReadsPerGene.out.tab ./STAR_ReadsPerGene 

#get rid of the start (Unmapped reads etc) from each Reads per Gene file - only reads left now. 
for i in `ls *ReadsPerGene.out.tab`;do grep "ENSOARG" $i > ${i}.clean;done

#put in a header file for all files in the folder 
for i in `ls *clean`;do echo Processing ${i}; filename=${i%_ReadsPerGene.out.tab.clean}; echo -e "Genename\t${filename}_all\t${filename}_strand1\t$filename" > tmp;cat tmp $i > ${i}.2;mv ${i}.2 $i;  done
# The first column after gene name seems to be both strands, in the gene counts forward and reverse strand both have a lot of reads. 
# This is for sample A23, A1 and A2 don't have much in column 3. The library is stranded so I will take column 4, but there is a lot of antisense reads so I might analyse them separately too. 

# Cut the column with gene names from one file 
cut -f 1 A10_ReadsPerGene.out.tab.clean > A0_Genename_ReadsPerGene.out.tab.clean.cut ;
# Cut the gene counts column from all files 
for i in `ls *clean`;do cut -f 2 ${i} > ${i}_col2.cut ; done ; 
for i in `ls *clean`;do cut -f 3 ${i} > ${i}_col3.cut ; done ; 
for i in `ls *clean`;do cut -f 4 ${i} > ${i}_col4.cut ; done ; 

mkdir all_reads ;
mkdir strand1_reads ;
mkdir strand2_MY_reads ;
cp A0_Genename_ReadsPerGene.out.tab.clean.cut ./all_reads;
cp A0_Genename_ReadsPerGene.out.tab.clean.cut ./strand1_reads;
cp A0_Genename_ReadsPerGene.out.tab.clean.cut ./strand2_MY_reads ;
mv *_col2.cut ./all_reads ;
mv *_col3.cut ./strand1_reads ;
mv *_col4.cut ./strand2_MY_reads ;

# Combine Gene counts - in their respective folders 1
cd strand2_MY_reads ;
paste * >> ovine_PBMC_GeneCounts.txt ; 
cd .. ;
cd all_reads ;
paste * >> ovine_PBMC_GeneCounts_unstranded.txt ; 
cd .. ;
cd strand1_reads ;
paste * >> ovine_PBMC_GeneCounts_antisense.txt ; 
# END CODE 

# Indexes for IGV #################
################################### 

samtools index A10_Aligned.sortedByCoord.out.bam A10_Aligned.sortedByCoord.out.bam.bai
for i in ls *.bam; do echo "samtools index $i $i.bai" >> /home/workspace/dniedziela/ovine_PBMC/scripts/samtools_index.sh;done 

# I looked at one file, alignments colored and sorted by read strand. Alignment view squished. Data looks unstranded. 
# Group by first-in-pair seems to group by strand not read? 

# QuanTISeq ########## 

# Singularity version 3.7.1 installed 

wget -c https://icbi.i-med.ac.at/software/quantiseq/doc/downloads/quanTIseq_pipeline.sh 
chmod +x quanTIseq_pipeline.sh

# TPM file made in R 

bash ./software/quanTIseq_pipeline.sh --inputfile=./Quantiseq/tpms_for_quantiseq.tab --outputdir=./Quantiseq --pipelinestart="decon" 

nohup bash ./software/quanTIseq_pipeline.sh --inputfile=./Quantiseq/tpms_for_quantiseq.tab --outputdir=./Quantiseq --pipelinestart="decon" > quantiseq.nohup & 
