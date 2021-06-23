##################################
# Ovine PBMC Liver Fluke RNA-seq #
##################################

# Author: Carolina N. Correia
# Last updated on: 16/06/2020

###################################
# Download raw data from Novogene #
###################################

# Create and enter working directory on Rodeo:
mkdir -p /home/workspace/ccorreia/ovinePBMC_RNAseq/fastq
cd !$

# Download batch 1 files:
screen -D -R ovine_batchF001_1
wget -c https://s3.eu-west-1.amazonaws.com/novogene-europe/HW/project/X204SC20032696-Z01-F001_1_20200612_LkMpsi.zip
# Detach the screen session by pressing Ctrl+A then d

# Download batch 2 files:
screen -D -R ovine_batchF001_2
wget -c https://s3.eu-west-1.amazonaws.com/novogene-europe/HW/project/X204SC20032696-Z01-F001_2_20200612_ZVI2CP.zip
# Detach the screen session by pressing Ctrl+A then d

# Download batch 3 files:
screen -D -R ovine_batchF001_3
wget -c https://s3.eu-west-1.amazonaws.com/novogene-europe/HW/project/X204SC20032696-Z01-F001_3_20200612_Xf9o1r.zip
# Detach the screen session by pressing Ctrl+A then d

# Check if download is complete:
cd /home/workspace/ccorreia/ovinePBMC_RNAseq/fastq
screen -D -R ovine_batchF001_1 # Detach the screen session by pressing Ctrl+A then d.
screen -D -R ovine_batchF001_2 # Detach the screen session by pressing Ctrl+A then d.
screen -D -R ovine_batchF001_3 # Detach the screen session by pressing Ctrl+A then d.

# Check if the downloading jobs are still active (screen sessions will still be
# active after the download is complete):
htop

# After you have confirmed all downloads have been completed,
# the screens can be terminated:
screen -X -S ovine_batchF001_1 kill
screen -X -S ovine_batchF001_2 kill
screen -X -S ovine_batchF001_3 kill

# Check that all screen sessions have been terminated:
htop

########################
# Perform MD5 checksum #
########################

# Enter working directory:
cd /home/workspace/ccorreia/ovinePBMC_RNAseq/fastq

# Create a file with the md5 hashes provided by Novogene via email:
printf "2ff82e3a56b512e11d70121c7f78a548  X204SC20032696-Z01-F001_1_20200612_LkMpsi.zip
02bcf3c083e0f284a84a2a28fbe32ad4  X204SC20032696-Z01-F001_2_20200612_ZVI2CP.zip
6741c7fc29710ada71e76cb1c85a69b6  X204SC20032696-Z01-F001_3_20200612_Xf9o1r.zip" \
> novogene.md5

# Perform md5sum check for .zip files:
md5sum -c novogene.md5 >> \
/home/workspace/ccorreia/ovinePBMC_RNAseq/fastq/md5check_UCD.txt

# Check that all .zip files passed the check:
grep -c 'OK' md5check_UCD.txt

# If all files passed the md5 check, unzip them:
unzip X204SC20032696-Z01-F001_1_20200612_LkMpsi.zip
unzip X204SC20032696-Z01-F001_2_20200612_ZVI2CP.zip
unzip X204SC20032696-Z01-F001_3_20200612_Xf9o1r.zip

# Modify .fq.gz file permissions to read and execute only:
for file in `find /home/workspace/ccorreia/ovinePBMC_RNAseq/fastq/ \
-name *.fq.gz`; \
do chmod -R 555 $file; \
done

# Create a bash script to perform md5sum check for all .fq.gz files:
for file in `find /home/workspace/ccorreia/ovinePBMC_RNAseq/fastq/ \
-name MD5.txt`; \
do folder=`dirname $file`; \
echo "cd $folder; \
md5sum -c $file >> ${folder}/md5check_UCD.txt" >> check_md5.sh; \
done

# Run script on Rodeo:
chmod 755 check_md5.sh
nohup ./check_md5.sh &

# Check that all .fq.gz files passed the check:
for file in `find /home/workspace/ccorreia/ovinePBMC_RNAseq/fastq/ \
-name md5check_UCD.txt`; \
do grep 'OK' $file >> passed.txt; \
done

wc -l passed.txt # 99 files (96 .fq.gz + 3 .zip)


### only delete .zip folders after grace backs the data up in an external hard drive
# If all .fq.gz files passed the md5 check, delete the .zip folders
cd /home/workspace/ccorreia/ovinePBMC_RNAseq/fastq
rm -r X204SC20032696*.zip



