#!/bin/bash

##################################################################
# raw reads must be in folder 'raw_reads' in the $workdir folder #
# processed reads will be in the filtered_reads folder           #
##################################################################

# Go to work folder
workdir = "/ssd4tb"
cd $workdir

# Initial quality check using FastQC and MultiQC
echo "Preprocessing pipeline started" > $workdir/pipeline_log.txt
echo Current Date and Time is: `date +"%Y-%m-%d %T"` >> $workdir/pipeline_log.txt

SECONDS=0
echo "Running FastQC on raw reads" | tee -a $workdir/pipeline_log.txt

# run FastQC on all files in the dir trimmed_reads
find $workdir/raw_reads/ -name '*.fastq.gz' | xargs fastqc

# combine reports with MultiQC in the folder 'raw_reads'
multiqc $workdir/raw_reads/ -o $workdir/raw_reads/ -n raw_reads_report

echo "Time needed to finish FastQC step on raw reads: $SECONDS seconds" >> $workdir/pipeline_log.txt

# Removing erroneous k-mers from Illumina paired-end reads
# Run with the highest possible number of cores
SECONDS=0
echo "Running Rcorrector on raw reads" | tee -a $workdir/pipeline_log.txt
mkdir cor_reads
echo "Corrected reads will be in cor_reads folder" | tee -a $workdir/pipeline_log.txt

# Raw reads are in folder 'raw_reads'
fqdir=$workdir/raw_reads
for r1 in "$fqdir"/*1.fastq.gz; do
    r2=${r1%1.fastq.gz}2.fastq.gz
    if [[ -f $r2 ]] ; then
        perl $HOME/rcorrector/run_rcorrector.pl -1 $r1 -2 $r2 -t 24 -od '$workdir/cor_reads/'
    else
        echo "$r2 not found" >&2
    fi
done
echo "Time needed to finish Rcorrector step on raw reads: $SECONDS seconds" >> $workdir/pipeline_log.txt

# Discard read pairs for which one of the reads is deemed unfixable
# Original python script from https://github.com/harvardinformatics/TranscriptomeAssemblyTools
# Translated to Python3
# Not memory intensive, therefore assign as many cores as possible with -P => parallel processes

SECONDS=0
echo "Running Filter Uncorrectable on corrected reads" | tee -a $workdir/pipeline_log.txt
echo "Filtered reads will be in folder clean_reads" | tee -a $workdir/pipeline_log.txt
mkdir clean_reads
cd clean_reads
ls $workdir/cor_reads/*1.cor.fq.gz | xargs -P30 -I@ bash -c 'python $HOME/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq_P3.py -1 "$1" -2 ${1%1.*.*}2.cor.fq.gz' _ @
echo "Time needed to finish Rcorrector step on raw reads: $SECONDS seconds" >> $workdir/pipeline_log.txt
cd

SECONDS=0
echo "Running Trim Galore" | tee -a pipeline_log.txt
echo Current Date and Time is: `date +"%Y-%m-%d %T"` >> $workdir/pipeline_log.txt

# adapter removal and read quality trimming of paired-read fastq-files
# trim_galore --paired --retain_unpaired --phred33 --output_dir trimmed_reads --length 50 -q 5 --stringency 1 -e 0.1 --cores 2 SAMPLE_R1.fastq.gz SAMPLE_R2.fastq.gz

# create dir for trimmed reads
mkdir trimmed_reads

# run on all pair-read fastq.gz files in a folder
# number of core fixed to 8! Read cutadapt manual for explanation

fqdir=$workdir/clean_reads
for r1 in "$fqdir"/*1.cor.fq; do
    r2=${r1%1.cor.fq}2.cor.fq
    if [[ -f $r2 ]] ; then
trim_galore --paired --retain_unpaired --phred33 \
--output_dir $workdir/trimmed_reads --length 50 -q 5 \
--stringency 1 -e 0.1 --cores 8 "$r1" "$r2"
    else
        echo "$r2 not found" >&2
    fi
done

echo "Time needed to finish Trim Galor step: $SECONDS seconds" >> $workdir/pipeline_log.txt

# Check against SILVA rRNA db
# The files you want are *_blacklist_paired_unaligned.fq.gz

fqdir=$workdir/trimmed_reads
for r1 in "$fqdir"/*1.cor_val_1.fq; do
    r2=${r1%1.cor_val_1.fq}2.cor_val_2.fq
    if [[ -f $r2 ]] ; then
        bowtie2 --quiet --very-sensitive-local \
--phred33  -x $HOME/SILVA/SILVA.db -1 "$r1" -2 "$r2" --threads 6 \
--met-file ${r1%.fq}_bowtie2_metrics.txt \
--al-conc-gz ${r1%.fq}_blacklist_paired_aligned.fq.gz \
--un-conc-gz ${r1%.fq}_blacklist_paired_unaligned.fq.gz  \
--al-gz ${r1%.fq}_blacklist_unpaired_aligned.fq.gz \
--un-gz ${r1%.fq}_blacklist_unpaired_unaligned.fq.gz
    else
        echo "$r2 not found" >&2
    fi
done

SECONDS=0
echo "Running FastQC on processed reads" | tee -a $workdir/pipeline_log.txt
echo Current Date and Time is: `date +"%Y-%m-%d %T"` >> $workdir/pipeline_log.txt

# Moving the files we want to filtered_reads folder

mkdir $workdir/filtered_reads
mv $workdir/trimmed_reads/*_paired_unaligned.fq* $workdir/filtered_reads
cd $workdir/filtered_reads

# Postprocessing quality check
# run FastQC on all files in the trimmed_reads folder

find -name '*_paired_unaligned.fq*.gz' | xargs fastqc

# Move Trim Galore and FastQC reports to endQC folder for MultiQC
mv $workdir/trimmed_reads/*.txt $workdir/endQC
mv $workdir/filtered_reads/*.html $workdir/endQC

# combine reports with MultiQC
multiqc $workdir/endQC/ -o $workdir/endQC/ -n filtered_reads_report

echo "Time needed to finish second FastQC step on processed reads: $SECONDS seconds" | tee -a $workdir/pipeline_log.txt
echo "FastQC on processed reads finished." | tee -a $workdir/pipeline_log.txt
echo "###################################" | tee -a $workdir/pipeline_log.txt
echo "Preprocessing pipeline finished." | tee -a $workdir/pipeline_log.txt
cd
