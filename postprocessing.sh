#!/bin/bash

# Reads from the preprocessing pipeline are in filtered_reads folder
# Trinity assembly is in $workdir folder
# Trinity tools are in '$HOME/trinityrnaseq-master/util' folder
# Define $pattern for automatic listing of samples
# Sample_files.txt is sample file used with Trinity
# Requires dplyr, data.table and stringr (v1.5.0) packages installed in R

pattern="Np_"
workdir="/ssd4tb"
export PATH="$PATH:$HOME/trinityrnaseq-master/util/"

################################################
# Removal of redundancy and quality estimation #
################################################

# Removing redundancy
echo "Removing redundancy from the assembly" | tee $workdir/postprocess_log.txt
cd-hit-est \
	-i trinity_out_dir.Trinity.fasta \
	-o Trinity_cdhit90.fasta -c 0.90 -n 9 -d 0 -M 0 -T 30 -s 0.9 -aS 0.9

echo "Redundancy removed, proceeding to basic statistics" | tee -a $workdir/postprocess_log.txt

# Basic statistics
TrinityStats.pl Trinity_cdhit90.fasta > Trinity_cdhit90_stats.txt

# Mapping reads to the assembly
## Creating bowtie2 db
bowtie2-build '$workdir/Trinity_cdhit90.fasta' \
	Trinity.db --threads 30

## Mapping the reads
bowtie2 -p 12 -q --no-unal -k 20 \
-x '$workdir/Trinity.db' \
-1 '$workdir/trinity_out_dir/insilico_read_normalization/left.norm.fq' \
-2 '$workdir/trinity_out_dir/insilico_read_normalization/right.norm.fq' \
2> '$workdir/align_stats.txt' | samtools view -@10 -Sb -o bowtie2.bam

# Preparing new gene_trans_map for cdhit90.fasta
# Extracting sequence names
awk 'sub(/^>/, "")' $workdir/Trinity_cdhit90.fasta \
> Trinity_cdhit90_headers.txt

# Removing the unwanted part
awk '{print $1}' $workdir/Trinity_cdhit90_headers.txt \
> Trinity_cdhit90_header_filtered.txt

# Filtering the original gene_trans_map
# Use -v to inverse the selection
grep -Fwf $workdir/Trinity_cdhit90_header_filtered.txt \
$workdir/trinity_out_dir.Trinity.fasta.gene_trans_map \
> Trinity_cdhit90.fasta.gene_trans_map

# Kallisto
align_and_estimate_abundance.pl \
--transcripts $workdir/Trinity_cdhit90.fasta \
--gene_trans_map $workdir/Trinity_cdhit90.fasta.gene_trans_map \
--seqType fq --samples_file $workdir/sample_files_compl.txt \
--est_method kallisto --aln_method bowtie2 --SS_lib_type RF \
--thread_count 24 --trinity_mode --prep_reference --output_dir $workdir/kallisto_outdir

# for Kallisto use this
find $pattern* -name "abundance.tsv" | tee quant_files.list

abundance_estimates_to_matrix.pl \
--est_method kallisto \
--gene_trans_map $workdir/Trinity_cdhit90.fasta.gene_trans_map \
--quant_files $workdir/quant_files.list --name_sample_by_basedir

# ExN50
export PATH="$PATH:$HOME/trinityrnaseq-master/util/misc"

contig_ExN50_statistic.pl $workdir/kallisto.isoform.TMM.EXPR.matrix \
$workdir/Trinity_cdhit90.fasta transcript | tee ExN50.transcript.stats

plot_ExN50_statistic.Rscript $workdir/ExN50.transcript.stats

xpdf $workdir/ExN50_plot.pdf

echo "Basic statistics finished. All output files generated." | tee -a $workdir/postprocess_log.txt
echo "Starting annotations using Trinotate" | tee -a $workdir/postprocess_log.txt
