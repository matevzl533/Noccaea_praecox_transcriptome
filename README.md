# Transcriptome de novo assembly protocol
Protocol roughly follows the protocol described at: https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html and uses Trinity for de novo assembly. Annotations of the contigs were done using Trinotate. Differential expression analysis was done using DESeq2, edgeR and limma/voom in R.

## Preprocessing

### Initial quality control

We use FastQC with MultiQC to assess the quality of raw reads and visualise the reports.

### Removing erroneous k-mers from Illumina paired-end reads

We use [rCorrector](https://github.com/mourisl/Rcorrector), a tool specifically designed for RNA-seq data. Besides being a top performer in side-by-side comparisons of kmer-based read error correction tools, it tags reads in the fastq output as corrected, or uncorrectable.

### Discard read pairs for which one of the reads is deemed unfixable

A python script to remove any read pair where at least one read has an unfixable error, can be obtained from the Harvard Informatics GitHub repository [TranscriptomeAssemblyTools](https://github.com/harvardinformatics/TranscriptomeAssemblyTools).

### Trim adapter and low quality bases from fastq files

We used [Trim Galore](https://github.com/FelixKrueger/TrimGalore) for adapter and quality trimming. Trim Galore is a wrapper around [Cutadapt](https://github.com/marcelm/cutadapt) and [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

### Map trimmed reads to a blacklist to remove unwanted (rRNA reads)

From SILVA, we download the SSUParc and LSUParc fasta files (https://ftp.arb-silva.de/?pk_vid=8352a8ccf0ead1d7168388545541b6c1), concatenating them, and replacing U characters with T, as our sequence reads are in DNA space.

### Run fastqc on your processed reads that pass qc and filtering from the above steps

Re-run the QC from step 1. Ideally, there will be no trend in adapter contamination by cycle, and there will be increased evenness in kmer distributions, GC content, no over-represented sequences, etc.

## Transcriptome assembly

### de novo assembly with Trinity

We used Trinity for transcriptome de novo assembly with default parameters with the addition of strand information (--SS_lib_type RF) and minimum contig length (--min_contig_length 200).

### Removing redundancy

To reduce redundancy, contigs from the assembly were clustered using CD-HIT with the following parameters: -s 0.9 -aS 0.9.

## Accessing assembly quality

### TrinityStats.pl

N50 statistics, and counts of the number of Trinity contigs was generated using the TrinityStats.pl script that comes with Trinity.

### Quantify read support for the assembly

To evaluate read support for the assembly is a three step process. First, you build a bowtie2 index for your assembly. Next, you map your reads and calculate alignment statistics. The align_stats.txt file will provide info on the percentage of read pairs that mapped concordantly, as well as an overall alignment rate.

### Quantifying completeness

To assess completeness, we use BUSCO. 
