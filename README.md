# Noccaea praecox transcriptome de novo assembly protocol
Protocol roughly follows the protocol described at: https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html and uses Trinity for de novo assembly. Annotations of the contigs were done using Trinotate.

## Preprocessing

The codes for the preprocessing steps are available in the preprocessing.sh file.

### Initial quality control

We use FastQC with MultiQC to assess the quality of raw reads and visualise the reports.

### Removing erroneous k-mers from Illumina paired-end reads

We use [rCorrector](https://github.com/mourisl/Rcorrector), a tool specifically designed for RNA-seq data. Besides being a top performer in side-by-side comparisons of kmer-based read error correction tools, it tags reads in the fastq output as corrected or uncorrectable.

### Discard read pairs for which one of the reads is deemed unfixable

A python script to remove any read pair where at least one read has an unfixable error can be obtained from the Harvard Informatics GitHub repository [TranscriptomeAssemblyTools](https://github.com/harvardinformatics/TranscriptomeAssemblyTools).

### Trim adapter and low-quality bases from fastq files

We used [Trim Galore](https://github.com/FelixKrueger/TrimGalore) for adapter and quality trimming. Trim Galore is a wrapper around [Cutadapt](https://github.com/marcelm/cutadapt) and [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

### Map trimmed reads to a blacklist to remove unwanted (rRNA reads)

From SILVA, we download the SSUParc and LSUParc fasta files (https://ftp.arb-silva.de/?pk_vid=8352a8ccf0ead1d7168388545541b6c1), concatenate them, and replacing U characters with T, as our sequence reads are in DNA space. This must be done before running samples, as it is not included in the main code.

### Run fastqc on your processed reads that pass qc and filtering from the above steps

Re-run the QC from step 1. Ideally, there will be no trend in adapter contamination by cycle, and there will be increased evenness in kmer distributions, GC content, no over-represented sequences, etc.

## Transcriptome assembly

### de novo assembly with Trinity

We used [Trinity](https://github.com/trinityrnaseq/trinityrnaseq) for transcriptome de novo assembly with default parameters with the addition of strand information (--SS_lib_type RF).
The codes for the steps after assembly with Trinity (redundancy removal and evaluation of assembly quality) are included in the postprocessing.sh file.

### Removing redundancy

To reduce redundancy, contigs from the assembly were clustered using [CD-HIT](https://github.com/weizhongli/cdhit) with the following parameters: -s 0.9 -aS 0.9.

## Accessing assembly quality

### TrinityStats.pl

N50 statistics and counts of the number of Trinity contigs were generated using the TrinityStats.pl script that comes with Trinity.

### Quantify read support for the assembly

To evaluate read support for the assembly, we built a bowtie2 index for the assembly, mapped the reads and calculated the alignment statistics. The align_stats.txt file provides info on the percentage of read pairs that mapped concordantly, as well as an overall alignment rate.

### Quantifying completeness

To assess completeness, we used BUSCO.

## Annotations

The assembled transcriptome was annotated using the TransDecoder (Trinity) and [Trinotate](https://github.com/Trinotate/Trinotate/wiki) as described in the wikis of both tools.

## Filtration of the contaminants

To find contigs originating outside of the N. praecox transcriptome, we used the [NCBI Foreign Contamination Screen (FCS) caller](https://github.com/ncbi/fcs)
