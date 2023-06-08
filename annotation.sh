# Preparing new gene_trans_map for cdhit90.fasta

## Extracting sequence names
awk 'sub(/^>/, "")' '/ssd4tb/Trinity_cdhit90.fasta' \
> '/ssd4tb/Trinity_cdhit90_headers.txt'

## Removing the unwanted part
awk '{print $1}' '/ssd4tb/Trinity_cdhit90_headers.txt' \
> Trinity_cdhit90_header_filtered.txt

## Filtering the original gene_trans_map
grep -Fwf '/ssd4tb/Trinity_cdhit90_header_filtered.txt' \
'/ssd4tb/trinity_out_dir.Trinity.fasta.gene_trans_map' \
> Trinity_cdhit90.fasta.gene_trans_map

## Check if both files have the same number of rows
grep ">" '/ssd4tb/Trinity_cdhit90.fasta' | wc -l
wc -l Trinity_cdhit90.fasta.gene_trans_map

# Run TransDecoder
TransDecoder.LongOrfs -t '/ssd4tb/Trinity_cdhit90.fasta' -S
TransDecoder.Predict -t '/ssd4tb/Trinity_cdhit90.fasta'

# Run Trinotate
## Create Trinotate db
Trinotate --create \
 --db Np_complete.sqlite \
 --trinotate_data_dir '/ssd4tb/Np_complete' \
 --use_diamond

## Initialize Trinotate db
Trinotate --db Np_cdhit90.sqlite --init \
           --gene_trans_map '/ssd4tb/Trinity_cdhit90.fasta.gene_trans_map' \
           --transcript_fasta '/ssd4tb/Trinity_cdhit90.fasta' \
           --transdecoder_pep '/ssd4tb/Trinity_cdhit90.fasta.transdecoder.pep'

## Annotate contigs
Trinotate --db <sqlite.db> --CPU 30 \
               --transcript_fasta '/ssd4tb/Trinity_cdhit90.fasta' \
               --transdecoder_pep '/ssd4tb/Trinity_cdhit90.fasta.transdecoder.pep' \
               --trinotate_data_dir '/ssd4tb/Np_cdhit90'
               --run "swissprot_blastp swissprot_blastx pfam signalp6 tmhmmv2 infernal EggnogMapper" \
               --use_diamond

## Generate Trinotate report
Trinotate --db Np_cdhit90.sqlite --report > Np_cdhit90.tsv
