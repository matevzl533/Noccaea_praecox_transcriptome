#!/usr/bin/python

# Script from: https://bioinformaticsworkbook.org/dataWrangling/fastaq-manipulations/calculate-sequence-lengths-in-a-fasta-file.html#gsc.tab=0
# Run with seq_length.py file.fasta > stats.txt

from Bio import SeqIO
import sys
cmdargs = str(sys.argv)
for seq_record in SeqIO.parse(str(sys.argv[1]), "fasta"):
 output_line = '%s\t%i' % \
(seq_record.id, len(seq_record))
 print(output_line)
