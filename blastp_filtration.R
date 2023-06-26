#!/usr/bin/Rscript

# require the needed packages
require("dplyr")
require("data.table")

# Set work directory
setwd("/ssd4tb")

# An original function from trinotateR (https://github.com/cstubben/trinotateR)
read_trinotate <- function(file, ...){
   ## suppress fread warning messages like Bumped column 5 to type character 
   x <- suppressWarnings( fread(file, sep="\t", na.strings=".", ...) )
   # message("Read ", nrow(x), " rows")
   # replace "#gene_id"  since names in row 1 start with # 
   names(x)[1] <- "gene_id"
   setkeyv(x, c("gene_id", "transcript_id") )
   x
}

# A modified function from trinotateR due to changes in the Trinotate output table
split_blastp <- function(x, hit = "sprot_Top_BLASTP_hit"){
  y <- x[!is.na( get(hit) ), .( get(hit), gene_id, transcript_id, prot_id) ]
  # split multiple annotations in backtick-delimited list
  z <- strsplit(y$V1, "`" )
  n <- sapply(z, length)
  # and then split BLAST annotation into 7 columns
  
  z <- strsplit( unlist(z), "\\^" )
  # first two columns are identical 
  if(  any(sapply(z, "[", 1) != sapply(z, "[", 2)  ) ) print("WARNING: check different values in columns 1 and 2")
  
  NAME <- gsub("^RecName: Full=", "", sapply(z, "[", 6) )
  NAME <- gsub("SubName: Full=", "", NAME )
  NAME <- gsub(";$", "", NAME )
  #   drop {ECO:0000256|RuleBase:RU361189} from names
  NAME <- gsub(" \\{[^}]+}", "", NAME)
  
  x1 <- data.frame(  gene = rep(y$gene_id, n), 
                     transcript = rep(y$transcript_id, n) ,  
                     protein = rep(gsub(".*\\|", "", y$prot_id), n),
                     uniprot = sapply(z, "[", 1) , 
                     align = sapply(z, "[", 3) , 
                     identity = as.numeric( gsub( "%ID", "", sapply(z, "[", 4))) , 
                     evalue = as.numeric(gsub("E:", "", sapply(z, "[", 5) )), 
                     name = NAME, 
                     lineage = sapply(z, "[", 7) , 
                     domain = gsub( "; .*", "", sapply(z, "[", 7) ) , 
                     genus = gsub( ".*; ", "", sapply(z, "[", 7) ) ,
                     stringsAsFactors=FALSE )
  message(nrow(x1), " ", hit, " annotations")
  data.table(x1)
}

# An original function from trinotateR (https://github.com/cstubben/trinotateR)
summary_blast <- function(x){
   uniq <-  c( UniProt = uniqueN(x$uniprot),  genes=  uniqueN(x$gene), transcripts = uniqueN(x$transcript), proteins = uniqueN(x$protein) )
   y <- x[,.( genes=uniqueN(gene), transcripts=uniqueN(transcript), proteins = uniqueN(protein) , total =.N ), by=.(uniprot, domain, genus, name) ]
 
   y <- y[order(-y$genes, tolower(y$name)),]

    annot <- colSums(y[, .(total, genes, transcripts, proteins) ])    
   z <-  cbind(unique=uniq, annotations= annot) 
   attr(y, "counts") <- z
   message(nrow(y), " rows")
   y
}

# Import the original Trinotate table
trinotate_data <- read_trinotate("Trinity_cdhit90.tsv")

# Split BLASTp results
blastp <- split_blastp(trinotate_data)

# Import qcovs not present in trinotate table
blast_add <- read.delim("Trinity_cdhit90_blastp_additional.txt", row.names = NULL, col.names=c("transcript","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qlen","slen","qcovs","qcohsp")) %>%
  mutate(transcript = stringr::str_split_i(transcript, "::", 2)) %>%
  select(transcript, qcovs)

# Join both tables 
blastp <- left_join(blastp, blast_add, by = "transcript")

# Apply the filters for removal of contaminants and create a filter object
blastp_filter <- blastp %>%
  # if hits fall into Primates, Muroidea or Bovidea
  mutate(contamination = ifelse(grepl("Muroidea;", lineage),"cont",ifelse(grepl("Bovidea", lineage),"cont", ifelse(grepl("Primates", lineage),"cont", ifelse(grepl("Fungi", lineage), "cont", "/")))), "/") %>%
  filter(contamination == "cont") %>%
  # if hits have low enough e-value
  filter(evalue < 1e-20) %>%
  # if hits have high enough identity
  filter(identity > 90) %>%
  # if hits have high enough cover
  filter(qcovs > 40) %>%
  select(transcript)

# Write blastp_filter to a file for filtration of the original fasta
write.table(blastp_filter, file="blastp_filter.txt", sep="\t", col.names = F, row.names = F, quote = FALSE)

# Filter the original trinotate table
trinotate_data_filtered <- trinotate_data %>%
  filter(., !transcript_id %in% blastp_filter$transcript) 

# Save filtered trinotate table as csv
write.csv(trinotate_data_filtered, file="trinotate_table_filtered.csv")
