# Corrected original functions from trinotateR library due to changes in Trinotate tables and additional needs during the analysis

split_GOp <- function(x, hit = "gene_ontology_BLASTP"){
  y <- x[!is.na( get(hit) ), .( get(hit), gene_id, transcript_id, prot_id) ]
  
  # split multiple annotations in backtick-delimited list
  z <- strsplit(y$V1, "`" )
  n <- sapply(z, length)
  # and then split GO annotation into 3 columns
  
  z <- strsplit( unlist(z), "\\^" )
  
  x1 <- data.frame(  gene = rep(y$gene_id, n), 
                     transcript = rep(y$transcript_id, n) ,  
                     protein = rep(gsub(".*\\|", "", y$prot_id), n),
                     go = sapply(z, "[", 1) , 
                     ontology = sapply(z, "[", 2) , 
                     name = sapply(z, "[", 3) , 
                     stringsAsFactors=FALSE )
  message(nrow(x1), " ", hit, " annotations")
  data.table(x1)
}

split_GOx <- function(x, hit = "gene_ontology_BLASTX"){
  y <- x[!is.na( get(hit) ), .( get(hit), gene_id, transcript_id, prot_id) ]
  
  # split multiple annotations in backtick-delimited list
  z <- strsplit(y$V1, "`" )
  n <- sapply(z, length)
  # and then split GO annotation into 3 columns
  
  z <- strsplit( unlist(z), "\\^" )
  
  x1 <- data.frame(  gene = rep(y$gene_id, n), 
                     transcript = rep(y$transcript_id, n) ,  
                     protein = rep(gsub(".*\\|", "", y$prot_id), n),
                     go = sapply(z, "[", 1) , 
                     ontology = sapply(z, "[", 2) , 
                     name = sapply(z, "[", 3) , 
                     stringsAsFactors=FALSE )
  message(nrow(x1), " ", hit, " annotations")
  data.table(x1)
}
split_blastx <- function(x, hit = "sprot_Top_BLASTX_hit"){
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
summary_blast_exp <- function(x){
  uniq <-  c( UniProt = uniqueN(x$uniprot),  genes=  uniqueN(x$gene), transcripts = uniqueN(x$transcript), proteins = uniqueN(x$protein) )
  y <- x[,.( genes=uniqueN(gene), transcripts=uniqueN(transcript), proteins = uniqueN(protein) , total =.N ), by=.(uniprot, domain, genus, name, tax) ]
  
  y <- y[order(-y$genes, tolower(y$name)),]
  
  annot <- colSums(y[, .(total, genes, transcripts, proteins) ])    
  z <-  cbind(unique=uniq, annotations= annot) 
  attr(y, "counts") <- z
  message(nrow(y), " rows")
  y
}
