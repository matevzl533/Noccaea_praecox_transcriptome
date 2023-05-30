library(trinotateR)
library(dplyr)
library(ggplot2)

# Modified split_GO functions due to Trinotate table changes
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

# Read Trinotate table
trinotate_data <- read_trinotate("data/Np_cdhit90.tsv")

# Generate GO object
go <- split_GOp(trinotate_data)

# Count and remove terms with less than 1000 counts
go_1000 <- go %>%
  count(ontology, name, sort=TRUE) %>%
  filter(n > 1000) %>%
  arrange(ontology)

# Generate a vector to override ggplot alphabetical ordering
go_1000$name <- factor(go_1000$name, levels = go_1000$name)

# Generate ggplot object
go_plot.gg <- go_1000 %>%
  ggplot(.) +
  geom_bar(aes(x=name, y=n, fill=ontology), stat = "identity") +
  theme_classic() +
  ylab("No. of transcripts") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank())

# Save ggplot object
ggsave(go_plot.gg, file="plots/GO_terms_plot.pdf", dpi = 300, width = 3000, height = 1500, units = "px")
