library(trinotateR)
library(dplyr)
library(ggplot2)

# Modified split_GO functions from trinotateR due to Trinotate table changes
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

##################################################################
# Figure 1
##################################################################

# Load length distribution data
trinity_length <- read.table("data/Trinity_length_distribution.txt")
trinity_cdhit90_length <- read.table("data/Trinity_cdhit90_length_distribution.txt")

# Separate into bins and plot
## Original assembly
trinity_ld.gg <- trinity_length %>%
  group_by(group = cut(V2, breaks = seq(0, 3000, 200))) %>%
  summarise(n = sum(V2)) %>%
  ggplot(.) +
  geom_bar(aes(x=group, y=n), stat = "identity", fill="#F8766D") +
  theme_classic() +
  ylab("No. of transcripts") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank(),
        legend.position = "none")

## Reduced redundancy assembly
trinity_cdhit90_ld.gg <- trinity_cdhit90_length %>%
  group_by(group = cut(V2, breaks = seq(0, 3000, 200))) %>%
  summarise(n = sum(V2)) %>%
  ggplot(.) +
  geom_bar(aes(x=group, y=n), stat = "identity", fill="#7CAE00") +
  theme_classic() +
  ylab("No. of unigenes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank(),
        legend.position = "none")

# Arrange plots into a grouped plot and save
length_distributions.gg <- ggpubr::ggarrange(trinity_ld.gg, trinity_cdhit90_ld.gg, ncol=2, nrow = 1, labels = c("a", "b"))
ggsave(length_distributions.gg, file="plots/length_distributions.png", dpi = 300, width = 3000, height = 1500, units = "px")


##################################################################
# Figure 2
##################################################################

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

##################################################################
# Figure 3
##################################################################

library(limma)

# Get Kegg IDs
kegg_genes <- getGeneKEGGLinks(species.KEGG = "ath", convert = FALSE)
kegg_links <- getKEGGPathwayNames(species.KEGG = "ath", remove.qualifier = FALSE)

# Join the Kegg IDs and gene IDs
kegg <- ath %>%
  select(Kegg) %>%
  rename(GeneID = Kegg) %>%
  left_join(., kegg_genes, by="GeneID", multiple="all") %>%
  mutate(PathwayID = stringr::str_remove(PathwayID, "path:")) %>%
  left_join(., kegg_links, by="PathwayID") %>%
  mutate(Description = stringr::str_split_i(Description, " - ",1)) %>%
  group_by(Description) %>%
  summarise(Count = n()) %>%
  arrange(., desc(Count))
  
# Filter for Kegg terms with more than 500 transcripts
kegg_500 <- kegg %>%
  filter(Count > 500)
  
# Generate a vector to override ggplot alphabetical ordering
kegg_500$Description <- factor(kegg_500$Description, levels = kegg_500$Description)

# Generate ggplot object
kegg_500.gg <- kegg_500 %>%
  filter(!is.na(Description)) %>%
  ggplot(.) +
  geom_bar(aes(x=Description, y=Count), stat = "identity") +
  theme_classic() +
  ylab("No. of transcripts") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank())
  
# Save ggplot object
ggsave(kegg_plot.gg, file="plots/kegg_terms_plot.pdf", dpi = 300, width = 3000, height = 1500, units = "px")  
