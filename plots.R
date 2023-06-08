##########################################################
# Folders structure                                      #
# data...folder with input files (eg. trinotate results  #
# plots...folder for graphical results                   #
# tables...folder for text results                       #
# customR...folder for external sources with R functions #
##########################################################

#devtools::install_github("cstubben/trinotateR")
library(trinotateR)
library(dplyr)
library(ggplot2)
library(palmerpenguins)
library(ggpubr)
library(limma)
source("customR/trinotateR_add.R")
source("customR/penguin_palette.R")

# Import Trinotate table
trinotate_data <- read_trinotate("data/Np_cdhit90.tsv")
# Import contigs' length distribution
trinity_cdhit90_length <- read.table("data/Trinity_cdhit90_length_distribution.txt", col.names = c("transcript", "length"))
# Import BUSCO results
busco <- read.table("data/BUSCO.txt", header = T)

# Plot Figure 1a
busco_transcripts.gg <- busco %>%
  ggplot(.,aes(x="", y=Transcripts, fill=BUSCOs)) +
  geom_bar(stat = "identity", width=1) +
  coord_polar("y", start=0, direction = 1) +
  theme_void() +
  ylab("No. of transcripts") +
  #labs(fill="Top 10 KEGG terms") +
  scale_fill_penguin(palette = "middle-main") +
  theme(axis.text = element_blank())

# Save ggplot object
ggsave(busco_transcripts.gg, file="plots/busco_transcripts_plot.png", dpi = 300, width = 3000, height = 1500, units = "px")  

# Plot Figure 1b
busco_proteins.gg <- busco %>%
  ggplot(.,aes(x="", y=Proteins, fill=BUSCOs)) +
  geom_bar(stat = "identity", width=1) +
  coord_polar("y", start=0, direction = 1) +
  theme_void() +
  ylab("No. of transcripts") +
  #labs(fill="Top 10 KEGG terms") +
  scale_fill_penguin(palette = "middle-main") +
  theme(axis.text = element_blank())

# Save ggplot object
ggsave(busco_proteins.gg, file="plots/busco_proteins_plot.png", dpi = 300, width = 3000, height = 1500, units = "px")  

# Complete Figure 1
busco_grouped.gg <- ggpubr::ggarrange(busco_transcripts.gg,busco_proteins.gg, ncol=2, nrow=1, align = "v",
                  labels=c("a", "b"),common.legend = T, legend = "bottom")

# Save ggplot object
ggsave(busco_grouped.gg, file="plots/busco_complete_plot.png", dpi = 300, width = 3000, height = 1500, units = "px")

# Cut contig lengths into bins and beautify for the plot
trinity_cdhit90_ld.plot <- trinity_cdhit90_length %>%
  group_by(group = cut(length, breaks = seq(0, 3000, 200))) %>%
  summarise(n = sum(length)) %>%
  mutate(group = stringr::str_remove(group, "]")) %>%
  mutate(group = stringr::str_remove(group, "[(]")) %>%
  mutate(group = stringr::str_replace(group, ",", "-")) %>%
  mutate(group = c(group[-n()],  ">3e+03"))

# Override the ggplot x-axis alphabethical order
trinity_cdhit90_ld.plot$group <- factor(trinity_cdhit90_ld.plot$group, levels = trinity_cdhit90_ld.plot$group)

# Generate Figure 2 plot
trinity_cdhit90_ld.gg <-  trinity_cdhit90_ld.plot %>%
  ggplot(.) +
  geom_bar(aes(x=group, y=n), stat = "identity", fill=penguin_corp_color("pistachio")) +
  theme_classic() +
  ylab("No. of unigenes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank(),
        legend.position = "none")

# Save the plot
ggsave(trinity_cdhit90_ld.gg, file="plots/length_distribution_cdhit90.png", dpi = 300, width = 3000, height = 1500, units = "px")

# Split BLASTp results
blastp <- split_blastp(trinotate_data)
# Import qcovs not present in trinotate table
blast_add <- read.delim("data/Trinity_cdhit90_blastp_additional.txt", row.names = NULL) %>%
  mutate(transcript = stringr::str_split_i(transcript, "::", 2)) %>%
  select(transcript, qcovs)
# Join both tables 
blastp <- left_join(blastp, blast_add, by = "transcript")
# Summarize BLASTp results
blastp_sum <- summary_blast(blastp)
# Print them out
blastp_sum %>%
  group_by(genus) %>%
  summarise(transcripts = sum(transcripts)) %>%
  arrange(desc(transcripts))

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

# Write blastp_filter to file for filtration of the original fasta
write.table(blastp_filter, file="tables/blastp_filter.txt", sep="\t", col.names = F, row.names = F, quote = FALSE)

# Filter the blastp results
blastp_filtered <- blastp %>%
  filter(., !transcript %in% blastp_filter$transcript)

# Summarize filetered BLASTp results
blastp_sum_filtered <- summary_blast(blastp_filtered)

# Print them out
blastp_sum_filtered %>%
  group_by(genus) %>%
  summarise(transcripts = sum(transcripts)) %>%
  arrange(desc(transcripts))

# Split the lineage and extract higher taxonomic groups
blastp_filtered <- blastp_filtered %>%
  filter(., !transcript %in% blastp_filter$transcript) %>%
  mutate(tax = case_when(
    grepl("Viridiplantae", lineage) == TRUE ~ "Viridiplantae",
    grepl("Bacteria", lineage) == TRUE ~ "Bacteria",
    grepl("Metazoa", lineage) == TRUE ~ "Metazoa",
    grepl("Fungi", lineage) == TRUE ~ "Fungi",
    .default = "Other"
  ))

# Generate Figure 3a plot
blastp_tax.gg <- summary_blast_exp(blastp_filtered) %>% 
  group_by(tax) %>%
  summarise(transcripts = sum(transcripts)) %>%
  dplyr::slice(1:5) %>%
  arrange(desc(transcripts)) %>%
  ggplot(., aes(x="", y=transcripts, fill=tax)) +
   geom_bar(stat = "identity", width=1) +
   coord_polar("y", start=0) +
   theme_void() +
   ylab("No. of transcripts") +
   labs(fill="Higher taxonomies") +
   scale_fill_penguin(palette = "fiver") +
   theme(axis.text = element_blank())
# Save the plot
ggsave(blastp_tax.gg, file="plots/blastp_tax.png", dpi=300, width=1500, height = 1500, units = "px")

# Prepare the results for Figure 2b
blastp_genus.plot <- blastp_sum_filtered %>%
  group_by(genus) %>%
  summarise(transcripts = sum(transcripts)) %>%
  arrange(desc(transcripts)) %>%
  dplyr::top_n(10, transcripts)

# Override the ggplot x-axis alphabethical order
blastp_genus.plot$genus <- factor(blastp_genus.plot$genus, levels = blastp_genus.plot$genus)

# Generate Figure 3b plot
blastp_genus.gg <- blastp_genus.plot %>%
  ggplot(., aes(x="", y=transcripts, fill=genus)) +
    geom_bar(stat = "identity", width=1) +
    coord_polar("y", start=0, direction = -1) +
    theme_void() +
    ylab("No. of transcripts") +
    labs(fill="Genera") +
    scale_fill_penguin(palette = "middle-main") +
    theme(axis.text = element_blank())

# Save the plot
ggsave(blastp_genus.gg, file="plots/blastp_genus.png", dpi=300, width=1500, height = 1500, units = "px")

# Save a combined plot
blastp_fig <- ggpubr::ggarrange(blastp_tax.gg, blastp_genus.gg, ncol=2, nrow = 1, align = "v", labels=c("a", "b"))
ggsave(blastp_fig, file="plots/blast_fig.png", dpi = 300, width = 3000, height = 1500, units = "px")

# Filter the original trinotate table
trinotate_data_filtered <- trinotate_data %>%
  filter(., !transcript_id %in% blastp_filter$transcript) %>%
  mutate(abundance = rowSums(.[,19:23])) 

# Save filtered trinotate table
write.csv(trinotate_data_filtered, file="tables/trinotate_table_filtered.csv")

# Print original and filtered summary side-by-side for comparison
summart_table <- cbind(summary_trinotate(trinotate_data[,1:17]), summary_trinotate(trinotate_data_filtered[,1:17]))
write.table(summart_table, file="intermediate_objects/trinotate_results.txt", sep="\t", quote = FALSE)

# Split GO results from filtered trinotate table
go <- split_GOp(trinotate_data_filtered)

# Count and remove terms with less than 1000 counts
go_1000 <- go %>%
  select(transcript, ontology, name) %>%
  distinct() %>%
  dplyr::count(ontology, name, sort=TRUE) %>%
  filter(n > 1000) %>%
  arrange(ontology)

# Override the ggplot x-axis alphabethical order
go_1000$name <- factor(go_1000$name, levels = go_1000$name)

# Generate Figure 4 plot
go_plot.gg <- go_1000 %>%
  ggplot(.) +
  geom_bar(aes(x=name, y=n, fill=ontology), stat = "identity") +
  theme_classic() +
  ylab("No. of transcripts") +
  scale_fill_penguin(palette = "fiver") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank())

# Save the ggplot object
ggsave(go_plot.gg, file="plots/GO_terms_plot.png", dpi = 300, width = 3000, height = 1500, units = "px")

# Count the number of ontology terms
go %>%
  select(transcript, ontology) %>%
  distinct() %>%
  dplyr::count(ontology) %>%
  mutate(percent = round(n*100/sum(n),1))

# Calculate the top biological process terms
go %>%
  filter(ontology == "biological_process") %>%
  select(transcript, name) %>%
  distinct() %>%
  dplyr::count(name) %>%
  arrange(desc(n)) %>%
  mutate(percent = round(n*100/sum(n),1))

# Calculate the top cell components terms
go %>%
  filter(ontology == "cellular_component") %>%
  select(transcript, name) %>%
  distinct() %>%
  dplyr::count(name) %>%
  arrange(desc(n)) %>%
  mutate(percent = round(n*100/sum(n),1))

# Calculate the top molecular functions terms
go %>%
  filter(ontology == "molecular_function") %>%
  select(transcript, name) %>%
  distinct() %>%
  dplyr::count(name) %>%
  arrange(desc(n)) %>%
  mutate(percent = round(n*100/sum(n),1))

# Prepare KEGG terms for annotation
ath <- trinotate_data_filtered %>%
  filter(grepl(":ath:", Kegg)) %>%
  mutate(Kegg = stringr::str_remove(Kegg, "KEGG:")) %>%
  mutate(Kegg = stringr::str_split_i(Kegg, ":", 2)) 

# Prepare Arabidopsis thaliana KEGG terms
kegg_genes <- getGeneKEGGLinks(species.KEGG = "ath", convert = FALSE)
kegg_links <- getKEGGPathwayNames(species.KEGG = "ath", remove.qualifier = FALSE)

# Identify and prepare KEGG terms for plotting
kegg <- ath %>%
  select(Kegg) %>%
  dplyr::rename(GeneID = Kegg) %>%
  left_join(., kegg_genes, by="GeneID", multiple="all") %>%
  mutate(PathwayID = stringr::str_remove(PathwayID, "path:")) %>%
  left_join(., kegg_links, by="PathwayID") %>%
  mutate(Description = stringr::str_split_i(Description, " - ",1)) %>%
  group_by(Description) %>%
  summarise(Count = n()) %>%
  arrange(., desc(Count))

# Filter top 10 kegg terms
kegg_top10 <- kegg %>%
  top_n(10)

# Override the ggplot x-axis alphabethical order
kegg_top10$Description <- factor(kegg_top10$Description, levels = kegg_top10$Description)

# Generate Figure 4 plot
kegg_top10.gg <- kegg_top10 %>%
  filter(!is.na(Description)) %>%
  ggplot(.,aes(x="", y=Count, fill=Description)) +
  geom_bar(stat = "identity", width=1) +
  coord_polar("y", start=0, direction = 1) +
  theme_void() +
  ylab("No. of transcripts") +
  labs(fill="Top 10 KEGG terms") +
  scale_fill_penguin(palette = "middle-main") +
  theme(axis.text = element_blank())

# Save ggplot object
ggsave(kegg_top10.gg, file="plots/kegg_terms_plot.png", dpi = 300, width = 3000, height = 1500, units = "px")  
