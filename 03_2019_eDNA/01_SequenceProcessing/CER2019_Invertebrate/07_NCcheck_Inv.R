####
#### CER rice 2019 eDNA analysis
#### No.7 CMR-011 Invertebrate: Negative & Positive control check
####

# Load workspace and data
load("06_ReadNSummary_InvOut/06_ReadNSummary_InvOut.RData")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder07 <- "07_NCcheck_InvOut"
dir.create(output_folder07)

# Load library and functions
library(phyloseq); packageVersion("phyloseq") # 1.32.0, 2021.1.14
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2021.1.13
library(lubridate); packageVersion("lubridate") # 1.7.9.2, 2021.1.14
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.1.13
library(ggsci); packageVersion("ggsci") # 2.9, 2021.1.13
theme_set(theme_cowplot())


# Preparetion to import to phyloseq
dim(sample_sheet); dim(tax_claident); dim(seqtab_nochim)
all(rownames(sample_sheet) == rownames(seqtab_nochim)) # sample name check

tax_claident2$seq <- colnames(seqtab_nochim) # save sequence info
tax_claident2$seqlen <- nchar(colnames(seqtab_nochim)) # calculate sequence length
colnames(seqtab_nochim) <- rownames(tax_claident2) # change col name
all(rownames(tax_claident2) == colnames(seqtab_nochim)) # taxa name check

tax_claident2$std_or_field <- "Field DNA"
tax_claident2[substr(tax_claident2$species, 1, nchar(std_seq_head)) == std_seq_head,]$std_or_field <- "Standard DNA"

# Import data to phyloseq
ps0 <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows=FALSE),
                sample_data(sample_sheet),
                tax_table(as.matrix(tax_claident2)))

# Generating table for manual NC check
sample_sheet_na <- matrix(rep(NaN, ncol(tax_claident2)*ncol(sample_sheet)), ncol=ncol(sample_sheet))
colnames(sample_sheet_na) <- colnames(sample_sheet)
sample_sheet_comb <- rbind(sample_sheet_na, as.matrix(sample_sheet))
tax_table_comb <- rbind(t(as.matrix(tax_claident2)), seqtab_nochim)
table_for_nc_check <- cbind(sample_sheet_comb, tax_table_comb)
table_for_nc_check[,"std_validity"][1:ncol(tax_claident2)] <- colnames(tax_claident2)
write.csv(table_for_nc_check, sprintf("%s/SummarizedTable_for_NCcheck.csv", output_folder07))

# Visualize pattern
sample_summary$sample_nc <- as.factor(sample_summary$sample_nc)
ps_m1 <- pivot_longer(sample_summary[,c("sample_nc", "STD_all", "NonSTD_all")],
                      cols = -sample_nc, names_to = "variable")
ps_m1$variable <- as.factor(ps_m1$variable)
p1 <- ggplot(ps_m1, aes(x = sample_nc, y = value, color = variable, group = variable:sample_nc)) +
  geom_boxplot(colour = "black", outlier.shape = NA) + #scale_y_log10() +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.3)) +
  scale_color_nejm(alpha = 0.5) + ylab("Sequence reads/sample") + xlab(NULL) +
  guides(color = guide_legend(title = "STD or non-STD")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_x_discrete(limits = c("sample", "std_nc", "field_nc", "pcr_nc"),
                   labels = c("field_nc" = "Field NC",
                              "pcr_nc" = "PCR NC",
                              "sample" = "Sample",
                              "std_nc" = "Standard NC"))
ggsave(sprintf("%s/SequenceReads_Barplot.pdf", output_folder07), plot = p1, width = 6, height = 6)


# Examination of ASVs that were detected in field NC, PCR NC and standard NC
# Extract taxa that are detected from NC samples
# Barplot of NC-detected taxa (Optional)
ps_nctax <- prune_taxa(taxa_sums(prune_samples(sample_data(ps0)$sample_nc != "sample", ps0)) > 0, ps0)
p2 <- plot_bar(ps_nctax, x = "sample_nc", fill = "phylum") +
  geom_bar(stat = "identity", colour = NA) + scale_fill_igv()

# Boxplot for NC-detected taxa
ps_m2 <- psmelt(ps_nctax) %>%
  filter(sample_nc != "sample") %>%
  filter(std_or_field == "Field DNA") %>%
  filter(Abundance > 0) %>% tibble()
ps_m2$phylum <- as.factor(ps_m2$phylum)
ps_m2$sample_nc <- as.factor(ps_m2$sample_nc)

p3 <- ggplot(ps_m2, aes(x = sample_nc, y = Abundance, group = phylum:sample_nc, fill = phylum)) +
  geom_boxplot(width = 0.6, position = position_dodge(width=0.7)) + scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_igv() + xlab(NULL) + ylab("Log10(Sequence reads)") +
  scale_x_discrete(limits = c("field_nc", "std_nc", "pcr_nc"),
                   labels = c("field_nc" = "Field NC",
                              "pcr_nc" = "PCR NC",
                              "std_nc" = "Standard NC"))
ggsave(sprintf("%s/SequenceReads_NCtaxaBoxplot.pdf", output_folder07), plot = p3, width = 8, height = 6)


# Examination of positive controls and field negative controls
# (Time series plot)
# (Individual sample plot: Field negative controls)
ps_nctax2 <- prune_taxa(taxa_sums(prune_samples(sample_data(ps0)$sample_nc == "field_nc", ps0)) > 0, ps0)
ps_m3 <- psmelt(ps_nctax2) %>%
  filter(sample_nc == "field_nc") %>%
  filter(std_or_field == "Field DNA") %>% tibble()
ps_m3$date <- ymd(ps_m3$date)
ps_m3$phylum <- as.factor(ps_m3$phylum)

p4 <- ggplot(ps_m3, aes(x = date, y = Abundance, fill = phylum)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90)) + scale_fill_igv()

ps_m4 <- psmelt(ps_nctax2) %>%
  filter(sample_nc == "field_nc") %>%
  filter(std_or_field == "Standard DNA") %>% tibble()
ps_m4$date <- ymd(ps_m4$date)
ps_m4$phylum <- as.factor(ps_m4$phylum)

p5 <- ggplot(ps_m4, aes(x = date, y = Abundance, fill = species)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90)) + scale_fill_igv()

field_nc_reads <- plot_grid(p4, p5, ncol = 2, labels = "auto", align = "hv")
ggsave(sprintf("%s/SequenceReads_FieldNCreads.pdf", output_folder07), plot = field_nc_reads, width = 12, height = 8)


# Examination of PCR negative controls and standard negative controls
# (Individual sample plot: PCR negative controls)
ps_nctax3 <- prune_taxa(taxa_sums(prune_samples(sample_data(ps0)$sample_nc == "std_nc" | sample_data(ps0)$sample_nc == "pcr_nc", ps0)) > 0, ps0)
ps_m5 <- psmelt(ps_nctax3) %>%
  filter(sample_nc == "pcr_nc") %>%
  filter(std_or_field == "Field DNA")
ps_m5$date <- ymd(ps_m5$date)
ps_m5$phylum <- as.factor(ps_m5$phylum)

p6 <- ggplot(ps_m5, aes(x = Sample, y = Abundance, fill = phylum)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_fill_igv() +
  ggtitle("PCR negative controls: Field DNA reads")

ps_m6 <- psmelt(ps_nctax3) %>%
  filter(sample_nc == "pcr_nc") %>%
  filter(std_or_field == "Standard DNA")
ps_m6$date <- ymd(ps_m6$date)
ps_m6$phylum <- as.factor(ps_m6$phylum)

p7 <- ggplot(ps_m6, aes(x = Sample, y = Abundance, fill = species)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_fill_igv() +
  ggtitle("PCR negative controls: Standard DNA reads")

pcr_nc_reads <- plot_grid(p6, p7, ncol = 1, labels = "auto", align = "hv")
ggsave(sprintf("%s/SequenceReads_PCRNCreads.pdf", output_folder07), plot = pcr_nc_reads, width = 8, height = 8)


# (Individual sample plot: Standard negative controls)
ps_m7 <- psmelt(ps_nctax3) %>%
  filter(sample_nc == "std_nc") %>%
  filter(std_or_field == "Field DNA")
ps_m7$date <- ymd(ps_m7$date)
ps_m7$phylum <- as.factor(ps_m7$phylum)

p8 <- ggplot(ps_m7, aes(x = Sample, y = Abundance, fill = phylum)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_fill_igv() +
  ggtitle("Standard negative controls: Field DNA reads")

ps_m8 <- psmelt(ps_nctax3) %>%
  filter(sample_nc == "std_nc") %>%
  filter(std_or_field == "Standard DNA")
ps_m8$date <- ymd(ps_m8$date)
ps_m8$phylum <- as.factor(ps_m8$phylum)

p9 <- ggplot(ps_m8, aes(x = Sample, y = Abundance, fill = species)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_fill_igv() +
  ggtitle("Standard negative controls: Standard DNA reads")

std_nc_reads <- plot_grid(p8, p9, ncol = 1, labels = "auto", align = "hv")
ggsave(sprintf("%s/SequenceReads_STDNCreads.pdf", output_folder07), plot = std_nc_reads, width = 10, height = 8)


# Calculating ratio
dna_summary <- aggregate(sample_summary$NonSTD_all, list(sample_summary$sample_nc), mean)
dna_summary$prop <- dna_summary$x/dna_summary[dna_summary$Group.1 == "sample","x"]
write.csv(dna_summary, sprintf("%s/DNAreads_ContaminationLevel.csv", output_folder07))

# Save and output results
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder07, output_folder07))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder07, substr(Sys.time(), 1, 10)))

