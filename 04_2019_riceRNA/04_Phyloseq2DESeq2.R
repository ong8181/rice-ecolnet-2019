####
#### CER rice 2019
#### Rice RNA expression: phyloseq to DESeq2
####

# Load workspace
load("03_RNAdata2PhyloseqOut/03_RNAdata2PhyloseqOut.RData")

# create directory
output_folder04 <- "04_Phyloseq2DESeq2Out"
dir.create(output_folder04)

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)

# Load library
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2019.12.9
library(lubridate); packageVersion("lubridate") # 1.7.9.2, 2021.1.18
library(phyloseq); packageVersion("phyloseq") # 1.34.0, 2021.1.18
library(cowplot); packageVersion("cowplot") # 1.1.0, 2021.1.18
library(ggsci); packageVersion("ggsci") # 2.9, 2019.12.9
theme_set(theme_cowplot())
## DEG packages
library(DESeq2); packageVersion("DESeq2") # 1.30.0, 2021.1.18
library(ggpubr); packageVersion("ggpubr") # 0.4.0, 2021.1.18

# Generate subset phyloseq objects
ps_rna_top

# Remove genes with a small number of sequence reads
quantile(taxa_sums(ps_rna_top), probs = seq(0, 1, 0.1))
all(sample_sums(ps_rna_top) > 0)

# Divide samples into each sampling event (S1-S4)
ps_rna_s1 <- subset_samples(ps_rna_top, sampling_event == "S1") %>% prune_taxa(taxa_sums(.) > 0, .)
ps_rna_s2 <- subset_samples(ps_rna_top, sampling_event == "S2") %>% prune_taxa(taxa_sums(.) > 0, .)
ps_rna_s3 <- subset_samples(ps_rna_top, sampling_event == "S3") %>% prune_taxa(taxa_sums(.) > 0, .)
ps_rna_s4 <- subset_samples(ps_rna_top, sampling_event == "S4") %>% prune_taxa(taxa_sums(.) > 0, .)

# Phyloseq to DESeq2
deseq_rna_s1 <- phyloseq_to_deseq2(ps_rna_s1, ~ treatment)
deseq_rna_s2 <- phyloseq_to_deseq2(ps_rna_s2, ~ treatment)
deseq_rna_s3 <- phyloseq_to_deseq2(ps_rna_s3, ~ treatment)
deseq_rna_s4 <- phyloseq_to_deseq2(ps_rna_s4, ~ treatment)

# Do DESeq normalization
## Sampling Event 1 (Before treatment: 2019/6/23)
deseq_rna_s1 <- DESeq(deseq_rna_s1, test = "Wald", fitType = "parametric")
deseq_res_s1_1 <- results(deseq_rna_s1, contrast = c("treatment", "PN", "CT"), cooksCutoff = FALSE)
deseq_res_s1_2 <- results(deseq_rna_s1, contrast = c("treatment", "RM", "CT"), cooksCutoff = FALSE)
deseq_res_s1_3 <- results(deseq_rna_s1, contrast = c("treatment", "PN", "RM"), cooksCutoff = FALSE)

## Sampling Event 2 (After treatment: 2019/6/29)
deseq_rna_s2 <- DESeq(deseq_rna_s2, test = "Wald", fitType = "parametric")
deseq_res_s2_1 <- results(deseq_rna_s2, contrast = c("treatment", "PN", "CT"), cooksCutoff = FALSE)
deseq_res_s2_2 <- results(deseq_rna_s2, contrast = c("treatment", "RM", "CT"), cooksCutoff = FALSE)
deseq_res_s2_3 <- results(deseq_rna_s2, contrast = c("treatment", "PN", "RM"), cooksCutoff = FALSE)

## Sampling Event 3 (After treatment: 2019/7/12)
deseq_rna_s3 <- DESeq(deseq_rna_s3, test = "Wald", fitType = "parametric")
deseq_res_s3_1 <- results(deseq_rna_s3, contrast = c("treatment", "PN", "CT"), cooksCutoff = FALSE)
deseq_res_s3_2 <- results(deseq_rna_s3, contrast = c("treatment", "RM", "CT"), cooksCutoff = FALSE)
deseq_res_s3_3 <- results(deseq_rna_s3, contrast = c("treatment", "PN", "RM"), cooksCutoff = FALSE)

## Sampling Event 4 (After treatment, before heading: 2019/8/5)
deseq_rna_s4 <- DESeq(deseq_rna_s4, test = "Wald", fitType = "parametric")
deseq_res_s4_1 <- results(deseq_rna_s4, contrast = c("treatment", "PN", "CT"), cooksCutoff = FALSE)
deseq_res_s4_2 <- results(deseq_rna_s4, contrast = c("treatment", "RM", "CT"), cooksCutoff = FALSE)
deseq_res_s4_3 <- results(deseq_rna_s4, contrast = c("treatment", "PN", "RM"), cooksCutoff = FALSE)


# Visualize results
g1_1 <- ggmaplot(deseq_res_s1_1, main = expression("Sampling event 1, All locations: PN v.s. CT"),
                 fdr = 0.05, size = 0.4, fc = 1.5, 
                 palette = c("#B31B21", "#1465AC", "darkgray"),
                 genenames = as.vector(rownames(deseq_res_s1_1)),
                 #legend = "top", top = 20,
                 font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                 ggtheme = theme_bw()) + ylim(-6, 6) + xlim(-1, 18)
g1_2 <- ggmaplot(deseq_res_s1_2, main = expression("Sampling event 1, All locations: RM v.s. CT"),
                 fdr = 0.05, size = 0.4, fc = 1.5,
                 palette = c("#B31B21", "#1465AC", "darkgray"),
                 genenames = as.vector(rownames(deseq_res_s1_2)),
                 #legend = "top", top = 20,
                 font.label = c("bold", 11), font.legend = "bold", font.main = "bold",
                 ggtheme = theme_bw()) + ylim(-6, 6) + xlim(-1, 18)
g2_1 <- ggmaplot(deseq_res_s2_1, main = expression("Sampling event 2, All locations: PN v.s. CT"),
                 fdr = 0.05, size = 0.6, fc = 1.5, 
                 palette = c("#B31B21", "#1465AC", "darkgray"),
                 genenames = as.vector(rownames(deseq_res_s2_1)),
                 #legend = "top", top = 20,
                 #font.label = c("bold", 11),
                 top = 0,
                 #font.legend = "bold",font.main = "bold",
                 ggtheme = theme_bw()) + ylim(-6, 6) + xlim(-1, 18)
g2_2 <- ggmaplot(deseq_res_s2_2, main = expression("Sampling event 2, All locations: RM v.s. CT"),
                 fdr = 0.05, size = 0.6, fc = 1.5,
                 palette = c("#B31B21", "#1465AC", "darkgray"),
                 genenames = as.vector(rownames(deseq_res_s2_2)),
                 #legend = "top", top = 20,
                 top = 0,
                 #font.label = c("bold", 11), font.legend = "bold", font.main = "bold",
                 ggtheme = theme_bw()) + ylim(-6, 6) + xlim(-1, 18)
g3_1 <- ggmaplot(deseq_res_s3_1, main = expression("Sampling event 3, All locations: PN v.s. CT"),
                 fdr = 0.05, size = 0.4, fc = 1.5, 
                 palette = c("#B31B21", "#1465AC", "darkgray"),
                 genenames = as.vector(rownames(deseq_res_s3_1)),
                 #legend = "top", top = 20,
                 font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                 ggtheme = theme_bw()) + ylim(-6, 6) + xlim(-1, 18)
g3_2 <- ggmaplot(deseq_res_s3_2, main = expression("Sampling event 3, All locations: RM v.s. CT"),
                 fdr = 0.05, size = 0.4, fc = 1.5,
                 palette = c("#B31B21", "#1465AC", "darkgray"),
                 genenames = as.vector(rownames(deseq_res_s3_2)),
                 #legend = "top", top = 20,
                 font.label = c("bold", 11), font.legend = "bold", font.main = "bold",
                 ggtheme = theme_bw()) + ylim(-6, 6) + xlim(-1, 18)
g4_1 <- ggmaplot(deseq_res_s4_1, main = expression("Sampling event 4, All locations: PN v.s. CT"),
                 fdr = 0.05, size = 0.4, fc = 1.5, 
                 palette = c("#B31B21", "#1465AC", "darkgray"),
                 genenames = as.vector(rownames(deseq_res_s4_1)),
                 #legend = "top", top = 20,
                 font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                 ggtheme = theme_bw()) + ylim(-6, 6) + xlim(-1, 18)
g4_2 <- ggmaplot(deseq_res_s4_2, main = expression("Sampling event 4, All locations: RM v.s. CT"),
                 fdr = 0.05, size = 0.4, fc = 1.5,
                 palette = c("#B31B21", "#1465AC", "darkgray"),
                 genenames = as.vector(rownames(deseq_res_s4_2)),
                 #legend = "top", top = 20,
                 font.label = c("bold", 11), font.legend = "bold", font.main = "bold",
                 ggtheme = theme_bw()) + ylim(-6, 6) + xlim(-1, 18)

# Save DEG figure
ggsave(sprintf("%s/DEG_S1_AllLoc.pdf", output_folder04), plot = plot_grid(g1_1, g1_2, nrow = 1), width = 12, height =7)
ggsave(sprintf("%s/DEG_S2_AllLoc.pdf", output_folder04), plot = plot_grid(g2_1, g2_2, nrow = 1), width = 12, height =7)
ggsave(sprintf("%s/DEG_S3_AllLoc.pdf", output_folder04), plot = plot_grid(g3_1, g3_2, nrow = 1), width = 12, height =7)
ggsave(sprintf("%s/DEG_S4_AllLoc.pdf", output_folder04), plot = plot_grid(g4_1, g4_2, nrow = 1), width = 12, height =7)

# Save and output results
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder04, output_folder04))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder04, substr(Sys.time(), 1, 10)))
