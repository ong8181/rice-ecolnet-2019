####
#### CER rice 2019
#### Rice RNA expression: Position-specific DESeq2 analysis
####

# Load workspace
load("04_Phyloseq2DESeq2Out/04_Phyloseq2DESeq2Out.RData")

# create directory
output_folder05 <- "05_DESeq2LocationOut"
dir.create(output_folder05)

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

# Divide samples into each rice position (1-3)
ps_rna_s1l1 <- subset_samples(ps_rna_s1, location_id == "L1") %>% prune_taxa(taxa_sums(.) > 0, .)
ps_rna_s1l2 <- subset_samples(ps_rna_s1, location_id == "L2") %>% prune_taxa(taxa_sums(.) > 0, .)
ps_rna_s1l3 <- subset_samples(ps_rna_s1, location_id == "L3") %>% prune_taxa(taxa_sums(.) > 0, .)

ps_rna_s2l1 <- subset_samples(ps_rna_s2, location_id == "L1") %>% prune_taxa(taxa_sums(.) > 0, .)
ps_rna_s2l2 <- subset_samples(ps_rna_s2, location_id == "L2") %>% prune_taxa(taxa_sums(.) > 0, .)
ps_rna_s2l3 <- subset_samples(ps_rna_s2, location_id == "L3") %>% prune_taxa(taxa_sums(.) > 0, .)

ps_rna_s3l1 <- subset_samples(ps_rna_s3, location_id == "L1") %>% prune_taxa(taxa_sums(.) > 0, .)
ps_rna_s3l2 <- subset_samples(ps_rna_s3, location_id == "L2") %>% prune_taxa(taxa_sums(.) > 0, .)
ps_rna_s3l3 <- subset_samples(ps_rna_s3, location_id == "L3") %>% prune_taxa(taxa_sums(.) > 0, .)

ps_rna_s4l1 <- subset_samples(ps_rna_s4, location_id == "L1") %>% prune_taxa(taxa_sums(.) > 0, .)
ps_rna_s4l2 <- subset_samples(ps_rna_s4, location_id == "L2") %>% prune_taxa(taxa_sums(.) > 0, .)
ps_rna_s4l3 <- subset_samples(ps_rna_s4, location_id == "L3") %>% prune_taxa(taxa_sums(.) > 0, .)

# Phyloseq to DESeq2
deseq_rna_s1l1 <- phyloseq_to_deseq2(ps_rna_s1l1, ~ treatment)
deseq_rna_s1l2 <- phyloseq_to_deseq2(ps_rna_s1l2, ~ treatment)
deseq_rna_s1l3 <- phyloseq_to_deseq2(ps_rna_s1l3, ~ treatment)

deseq_rna_s2l1 <- phyloseq_to_deseq2(ps_rna_s2l1, ~ treatment)
deseq_rna_s2l2 <- phyloseq_to_deseq2(ps_rna_s2l2, ~ treatment)
deseq_rna_s2l3 <- phyloseq_to_deseq2(ps_rna_s2l3, ~ treatment)

deseq_rna_s3l1 <- phyloseq_to_deseq2(ps_rna_s3l1, ~ treatment)
deseq_rna_s3l2 <- phyloseq_to_deseq2(ps_rna_s3l2, ~ treatment)
deseq_rna_s3l3 <- phyloseq_to_deseq2(ps_rna_s3l3, ~ treatment)

deseq_rna_s4l1 <- phyloseq_to_deseq2(ps_rna_s4l1, ~ treatment)
deseq_rna_s4l2 <- phyloseq_to_deseq2(ps_rna_s4l2, ~ treatment)
deseq_rna_s4l3 <- phyloseq_to_deseq2(ps_rna_s4l3, ~ treatment)

# Do DESeq normalization
## Sampling Event 1 (Before treatment: 2019/6/23)
deseq_rna_s1l1 <- DESeq(deseq_rna_s1l1, test = "Wald", fitType = "parametric")
deseq_res_s1l1_1 <- results(deseq_rna_s1l1, contrast = c("treatment", "PN", "CT"), cooksCutoff = FALSE)
deseq_res_s1l1_2 <- results(deseq_rna_s1l1, contrast = c("treatment", "RM", "CT"), cooksCutoff = FALSE)
deseq_rna_s1l2 <- DESeq(deseq_rna_s1l2, test = "Wald", fitType = "parametric")
deseq_res_s1l2_1 <- results(deseq_rna_s1l2, contrast = c("treatment", "PN", "CT"), cooksCutoff = FALSE)
deseq_res_s1l2_2 <- results(deseq_rna_s1l2, contrast = c("treatment", "RM", "CT"), cooksCutoff = FALSE)
deseq_rna_s1l3 <- DESeq(deseq_rna_s1l3, test = "Wald", fitType = "parametric")
deseq_res_s1l3_1 <- results(deseq_rna_s1l3, contrast = c("treatment", "PN", "CT"), cooksCutoff = FALSE)
deseq_res_s1l3_2 <- results(deseq_rna_s1l3, contrast = c("treatment", "RM", "CT"), cooksCutoff = FALSE)

## Sampling Event 2 (After treatment: 2019/6/29)
deseq_rna_s2l1 <- DESeq(deseq_rna_s2l1, test = "Wald", fitType = "parametric")
deseq_res_s2l1_1 <- results(deseq_rna_s2l1, contrast = c("treatment", "PN", "CT"), cooksCutoff = FALSE)
deseq_res_s2l1_2 <- results(deseq_rna_s2l1, contrast = c("treatment", "RM", "CT"), cooksCutoff = FALSE)
deseq_rna_s2l2 <- DESeq(deseq_rna_s2l2, test = "Wald", fitType = "parametric")
deseq_res_s2l2_1 <- results(deseq_rna_s2l2, contrast = c("treatment", "PN", "CT"), cooksCutoff = FALSE)
deseq_res_s2l2_2 <- results(deseq_rna_s2l2, contrast = c("treatment", "RM", "CT"), cooksCutoff = FALSE)
deseq_rna_s2l3 <- DESeq(deseq_rna_s2l3, test = "Wald", fitType = "parametric")
deseq_res_s2l3_1 <- results(deseq_rna_s2l3, contrast = c("treatment", "PN", "CT"), cooksCutoff = FALSE)
deseq_res_s2l3_2 <- results(deseq_rna_s2l3, contrast = c("treatment", "RM", "CT"), cooksCutoff = FALSE)

## Sampling Event 3 (After treatment: 2019/7/12)
deseq_rna_s3l1 <- DESeq(deseq_rna_s3l1, test = "Wald", fitType = "parametric")
deseq_res_s3l1_1 <- results(deseq_rna_s3l1, contrast = c("treatment", "PN", "CT"), cooksCutoff = FALSE)
deseq_res_s3l1_2 <- results(deseq_rna_s3l1, contrast = c("treatment", "RM", "CT"), cooksCutoff = FALSE)
deseq_rna_s3l2 <- DESeq(deseq_rna_s3l2, test = "Wald", fitType = "parametric")
deseq_res_s3l2_1 <- results(deseq_rna_s3l2, contrast = c("treatment", "PN", "CT"), cooksCutoff = FALSE)
deseq_res_s3l2_2 <- results(deseq_rna_s3l2, contrast = c("treatment", "RM", "CT"), cooksCutoff = FALSE)
deseq_rna_s3l3 <- DESeq(deseq_rna_s3l3, test = "Wald", fitType = "parametric")
deseq_res_s3l3_1 <- results(deseq_rna_s3l3, contrast = c("treatment", "PN", "CT"), cooksCutoff = FALSE)
deseq_res_s3l3_2 <- results(deseq_rna_s3l3, contrast = c("treatment", "RM", "CT"), cooksCutoff = FALSE)

## Sampling Event 4 (After treatment, before heading: 2019/8/5)
deseq_rna_s4l1 <- DESeq(deseq_rna_s4l1, test = "Wald", fitType = "parametric")
deseq_res_s4l1_1 <- results(deseq_rna_s4l1, contrast = c("treatment", "PN", "CT"), cooksCutoff = FALSE)
deseq_res_s4l1_2 <- results(deseq_rna_s4l1, contrast = c("treatment", "RM", "CT"), cooksCutoff = FALSE)
deseq_rna_s4l2 <- DESeq(deseq_rna_s4l2, test = "Wald", fitType = "parametric")
deseq_res_s4l2_1 <- results(deseq_rna_s4l2, contrast = c("treatment", "PN", "CT"), cooksCutoff = FALSE)
deseq_res_s4l2_2 <- results(deseq_rna_s4l2, contrast = c("treatment", "RM", "CT"), cooksCutoff = FALSE)
deseq_rna_s4l3 <- DESeq(deseq_rna_s4l3, test = "Wald", fitType = "parametric")
deseq_res_s4l3_1 <- results(deseq_rna_s4l3, contrast = c("treatment", "PN", "CT"), cooksCutoff = FALSE)
deseq_res_s4l3_2 <- results(deseq_rna_s4l3, contrast = c("treatment", "RM", "CT"), cooksCutoff = FALSE)


# Visualize results
## Sampling event 1, PN v.s. CT, Location 1-3
s1_pn1 <- ggmaplot(deseq_res_s1l1_1, main = expression("Sampling event 1, Loc.1: PN v.s. CT"),
                   fdr = 0.05, size = 0.4, fc = 1.5, 
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(rownames(deseq_res_s1l1_1)),
                   font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                   ggtheme = theme_bw()) + ylim(-8, 8) + xlim(-0.5, 19)
s1_pn2 <- ggmaplot(deseq_res_s1l2_1, main = expression("Sampling event 1, Loc.2: PN v.s. CT"),
                   fdr = 0.05, size = 0.4, fc = 1.5, 
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(rownames(deseq_res_s1l2_1)),
                   font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                   ggtheme = theme_bw()) + ylim(-8, 8) + xlim(-0.5, 19)
s1_pn3 <- ggmaplot(deseq_res_s1l3_1, main = expression("Sampling event 1, Loc.3: PN v.s. CT"),
                   fdr = 0.05, size = 0.4, fc = 1.5, 
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(rownames(deseq_res_s1l3_1)),
                   font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                   ggtheme = theme_bw()) + ylim(-8, 8) + xlim(-0.5, 19)
## Sampling event 1, RM v.s. CT, Location 1-3
s1_rm1 <- ggmaplot(deseq_res_s1l1_2, main = expression("Sampling event 1, Loc.1: RM v.s. CT"),
                   fdr = 0.05, size = 0.4, fc = 1.5, 
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(rownames(deseq_res_s1l1_2)),
                   font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                   ggtheme = theme_bw()) + ylim(-8, 8) + xlim(-0.5, 19)
s1_rm2 <- ggmaplot(deseq_res_s1l2_2, main = expression("Sampling event 1, Loc.2: RM v.s. CT"),
                   fdr = 0.05, size = 0.4, fc = 1.5, 
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(rownames(deseq_res_s1l2_2)),
                   font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                   ggtheme = theme_bw()) + ylim(-8, 8) + xlim(-0.5, 19)
s1_rm3 <- ggmaplot(deseq_res_s1l3_2, main = expression("Sampling event 1, Loc.3: RM v.s. CT"),
                   fdr = 0.05, size = 0.4, fc = 1.5, 
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(rownames(deseq_res_s1l3_2)),
                   font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                   ggtheme = theme_bw()) + ylim(-8, 8) + xlim(-0.5, 19)


## Sampling event 2, PN v.s. CT, Location 1-3
s2_pn1 <- ggmaplot(deseq_res_s2l1_1, main = expression("Sampling event 2, Loc.1: PN v.s. CT"),
                   fdr = 0.05, size = 0.4, fc = 1.5, 
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(rownames(deseq_res_s2l1_1)),
                   font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                   ggtheme = theme_bw()) + ylim(-8, 8) + xlim(-0.5, 19)
s2_pn2 <- ggmaplot(deseq_res_s2l2_1, main = expression("Sampling event 2, Loc.2: PN v.s. CT"),
                   fdr = 0.05, size = 0.4, fc = 1.5, 
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(rownames(deseq_res_s2l2_1)),
                   font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                   ggtheme = theme_bw()) + ylim(-8, 8) + xlim(-0.5, 19)
s2_pn3 <- ggmaplot(deseq_res_s2l3_1, main = expression("Sampling event 2, Loc.3: PN v.s. CT"),
                   fdr = 0.05, size = 0.4, fc = 1.5, 
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(rownames(deseq_res_s2l3_1)),
                   font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                   ggtheme = theme_bw()) + ylim(-8, 8) + xlim(-0.5, 19)
## Sampling event 2, RM v.s. CT, Location 1-3
s2_rm1 <- ggmaplot(deseq_res_s2l1_2, main = expression("Sampling event 2, Loc.1: RM v.s. CT"),
                   fdr = 0.05, size = 0.4, fc = 1.5, 
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(rownames(deseq_res_s2l1_2)),
                   font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                   ggtheme = theme_bw()) + ylim(-8, 8) + xlim(-0.5, 19)
s2_rm2 <- ggmaplot(deseq_res_s2l2_2, main = expression("Sampling event 2, Loc.2: RM v.s. CT"),
                   fdr = 0.05, size = 0.4, fc = 1.5, 
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(rownames(deseq_res_s2l2_2)),
                   font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                   ggtheme = theme_bw()) + ylim(-8, 8) + xlim(-0.5, 19)
s2_rm3 <- ggmaplot(deseq_res_s2l3_2, main = expression("Sampling event 2, Loc.3: RM v.s. CT"),
                   fdr = 0.05, size = 0.4, fc = 1.5, 
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(rownames(deseq_res_s2l3_2)),
                   font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                   ggtheme = theme_bw()) + ylim(-8, 8) + xlim(-0.5, 19)


## Sampling event 3, PN v.s. CT, Location 1-3
s3_pn1 <- ggmaplot(deseq_res_s3l1_1, main = expression("Sampling event 3, Loc.1: PN v.s. CT"),
                   fdr = 0.05, size = 0.4, fc = 1.5, 
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(rownames(deseq_res_s3l1_1)),
                   font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                   ggtheme = theme_bw()) + ylim(-8, 8) + xlim(-0.5, 19)
s3_pn2 <- ggmaplot(deseq_res_s3l2_1, main = expression("Sampling event 3, Loc.2: PN v.s. CT"),
                   fdr = 0.05, size = 0.4, fc = 1.5, 
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(rownames(deseq_res_s3l2_1)),
                   font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                   ggtheme = theme_bw()) + ylim(-8, 8) + xlim(-0.5, 19)
s3_pn3 <- ggmaplot(deseq_res_s3l3_1, main = expression("Sampling event 3, Loc.3: PN v.s. CT"),
                   fdr = 0.05, size = 0.4, fc = 1.5, 
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(rownames(deseq_res_s3l3_1)),
                   font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                   ggtheme = theme_bw()) + ylim(-8, 8) + xlim(-0.5, 19)
## Sampling event 3, RM v.s. CT, Location 1-3
s3_rm1 <- ggmaplot(deseq_res_s3l1_2, main = expression("Sampling event 3, Loc.1: RM v.s. CT"),
                   fdr = 0.05, size = 0.4, fc = 1.5, 
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(rownames(deseq_res_s3l1_2)),
                   font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                   ggtheme = theme_bw()) + ylim(-8, 8) + xlim(-0.5, 19)
s3_rm2 <- ggmaplot(deseq_res_s3l2_2, main = expression("Sampling event 3, Loc.2: RM v.s. CT"),
                   fdr = 0.05, size = 0.4, fc = 1.5, 
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(rownames(deseq_res_s3l2_2)),
                   font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                   ggtheme = theme_bw()) + ylim(-8, 8) + xlim(-0.5, 19)
s3_rm3 <- ggmaplot(deseq_res_s3l3_2, main = expression("Sampling event 3, Loc.3: RM v.s. CT"),
                   fdr = 0.05, size = 0.4, fc = 1.5, 
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(rownames(deseq_res_s3l3_2)),
                   font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                   ggtheme = theme_bw()) + ylim(-8, 8) + xlim(-0.5, 19)

## Sampling event 4, PN v.s. CT, Location 1-3
s4_pn1 <- ggmaplot(deseq_res_s4l1_1, main = expression("Sampling event 4, Loc.1: PN v.s. CT"),
                   fdr = 0.05, size = 0.4, fc = 1.5, 
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(rownames(deseq_res_s4l1_1)),
                   font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                   ggtheme = theme_bw()) + ylim(-8, 8) + xlim(-0.5, 19)
s4_pn2 <- ggmaplot(deseq_res_s4l2_1, main = expression("Sampling event 4, Loc.2: PN v.s. CT"),
                   fdr = 0.05, size = 0.4, fc = 1.5, 
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(rownames(deseq_res_s4l2_1)),
                   font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                   ggtheme = theme_bw()) + ylim(-8, 8) + xlim(-0.5, 19)
s4_pn3 <- ggmaplot(deseq_res_s4l3_1, main = expression("Sampling event 4, Loc.3: PN v.s. CT"),
                   fdr = 0.05, size = 0.4, fc = 1.5, 
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(rownames(deseq_res_s4l3_1)),
                   font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                   ggtheme = theme_bw()) + ylim(-8, 8) + xlim(-0.5, 19)
## Sampling event 4, RM v.s. CT, Location 1-3
s4_rm1 <- ggmaplot(deseq_res_s4l1_2, main = expression("Sampling event 4, Loc.1: RM v.s. CT"),
                   fdr = 0.05, size = 0.4, fc = 1.5, 
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(rownames(deseq_res_s4l1_2)),
                   font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                   ggtheme = theme_bw()) + ylim(-8, 8) + xlim(-0.5, 19)
s4_rm2 <- ggmaplot(deseq_res_s4l2_2, main = expression("Sampling event 4, Loc.2: RM v.s. CT"),
                   fdr = 0.05, size = 0.4, fc = 1.5, 
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(rownames(deseq_res_s4l2_2)),
                   font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                   ggtheme = theme_bw()) + ylim(-8, 8) + xlim(-0.5, 19)
s4_rm3 <- ggmaplot(deseq_res_s4l3_2, main = expression("Sampling event 4, Loc.3: RM v.s. CT"),
                   fdr = 0.05, size = 0.4, fc = 1.5, 
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(rownames(deseq_res_s4l3_2)),
                   font.label = c("bold", 11),　font.legend = "bold",font.main = "bold",
                   ggtheme = theme_bw()) + ylim(-8, 8) + xlim(-0.5, 19)



# Summarize plot
s1_all <- plot_grid(plot_grid(g1_1 + xlim(-0.5, 19), g1_2 + xlim(-0.5, 19), ncol = 1, align = "hv"),
                    plot_grid(s1_pn1, s1_rm1, ncol = 1, align = "hv"),
                    plot_grid(s1_pn2, s1_rm2, ncol = 1, align = "hv"),
                    plot_grid(s1_pn3, s1_rm3, ncol = 1, align = "hv"),
                    ncol = 4, align = "hvlr", rel_widths = c(1,.6,.6,.6))
s2_all <- plot_grid(plot_grid(g2_1 + xlim(-0.5, 19), g2_2 + xlim(-0.5, 19), ncol = 1, align = "hv"),
                    plot_grid(s2_pn1, s2_rm1, ncol = 1, align = "hv"),
                    plot_grid(s2_pn2, s2_rm2, ncol = 1, align = "hv"),
                    plot_grid(s2_pn3, s2_rm3, ncol = 1, align = "hv"),
                    ncol = 4, align = "hvlr", rel_widths = c(1,.6,.6,.6))
s3_all <- plot_grid(plot_grid(g3_1 + xlim(-0.5, 19), g3_2 + xlim(-0.5, 19), ncol = 1, align = "hv"),
                    plot_grid(s3_pn1, s3_rm1, ncol = 1, align = "hv"),
                    plot_grid(s3_pn2, s3_rm2, ncol = 1, align = "hv"),
                    plot_grid(s3_pn3, s3_rm3, ncol = 1, align = "hv"),
                    ncol = 4, align = "hvlr", rel_widths = c(1,.6,.6,.6))
s4_all <- plot_grid(plot_grid(g4_1 + xlim(-0.5, 19), g4_2 + xlim(-0.5, 19), ncol = 1, align = "hv"),
                    plot_grid(s4_pn1, s4_rm1, ncol = 1, align = "hv"),
                    plot_grid(s4_pn2, s4_rm2, ncol = 1, align = "hv"),
                    plot_grid(s4_pn3, s4_rm3, ncol = 1, align = "hv"),
                    ncol = 4, align = "hvlr", rel_widths = c(1,.6,.6,.6))

  
# Save DEG figure
ggsave(sprintf("%s/DEG_S1_All.jpg", output_folder05), plot = s1_all, width = 22, height = 12)
ggsave(sprintf("%s/DEG_S2_All.jpg", output_folder05), plot = s2_all, width = 22, height = 12)
ggsave(sprintf("%s/DEG_S3_All.jpg", output_folder05), plot = s3_all, width = 22, height = 12)
ggsave(sprintf("%s/DEG_S4_All.jpg", output_folder05), plot = s4_all, width = 22, height = 12)

ggsave(sprintf("%s/DEG_S1_All.pdf", output_folder05), plot = s1_all, width = 22, height = 12)
ggsave(sprintf("%s/DEG_S2_All.pdf", output_folder05), plot = s2_all, width = 22, height = 12)
ggsave(sprintf("%s/DEG_S3_All.pdf", output_folder05), plot = s3_all, width = 22, height = 12)
ggsave(sprintf("%s/DEG_S4_All.pdf", output_folder05), plot = s4_all, width = 22, height = 12)

# Delete temporal objects
rm(g1_1); rm(g1_2); rm(g2_1); rm(g2_2); rm(g3_1); rm(g3_2)
rm(ps_rna_s1l1);rm(ps_rna_s1l2);rm(ps_rna_s1l3)
rm(ps_rna_s2l1);rm(ps_rna_s2l2);rm(ps_rna_s2l3)
rm(ps_rna_s3l1);rm(ps_rna_s3l2);rm(ps_rna_s3l3)
rm(ps_rna_s4l1);rm(ps_rna_s4l2);rm(ps_rna_s4l3)
rm(ps_rna_s1);rm(ps_rna_s2);rm(ps_rna_s3);rm(ps_rna_s4)

# Save and output results
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder05, output_folder05))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder05, substr(Sys.time(), 1, 10)))
