####
#### CER rice 2019
#### Figure code for rice RNA expression: phyloseq to DESeq2
####

# Load workspace
result_folder01 <- "../04_2019_riceRNA"
load(sprintf("%s/04_Phyloseq2DESeq2Out/04_Phyloseq2DESeq2Out.RData", result_folder01))

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

# Create output folder
fig_output <- "00_FigRaw"
#dir.create(fig_output)

# Revise the treatment names (2023.2.7)
treatment_ori <- c("CT", "PN", "RM")
treatment_rev <- c("CT", "GN", "CK")
sample_data(ps_rna_top) <- sample_data(ps_rna_top) %>% data.frame %>% 
  mutate(treatment = factor(treatment, levels = treatment_ori, labels = treatment_rev)) %>% sample_data


# <----------------------------------------------------> #
#                   Reload figures
# <----------------------------------------------------> #
## General overview
g1 <- ggplot(sample_data(ps_rna_top), aes(x = total_read00, fill = sampling_event)) +
  geom_histogram(alpha = 0.8) + scale_x_log10() + scale_fill_startrek() +
  xlab("Total reads") + ylab("Count") + ggtitle("Reads count")
genes_df <- data.frame(gene_name = names(taxa_sums(ps_rna_top)),
                       gene_exp = taxa_sums(ps_rna_top))
g2 <- ggplot(genes_df, aes(x = gene_exp)) +
  geom_histogram(alpha = 0.8) + scale_y_log10() +
  xlab("Reads count") + ylab("Count") + ggtitle("Genes expressed")


# <----------------------------------------------------> #
#  Revise titles
# <----------------------------------------------------> #
g1_1 <- g1_1 + ggtitle(expression(paste("Sampling event 1, All locations: CT v.s. GN")))
g1_2 <- g1_2 + ggtitle(expression(paste("Sampling event 1, All locations: CT v.s. CK")))
g2_1 <- g2_1 + ggtitle(expression(paste("Sampling event 2, All locations: CT v.s. GN")))
g2_2 <- g2_2 + ggtitle(expression(paste("Sampling event 2, All locations: CT v.s. CK")))
g3_1 <- g3_1 + ggtitle(expression(paste("Sampling event 3, All locations: CT v.s. GN")))
g3_2 <- g3_2 + ggtitle(expression(paste("Sampling event 3, All locations: CT v.s. CK")))
g4_1 <- g4_1 + ggtitle(expression(paste("Sampling event 4, All locations: CT v.s. GN")))
g4_2 <- g4_2 + ggtitle(expression(paste("Sampling event 4, All locations: CT v.s. CK")))


# <----------------------------------------------------> #
#                     Save figures
# <----------------------------------------------------> #
rna_hist_figs <- list(g1, g2)
maplot_all <- list(g1_1, g1_2, g2_1, g2_2, g3_1, g3_2, g4_1, g4_2)
saveRDS(rna_hist_figs, sprintf("%s/Fig_RNAhistogram.obj", fig_output))
saveRDS(maplot_all, sprintf("%s/Fig_RNAmaplot.obj", fig_output))


# <----------------------------------------------------> #
#           Refresh and load new workspace
# <----------------------------------------------------> #
# Delete all objects but "fig_output" and "result_folder01"
rm(list = setdiff(ls(), c("fig_output", "result_folder01")))
# Load workspace
load(sprintf("%s/05_DESeq2LocationOut/05_DESeq2LocationOut.RData", result_folder01))

# Revise titles
## G. nunn treatment
s1_pn1 <- s1_pn1 + ggtitle(expression(paste("Sampling event 1, Loc.1: CT v.s. GN")))
s1_pn2 <- s1_pn2 + ggtitle(expression(paste("Sampling event 1, Loc.2: CT v.s. GN")))
s1_pn3 <- s1_pn3 + ggtitle(expression(paste("Sampling event 1, Loc.3: CT v.s. GN")))
s2_pn1 <- s2_pn1 + ggtitle(expression(paste("Sampling event 2, Loc.1: CT v.s. GN")))
s2_pn2 <- s2_pn2 + ggtitle(expression(paste("Sampling event 2, Loc.2: CT v.s. GN")))
s2_pn3 <- s2_pn3 + ggtitle(expression(paste("Sampling event 2, Loc.3: CT v.s. GN")))
s3_pn1 <- s3_pn1 + ggtitle(expression(paste("Sampling event 3, Loc.1: CT v.s. GN")))
s3_pn2 <- s3_pn2 + ggtitle(expression(paste("Sampling event 3, Loc.2: CT v.s. GN")))
s3_pn3 <- s3_pn3 + ggtitle(expression(paste("Sampling event 3, Loc.3: CT v.s. GN")))
s4_pn1 <- s4_pn1 + ggtitle(expression(paste("Sampling event 4, Loc.1: CT v.s. GN")))
s4_pn2 <- s4_pn2 + ggtitle(expression(paste("Sampling event 4, Loc.2: CT v.s. GN")))
s4_pn3 <- s4_pn3 + ggtitle(expression(paste("Sampling event 4, Loc.3: CT v.s. GN")))
## C. kiiensis treatment
s1_rm1 <- s1_rm1 + ggtitle(expression(paste("Sampling event 1, Loc.1: CT v.s. CK")))
s1_rm2 <- s1_rm2 + ggtitle(expression(paste("Sampling event 1, Loc.2: CT v.s. CK")))
s1_rm3 <- s1_rm3 + ggtitle(expression(paste("Sampling event 1, Loc.3: CT v.s. CK")))
s2_rm1 <- s2_rm1 + ggtitle(expression(paste("Sampling event 2, Loc.1: CT v.s. CK")))
s2_rm2 <- s2_rm2 + ggtitle(expression(paste("Sampling event 2, Loc.2: CT v.s. CK")))
s2_rm3 <- s2_rm3 + ggtitle(expression(paste("Sampling event 2, Loc.3: CT v.s. CK")))
s3_rm1 <- s3_rm1 + ggtitle(expression(paste("Sampling event 3, Loc.1: CT v.s. CK")))
s3_rm2 <- s3_rm2 + ggtitle(expression(paste("Sampling event 3, Loc.2: CT v.s. CK")))
s3_rm3 <- s3_rm3 + ggtitle(expression(paste("Sampling event 3, Loc.3: CT v.s. CK")))
s4_rm1 <- s4_rm1 + ggtitle(expression(paste("Sampling event 4, Loc.1: CT v.s. CK")))
s4_rm2 <- s4_rm2 + ggtitle(expression(paste("Sampling event 4, Loc.2: CT v.s. CK")))
s4_rm3 <- s4_rm3 + ggtitle(expression(paste("Sampling event 4, Loc.3: CT v.s. CK")))

# Summarize plot
pn_all <- list(s1_pn1, s1_pn2, s1_pn3, s2_pn1, s2_pn2, s2_pn3,
               s3_pn1, s3_pn2, s3_pn3, s4_pn1, s4_pn2, s4_pn3)
rm_all <- list(s1_rm1, s1_rm2, s1_rm3, s2_rm1, s2_rm2, s2_rm3,
               s3_rm1, s3_rm2, s3_rm3, s4_rm1, s4_rm2, s4_rm3)


# <----------------------------------------------------> #
#                     Save figures
# <----------------------------------------------------> #
#fig_output <- "00_FigRaw"
saveRDS(pn_all, sprintf("%s/Fig_RNAmaplot_PN.obj", fig_output))
saveRDS(rm_all, sprintf("%s/Fig_RNAmaplot_RM.obj", fig_output))



# <----------------------------------------------------> #
#               Save session information
# <----------------------------------------------------> #
#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/SessionInfo_07_RNA2019a_%s.txt", substr(Sys.time(), 1, 10)))

