####
#### CER Rice 2019 eDNA study
#### No.2 Filtering time series preparation (DNA concentration)
####

# Load workspace
load("01_PhyloseqImportOut/01_PhyloseqImportOut.RData")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder02 <- "02_TSfilterPrepOut"
dir.create(output_folder02)

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2021.1.14
library(phyloseq); packageVersion("phyloseq") # 1.32.0, 2021.1.14
library(cowplot); packageVersion("cowplot") # 1.1.0, 2021.1.14
theme_set(theme_cowplot())

# Extracting data as time series
ps_pro_sample1 <- prune_taxa(taxa_sums(ps_pro_sample) > 0, ps_pro_sample)
ps_fun_sample1 <- prune_taxa(taxa_sums(ps_fun_sample) > 0, ps_fun_sample)
ps_inv_sample1 <- prune_taxa(taxa_sums(ps_inv_sample) > 0, ps_inv_sample)
ps_euk_sample1 <- prune_taxa(taxa_sums(ps_euk_sample) > 0, ps_euk_sample)

# Check the number of taxa detected
(tax_n_pro <- ncol(otu_table(ps_pro_sample1))) # 2669 (14 Jan 2021)
(tax_n_fun <- ncol(otu_table(ps_fun_sample1))) # 2255 (14 Jan 2021)
(tax_n_inv <- ncol(otu_table(ps_inv_sample1))) # 1094 (14 Jan 2021)
(tax_n_euk <- ncol(otu_table(ps_euk_sample1))) # 873 (14 Jan 2021)
sum(tax_n_pro, tax_n_fun, tax_n_inv, tax_n_euk)
# 6891 taxa from 324 samples * 30-100 ml water!

# Set the lowest DNA copy numbers of standard DNAs to illustrate them
# PRO: max. 50000, 25000, 12500, 5000, 2500
# EUK: max. 20000, 10000, 5000, 2500, 500
# FUN: max. 1000, 500, 250, 100, 50
# INV: max. 400, 200, 100, 50, 10
ll_pro <- 2500
ll_fun <- 50
ll_inv <- 10
ll_euk <- 500

# Compile entropy information
lower_limit_coef <- 100
p1 <- ggplot(NULL, aes(x = taxa_sums(ps_pro_sample1)+0.5)) +
  geom_histogram(alpha = 0.8) + scale_x_log10() +
  geom_vline(xintercept = c(ll_pro/lower_limit_coef)*nrow(sample_data(ps_pro_sample1)), linetype = 2, size = 1, color = "red3") +
  xlab(expression(paste(Log[10], "(Sum of DNA copy estimate + 0.5)"))) + ylab("Count") +
  ggtitle(sprintf("Prokaryote 16S (No of ASV = %s)", tax_n_pro))

p2 <- ggplot(NULL, aes(x = taxa_sums(ps_fun_sample1)+0.5)) +
  geom_histogram(alpha = 0.8) + scale_x_log10() +
  geom_vline(xintercept = c(ll_fun/lower_limit_coef)*nrow(sample_data(ps_fun_sample1)), linetype = 2, size = 1, color = "red3") +
  xlab(expression(paste(Log[10], "(Sum of DNA copy estimate + 0.5)"))) + ylab("Count") +
  ggtitle(sprintf("Fungal ITS (No of ASV = %s)", tax_n_fun))

p3 <- ggplot(NULL, aes(x = taxa_sums(ps_inv_sample1)+0.5)) +
  geom_histogram(alpha = 0.8) + scale_x_log10() +
  geom_vline(xintercept = c(ll_inv/lower_limit_coef)*nrow(sample_data(ps_inv_sample1)), linetype = 2, size = 1, color = "red3") +
  xlab(expression(paste(Log[10], "(Sum of DNA copy estimate + 0.5)"))) + ylab("Count") +
  ggtitle(sprintf("Animal COI (No of ASV = %s)", tax_n_inv))

p4 <- ggplot(NULL, aes(x = taxa_sums(ps_euk_sample1)+0.5)) +
  geom_histogram(alpha = 0.8) + scale_x_log10() +
  geom_vline(xintercept = c(ll_euk/lower_limit_coef)*nrow(sample_data(ps_euk_sample1)), linetype = 2, size = 1, color = "red3") +
  xlab(expression(paste(Log[10], "(Sum of DNA copy estimate + 0.5)"))) + ylab("Count") +
  ggtitle(sprintf("Eukaryote 18S (No of ASV = %s)", tax_n_euk))


p_all <- plot_grid(p1, p2, p3, p4, ncol = 2, labels = "auto", align = "hv")
ggsave(sprintf("%s/DNAcopyDistribution.pdf", output_folder02), plot = p_all, width = 10, height = 10)

# Save and output results
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder02, output_folder02))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder02, substr(Sys.time(), 1, 10)))
