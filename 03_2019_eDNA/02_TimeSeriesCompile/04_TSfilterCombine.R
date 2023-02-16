####
#### CER Rice 2019 eDNA study
#### No.4 Combine all phyloseq objects
#### 2021.1.14 Ushio
#### R 4.0.2
####

# Load workspace
load("03_TSfilter01Out/03_TSfilter01Out.RData")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder04 <- "04_TSfilterCombineOut"
dir.create(output_folder04)

# Load library and functions
library(phyloseq); packageVersion("phyloseq") # 1.32.0, 2021.1.14

# Combine all phyloseq object
# Prepare sample sheet combined
sample_meta_all <- as.data.frame(sample_data(ps_pro_sample2)[,c(3,10:16)])
sample_meta_pro <- as.data.frame(sample_data(ps_pro_sample2)[,17:25])
sample_meta_fun <- as.data.frame(sample_data(ps_fun_sample2)[,17:25])
sample_meta_inv <- as.data.frame(sample_data(ps_inv_sample2)[,17:25])
sample_meta_euk <- as.data.frame(sample_data(ps_euk_sample2)[,17:25])

meta_colnames <- colnames(sample_sheet_prok)[17:25]
colnames(sample_meta_pro) <- sprintf("Prok_%s", meta_colnames)
colnames(sample_meta_fun) <- sprintf("Fungi_%s", meta_colnames)
colnames(sample_meta_inv) <- sprintf("Inv_%s", meta_colnames)
colnames(sample_meta_euk) <- sprintf("Euk_%s", meta_colnames)

sample_combined <- cbind(sample_meta_all, sample_meta_pro, sample_meta_fun, sample_meta_inv, sample_meta_euk)

# Pre-combined
ps_combined0 <- merge_phyloseq(ps_pro_sample2,
                               ps_fun_sample2,
                               ps_inv_sample2,
                               ps_euk_sample2)

# Re-order taxtable information
potential_taxcol_names <- c("query",
                            "superkingdom", "kingdom", "subkingdom",
                            "phylum", "class", "subclass", "infraclass",
                            "superorder", "order", "suborder", "infraorder", "parvorder",
                            "superfamily", "family","subfamily",
                            "tribe", "genus", "subgenus",
                            "species", "subspecies",
                            "seq", "seqlen", "miseq_run")
tax_table(ps_combined0) <- tax_table(ps_combined0)[,potential_taxcol_names]

# Re-merge phyloseq object (this process is done to keep sample information)
ps_combined <- phyloseq(otu_table(ps_combined0),
                        sample_data(sample_combined), # Keep all sample information
                        tax_table(ps_combined0))

# Save and output results
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder04, output_folder04))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder04, substr(Sys.time(), 1, 10)))
