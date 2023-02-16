####
#### CER rice 2019
#### Rice RNA expression: RNA data compile (as phyloseq object)
####

# Load workspace
load("02_RSEMcountCompileOut/02_RSEMcountCompileOut.RData")

# create directory
output_folder03 <- "03_RNAdata2PhyloseqOut"
dir.create(output_folder03)

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)

# Load library
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2019.12.9
library(lubridate); packageVersion("lubridate") # 1.7.9.2, 2021.1.18
library(phyloseq); packageVersion("phyloseq") # 1.34.0, 2021.1.18

# Load annotation data
load("seqdata_annot/180604_des")

# Check data structure
dim(rna_count)
dim(sample_sheet)
dim(gene_info)

# Compiling sample data and RNA data
#(Extracting only monitoring data)
sample_sheet$date <- ymd(sample_sheet$date)
sample_sheet$plot <- as.factor(sample_sheet$plot)
sample_sheet$sampling_event <- as.factor(sample_sheet$sampling_event)
sample_sheet$treatment <- as.factor(sample_sheet$treatment)

# Adjust rownames and colnames of data sheet
all(sample_sheet$index_1st_run == colnames(rna_count))
rownames(sample_sheet) <- colnames(rna_count) <- sample_sheet$Sample_Name2
dim(rna_count); dim(gene_info); dim(des)
rownames(rna_count) <- rownames(des)

# Extract RNA with "data" category
data_cat_id <- which(des$NormalizationGroup == "data")
des2 <- des[data_cat_id,] # dim(des2)
rna_count2 <- rna_count[data_cat_id,]

# Remove rare RNA (include all RNAs)
quant <- 0
abundant_rna_index <- which(rowSums(rna_count2) > quant)
rna_data_top <- rna_count2[abundant_rna_index,]
gene_info_top <- des2[abundant_rna_index,]

# Check dimension and revise column and row names
dim(sample_sheet); dim(rna_data_top); dim(gene_info_top)
all(rownames(sample_sheet) == colnames(rna_data_top))
all(rownames(rna_data_top) == rownames(gene_info_top))

# Import to phyloseq
ps_rna_top_all <- phyloseq(otu_table(rna_data_top, taxa_are_rows = TRUE),
                           sample_data(sample_sheet),
                           tax_table(as.matrix(gene_info_top)))
ps_rna_top <- subset_samples(ps_rna_top_all, sample_type == "sample")
ps_rna_top <- prune_taxa(taxa_sums(ps_rna_top) > 0, ps_rna_top)

# Output RNA expression table, sample data and taxa information
rrna_all0 <- data.frame(otu_table(ps_rna_top, taxa_are_rows = TRUE)) %>% rownames_to_column()
rrna_tax0 <- data.frame(tax_table(ps_rna_top)@.Data) %>% rownames_to_column()
rrna_dat0 <- data.frame(sample_data(ps_rna_top)) %>% rownames_to_column()

# Output important objects
saveRDS(ps_rna_top, file = sprintf("%s/phyloseq_RNA_top.obj", output_folder03))
write_csv(rrna_all0, file = sprintf("%s/RNA_otu_table.csv", output_folder03))
write_csv(rrna_tax0, file = sprintf("%s/RNA_gene_table.csv", output_folder03))
write_csv(rrna_dat0, file = sprintf("%s/RNA_sample_table.csv", output_folder03))

# Delete temporal objects
rm(rrna_all0); rm(rrna_tax0); rm(rrna_dat0)
rm(rna_count2); rm(des2)
rm(index_name1); rm(index_name2)
rm(rsem_res1); rm(rsem_res2)
rm(abundant_rna_index)
rm(ps_rna_top_all)

# Save and output results
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder03, output_folder03))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder03, substr(Sys.time(), 1, 10)))

