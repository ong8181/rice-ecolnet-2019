####
#### CER rice 2019: Rice RNA expression analysis
#### No. 2: Load and merge count data
####

# create directory
output_folder02 <- "02_RSEMcountCompileOut"
dir.create(output_folder02)
dir.create("00_SessionInfo")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)

# Load library
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2021.1.18

# Import sample data
sample_sheet <- read.csv("sampledata/rna_sample_data_2run_combined.csv")

# Import RSEM results
rsem_res1 <- sort(list.files("01_RNAseq_CER2019_1stRun/final", pattern=".genes.results", full.names = T)) # Forward read files
rsem_res2 <- sort(list.files("01_RNAseq_CER2019_2ndRun/final", pattern=".genes.results", full.names = T)) # Forward read files
index_name1 <- sapply(str_split(sapply(str_split(rsem_res1, pattern = "/"), '[', 3), ".genes.results"), '[', 1)
index_name2 <- sapply(str_split(sapply(str_split(rsem_res2, pattern = "/"), '[', 3), ".genes.results"), '[', 1)
gene_info1 <- read.table("01_RNAseq_CER2019_1stRun/final/genes.list", header = F)
gene_info2 <- read.table("01_RNAseq_CER2019_2ndRun/final/genes.list", header = F)
## Merge two gene information if they are identical (delete each info object)
if(all(gene_info1 == gene_info2)) gene_info <- gene_info1; rm(gene_info1); rm(gene_info2)

# Analyzed variable
SELECTED_COL <- "expected_count"

# Read 1st run data
rna_count1 <- data.frame(col1 = read_tsv(rsem_res1[1], skip = 8, col_names = T, col_types = cols())[SELECTED_COL])
colnames(rna_count1) <- index_name1[1]
for(i in 2:length(rsem_res1)){
  start_row <- as.numeric(which("expected_count" == read_tsv(rsem_res1[i])[1:9,1]))
  df_tmp <- data.frame(read_tsv(rsem_res1[i], skip = start_row, col_names = T, col_types = cols())[SELECTED_COL])
  colnames(df_tmp) <- index_name1[i]
  rna_count1 <- cbind(rna_count1, df_tmp)
  rm(df_tmp)
}
dim(rna_count1); sum(!is.na(index_name1))

# Read 2nd run data
rna_count2 <- data.frame(col1 = read_tsv(rsem_res2[1], skip = 8, col_names = T, col_types = cols())[SELECTED_COL])
colnames(rna_count2) <- index_name2[1]
for(i in 2:length(rsem_res2)){
  if(any(index_name2[i] == na.omit(sample_sheet$index_2nd_run))){
    start_row <- as.numeric(which("expected_count" == read_tsv(rsem_res2[i])[1:9,1]))
    df_tmp <- data.frame(read_tsv(rsem_res2[i], skip = start_row, col_names = T, col_types = cols())[SELECTED_COL])
    colnames(df_tmp) <- index_name2[i]
    rna_count2 <- cbind(rna_count2, df_tmp)
    rm(df_tmp)
  }
}
dim(rna_count2); sum(!is.na(index_name2))

# Merge two
reseq_sample_id <- match(colnames(rna_count2), colnames(rna_count1))
rna_count <- rna_count1
rna_count[reseq_sample_id] <- rna_count1[reseq_sample_id] + rna_count2
dim(rna_count); rm(rna_count1); rm(rna_count2)
hist(colSums(rna_count))

# Check data
dim(rna_count); dim(gene_info)
sample_sheet$sum_expected_count <- colSums(rna_count)

# Save and output results
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder02, output_folder02))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder02, substr(Sys.time(), 1, 10)))
