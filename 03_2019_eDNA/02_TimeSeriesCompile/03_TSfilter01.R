####
#### CER Rice 2019 eDNA study
#### No.3 Filtering time series based on DNA concentration
####

# Load workspace
load("02_TSfilterPrepOut/02_TSfilterPrepOut.RData")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder03 <- "03_TSfilter01Out"
dir.create(output_folder03)

# Load library and functions
library(phyloseq); packageVersion("phyloseq") # 1.32.0, 2021.1.14

# Filter time series based on DNA copy numbers
ll_pro; ll_fun; ll_inv; ll_euk

# Exclude taxa of which average DNA copy number is lower than the lower limit
lower_limit_coef = 100
pro_cond1 <- taxa_sums(ps_pro_sample1) > nrow(sample_data(ps_pro_sample1))*ll_pro/lower_limit_coef
fun_cond1 <- taxa_sums(ps_fun_sample1) > nrow(sample_data(ps_fun_sample1))*ll_fun/lower_limit_coef
inv_cond1 <- taxa_sums(ps_inv_sample1) > nrow(sample_data(ps_inv_sample1))*ll_inv/lower_limit_coef
euk_cond1 <- taxa_sums(ps_euk_sample1) > nrow(sample_data(ps_euk_sample1))*ll_euk/lower_limit_coef

# Exclude taxa of which appearence is less than 10 times
pro_cond2 <- colSums(otu_table(ps_pro_sample1) > 0) > 10
fun_cond2 <- colSums(otu_table(ps_fun_sample1) > 0) > 10
inv_cond2 <- colSums(otu_table(ps_inv_sample1) > 0) > 10
euk_cond2 <- colSums(otu_table(ps_euk_sample1) > 0) > 10

# Include Pythiales and Chironomus kiiensis because they are manipulated taxa.
pro_cond3 <- c(tax_table(ps_pro_sample1)[,"order"] == "Pythiales") | c(tax_table(ps_pro_sample1)[,"species"] == "Chironomus kiiensis")
fun_cond3 <- c(tax_table(ps_fun_sample1)[,"order"] == "Pythiales") | c(tax_table(ps_fun_sample1)[,"species"] == "Chironomus kiiensis")
inv_cond3 <- c(tax_table(ps_inv_sample1)[,"order"] == "Pythiales") | c(tax_table(ps_inv_sample1)[,"species"] == "Chironomus kiiensis")
euk_cond3 <- c(tax_table(ps_euk_sample1)[,"order"] == "Pythiales") | c(tax_table(ps_euk_sample1)[,"species"] == "Chironomus kiiensis")

sum(pro_cond1); sum(fun_cond1); sum(inv_cond1); sum(euk_cond1) # 1146; 986; 556; 306
sum(pro_cond2); sum(fun_cond2); sum(inv_cond2); sum(euk_cond2) # 736; 238; 222; 271
sum(pro_cond3); sum(fun_cond3); sum(inv_cond3); sum(euk_cond3) # 7; 6; 37; 2
sum((pro_cond1 & pro_cond2) | pro_cond3,
    (fun_cond1 & fun_cond2) | fun_cond3,
    (inv_cond1 & inv_cond2) | inv_cond3,
    (euk_cond1 & euk_cond2) | euk_cond3) # 1467
sum(tax_n_pro, tax_n_fun, tax_n_inv, tax_n_euk) # 6891

ps_pro_sample2 <- prune_taxa((pro_cond1 & pro_cond2) | pro_cond3, ps_pro_sample1)
ps_fun_sample2 <- prune_taxa((fun_cond1 & fun_cond2) | fun_cond3, ps_fun_sample1)
ps_inv_sample2 <- prune_taxa((inv_cond1 & inv_cond2) | inv_cond3, ps_inv_sample1)
ps_euk_sample2 <- prune_taxa((euk_cond1 & euk_cond2) | euk_cond3, ps_euk_sample1)

# Save and output results
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder03, output_folder03))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder03, substr(Sys.time(), 1, 10)))
