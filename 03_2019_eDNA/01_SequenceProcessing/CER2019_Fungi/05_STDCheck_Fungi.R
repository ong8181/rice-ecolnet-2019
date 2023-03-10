####
#### CER rice 2019 eDNA analysis
#### No.5 Check standard curve (Fungi)
####

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder05 <- "05_STDCheck_FungiOut"
dir.create(output_folder05)

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2021.1.13
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.1.13
library(ggsci); packageVersion("ggsci") # 2.9, 2021.1.13
theme_set(theme_cowplot())

# Load helper functions
source("../functions/F01_HelperFunctions.R")

# Load workspace and data
load("01_CERrice2019_FungiSeqOut/01_CERrice2019_FungiSeqOut.RData")
tax_claident <- read.csv("04_TaxaSTDcombine_FungiOut/claident_tax_revise.csv")
sample_sheet <- read.csv("../sampledata_all/20200323_CMR-011_SampleSheet_FUN.csv")

# Set parameters
# PRO: max. 50000, 25000, 12500, 5000, 2500
# EUK: max. 20000, 10000, 5000, 2500, 500
# FUN: max. 1000, 500, 250, 100, 50
# INV: max. 400, 200, 100, 50, 10
std_copy_n <- c(1000, 500, 250, 100, 50)
col_name_level <- "Family"
seq_summary_fig_file <- sprintf("%s/Seq_Summary.pdf", output_folder05)
std_seq_head <- "STD_ITSasco"

# DADA2 summary objects
track

# Sample names check and sort if necessary
all(rownames(seqtab) == rownames(seqtab2)) # Should be TRUE
all(rownames(seqtab_nochim) == rownames(seqtab2)) # Should be TRUE
all(sample_sheet$Sample_Name2 == rownames(seqtab_nochim)) # Should be TRUE; For RMR-078, FALSE
nrow(sample_sheet) == nrow(seqtab_nochim) # Should be TRUE; For CMR-011_Fungi, FALSE
all(order(rownames(seqtab_nochim)) == 1:nrow(seqtab_nochim)) # Order check, should be TRUE
all(order(sample_sheet$Sample_Name2) == 1:nrow(sample_sheet)) # Order check, should be TRUE

rownames(sample_sheet) <- sample_sheet$Sample_Name2 # Replace rownames of sample_sheet


# Check whether there is 0 sequences samples
(zero_sample <- sample_sheet$Sample_Name2[is.na(match(sample_sheet$Sample_Name2, rownames(seqtab_nochim)))])
if(length(zero_sample) > 0){
  # If there is 0 sequence sample, run the code below
  # Generate dummy data frame
  seqtab_nochim_tmp <- matrix(0, ncol = ncol(seqtab_nochim), nrow = nrow(sample_sheet))
  rownames(seqtab_nochim_tmp) <- sample_sheet$Sample_Name2
  colnames(seqtab_nochim_tmp) <- colnames(seqtab_nochim)
  # Check dim
  dim(seqtab_nochim); dim(seqtab_nochim_tmp)
  # Add object
  seqtab_nochim_tmp[match(rownames(seqtab_nochim), rownames(seqtab_nochim_tmp)), 1:ncol(seqtab_nochim_tmp)] <-
    seqtab_nochim[,1:ncol(seqtab_nochim_tmp)]
  # Replace seq tables
  seqtab_nochim <- seqtab_nochim_tmp
}

# Extract standard sequeces
detected_std_name <- unique(tax_claident[which(substr(tax_claident[,"species"],1,nchar(std_seq_head)) == std_seq_head),"species"])
n_std_seq <- which(substr(tax_claident[,"species"],1,nchar(std_seq_head)) == std_seq_head)
std_table <- seqtab_nochim[,n_std_seq]
std_taxa <- tax_claident[n_std_seq, "species"]

# STD reads - copy number relationship
# Rename colnames
colnames(std_table) <- std_taxa
# Merge the same colnames
new_std_table <- data.frame(std_rank1 = MergeSTD(detected_std_name[1], std_data = std_table),
                            std_rank2 = MergeSTD(detected_std_name[2], std_data = std_table),
                            std_rank3 = MergeSTD(detected_std_name[3], std_data = std_table),
                            std_rank4 = MergeSTD(detected_std_name[4], std_data = std_table),
                            std_rank5 = MergeSTD(detected_std_name[5], std_data = std_table))
#nrow(new.std.table); length(sample.sheet$Sample_Name2)
new_std_table2 <- new_std_table[sample_sheet$sample_nc != "pcr_nc",]

# Linear regression
adj_r_fun <- function(x) summary(lm(as.numeric(x) ~ std_copy_n + 0))$adj.r.squared
lm_coef_fun <- function(x) summary(lm(as.numeric(x) ~ std_copy_n + 0))$coefficients[1]
r2_summary <- apply(new_std_table2, 1, adj_r_fun)
coef_summary <- apply(new_std_table2, 1, lm_coef_fun)
new_seqtab <- as.data.frame(seqtab_nochim[,-n_std_seq]) # Make seq table without standard DNA
new_seqtab2 <- new_seqtab[sample_sheet$sample_nc != "pcr_nc",] # Make seq table without negative control samples

## Visualize regression results
# 1. R2 value distribution
g1 <- ggplot(data.frame(values = r2_summary), aes(x = values)) +
  geom_histogram() + geom_hline(yintercept = 0, linetype = 2) + #+ xlim(0.9,1) + ylim(0,13)
  xlab(expression(paste("R"^{2}, " values"))) + ylab("Count") #+ scale_x_continuous(breaks=seq(0.8,1,by=0.05))

# 2. Slope distribution
g2 <- ggplot(data.frame(values = coef_summary), aes(x = values)) +
  geom_histogram() + geom_hline(yintercept = 0, linetype = 2) +
  xlab("Regression slope") + ylab("Count")

# 3. Regression examples
max_slope <- as.numeric(c(new_std_table2[which.max(coef_summary),], 0))
med_slope <- as.numeric(c(new_std_table2[which.min(abs(coef_summary - median(coef_summary))),], 0))
min_slope <- as.numeric(c(new_std_table2[which.min(coef_summary),], 0))
slope_summary <- pivot_longer(data.frame(copy =c(std_copy_n, 0),
                                         max_slope = max_slope,
                                         med_slope = med_slope,
                                         min_slope = min_slope), cols = -copy, names_to = "slope")
g3 <- ggplot(slope_summary, aes(x = copy, y = value, group = slope, colour = slope)) +
  geom_point(size = 3) + scale_color_manual(name = "Regression slope", values = c("red3", "darkred", "black")) +
  geom_abline(intercept = 0, slope = c(max(coef_summary), median(coef_summary), min(coef_summary)), color = c("red3", "darkred", "black")) +
  xlab(expression(paste("Standard DNA copies (", mu, l^{-1}, ")"))) + ylab("Standard DNA reads") +
  theme(legend.position = c(0.25, 0.75), legend.text = element_text(size = 7), legend.title = element_text(size = 8))

# 4. Read visualization
new_seqtab2$sample <- sample_sheet$Sample_Name2[sample_sheet$sample_nc != "pcr_nc"]
seqtab_plot <- pivot_longer(new_seqtab2, cols = -sample)
seqtab_plot <- seqtab_plot[seqtab_plot$value > 0,]
seqtab_plot$value2 <- seqtab_plot$value + 0.5
std_table_plot <- data.frame(sample = sample_sheet$Sample_Name2[sample_sheet$sample_nc != "pcr_nc"],
                             max = apply(new_std_table2, 1, max) + 0.5,
                             min = apply(new_std_table2, 1, min) + 0.5) # data.frame for overwrite standard DNA reads
# Box plot + jitter plot version
g4 <- ggplot(seqtab_plot, aes(x = sample, y = value2)) +
  geom_boxplot(shape = 16, alpha = 0.5) +
  geom_jitter(shape = 16, size = 1.5, alpha = 0.5, position = position_jitter(0.2)) +
  scale_y_log10() + ylab("Sequence reads") + xlab("Sample ID") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
  geom_segment(data = std_table_plot, aes(x = sample, xend = sample, y = min, yend = max),
               colour = "red3", size = 1.5, alpha = 0.5) + ggtitle("CMR-011, CERrice2019 Fungi seqs")

# Summarize visualization
top_row <- plot_grid(g1, g2, g3, ncol = 3, align = "h", labels = c("a", "b", "c"), rel_widths = c(1,1,1.3))
Fig_std <- plot_grid(top_row, g4, ncol = 1, labels = c("", "d"))

# save file
pdf(file = seq_summary_fig_file, width = 36, height = 20)
Fig_std; dev.off()

## Temporal workspace
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder05, output_folder05))


# Raw data check
new_seqtab_print <- seqtab_nochim
colnames(new_seqtab_print) <- sprintf("Taxa%05d", 1:ncol(seqtab_nochim))
tax_claident2 <- as.data.frame(tax_claident)
rownames(tax_claident2) <- sprintf("Taxa%05d", 1:ncol(seqtab_nochim))
write.csv(new_seqtab_print, sprintf("%s/allseq.csv", output_folder05))
write.csv(new_std_table, sprintf("%s/stdseq.csv", output_folder05))
write.csv(tax_claident2, sprintf("%s/tax.csv", output_folder05))

# Conversion of sample reads to calculated copy numbers
# Collect valid samples
valid_samples <- coef_summary[!(new_std_table2[,1] < 1 | r2_summary < 0.4)]
seqtab_valid <- seqtab_nochim[match(names(valid_samples), rownames(seqtab_nochim)),]
seqtab_conv0 <- seqtab_valid/valid_samples # Conversion
# Collect invalid samples (tentatively use raw values)
# (There are 1 speculative samples; S050)
(invalid_samples <- coef_summary[new_std_table2[,1] < 1 | r2_summary < 0.4])
seqtab_invalid <- seqtab_nochim[match(names(invalid_samples), rownames(seqtab_nochim)),]
if(is.null(dim(seqtab_invalid))) names(seqtab_invalid) <- names(invalid_samples)
# Collect PCR NC
seqtab_pcrnc <- seqtab_nochim[sample_sheet$sample_nc == "pcr_nc",]


# Combine valid, invalid & PCR NC samples
# Add 0 for PCR NC samples
seqtab_conv_w_std0 <- rbind(seqtab_conv0, seqtab_invalid, seqtab_pcrnc)
rownames(seqtab_conv_w_std0)
# Sort seqtab_conv_w_std
seqtab_conv_w_std0 <- seqtab_conv_w_std0[order(rownames(seqtab_conv_w_std0)),]
# Linear interpolations for invalid samples
before_invalid <- match(names(invalid_samples), rownames(seqtab_conv_w_std0)) - 1
invalid_index <- match(names(invalid_samples), rownames(seqtab_conv_w_std0))
after_invalid <- match(names(invalid_samples), rownames(seqtab_conv_w_std0)) + 1
seqtab_conv_w_std0[invalid_index,] <- (seqtab_conv_w_std0[before_invalid,] + seqtab_conv_w_std0[after_invalid,])/2

# Add 0 to PCR NC samples
seqtab_conv_w_std0[sample_sheet$sample_nc == "pcr_nc",] <- 0

# Coverted to "DNA copy numbers/water filetered"
water_vol <- sample_sheet$filt022_ml
water_vol[is.na(water_vol)] <- 1
seqtab_conv_w_std <- seqtab_conv_w_std0*100/water_vol # 1ul used in PCR * DNA extracts 100 ul/water volume (ml)

taxa_w_std <- tax_claident2
dim(seqtab_conv_w_std)
dim(taxa_w_std)

# Exclude standard DNAs from converted tables
seqtab_conv <- seqtab_conv_w_std[,-which(substr(tax_claident[,"species"],1,nchar(std_seq_head)) == std_seq_head)]
taxa_wo_std <- tax_claident2[-which(substr(tax_claident[,"species"],1,nchar(std_seq_head)) == std_seq_head),]
dim(seqtab_conv)
dim(taxa_wo_std)

# Add valid/invalid information to sample sheet
sample_sheet$std_validity <- "Valid"
sample_sheet$std_validity[invalid_index] <- "Invalid"
sample_sheet$std_validity[sample_sheet$sample_nc == "pcr_nc"] <- "PCR_NC"

# Save and output results
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder05, output_folder05))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder05, substr(Sys.time(), 1, 10)))
