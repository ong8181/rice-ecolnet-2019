####
#### CER rice 2019 eDNA analysis
#### No.6 CMR-010 Prokaryote: Read number visualization
####

# Load workspace and data
load("05_STDCheck_ProkOut/05_STDCheck_ProkOut.RData")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder06 <- "06_ReadNSummary_ProkOut"
dir.create(output_folder06)

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2021.1.13
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.1.13
library(ggsci); packageVersion("ggsci") # 2.9, 2021.1.13
library(scales); packageVersion("scales") #1.1.1, 2021.1.13
theme_set(theme_cowplot())

# Visualize data
fig_header <- "FigProk"

# Track data
track_df <- data.frame(track)
track_df$sample <- rownames(track_df)
track_df2 <- pivot_longer(track_df[,c(1:7,9)], cols = -sample, names_to = "variable")
v1 <- ggplot(track_df2, aes(x = variable, y = sample, fill = value)) +
  geom_tile(colour = "black", size = 0.001) + scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", midpoint = 10000, limits = c(0, 30000), oob = squish) +
  theme_gray(base_size = 8) + theme(axis.text.y = element_text(size = 1)) +
  geom_hline(yintercept = which(sample_sheet$sample_nc == "pcr_nc"), size = 0.1, alpha = 0.5)

ggsave(plot = v1, filename = sprintf("%s/%s_Track.pdf", output_folder06, fig_header), width = 10, height = 10)

# Sequence table
n_taxa <- 1:200
seqtab01 <- data.frame(seqtab_nochim)[,n_taxa]
colnames(seqtab01) <- sprintf("%s_%s_%s", tax_claident2$query[n_taxa], tax_claident2$phylum[n_taxa], tax_claident2$species[n_taxa])
seqtab01$sample <- rownames(seqtab01)
seqtab01_df <- pivot_longer(seqtab01, cols = -sample, names_to = "variable")
v2 <- ggplot(seqtab01_df, aes(x = variable, y = sample, fill = value)) +
  geom_tile() + scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", midpoint = 2000, oob = squish) +
  theme_gray(base_size = 8) + theme(axis.text.y = element_text(size = 1),
                                    axis.text.x = element_text(angle = 90, size = 3)) +
  geom_hline(yintercept = which(sample_sheet$sample_nc == "pcr_nc"), size = 0.1, alpha = 0.5)

ggsave(plot = v2, filename = sprintf("%s/%s_SeqtabNochim.pdf", output_folder06, fig_header), width = 10, height = 10)


seqtab02 <- data.frame(seqtab_conv)[,n_taxa]
colnames(seqtab02) <- sprintf("%s_%s_%s", tax_claident2$query[-n_std_seq][n_taxa], tax_claident2$phylum[-n_std_seq][n_taxa], tax_claident2$species[-n_std_seq][n_taxa])
seqtab02$sample <- rownames(seqtab02)
seqtab02_df <- pivot_longer(seqtab02, cols = -sample, names_to = "variable")
v3 <- ggplot(seqtab02_df, aes(x = variable, y = sample, fill = value)) +
  geom_tile() + scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", midpoint = 200000, oob = squish) +
  theme_gray(base_size = 8) + theme(axis.text.y = element_text(size = 1),
                                    axis.text.x = element_text(angle = 90, size = 3)) +
  geom_hline(yintercept = which(sample_sheet$sample_nc == "pcr_nc"), size = 0.1, alpha = 0.5)

ggsave(plot = v3, filename = sprintf("%s/%s_SeqtabConvNoSTD.pdf", output_folder06, fig_header), width = 10, height = 10)


## Summary Table
sample_summary <- sample_sheet

track2 <- data.frame(matrix(NA, nrow = nrow(sample_summary), ncol = ncol(track))) # 
rownames(track2) <- sample_summary$Sample_Name2
colnames(track2) <- colnames(track)
track2[match(rownames(track), rownames(track2)),] <- track

sample_summary$dada2_input <- track2$input
sample_summary$dada2_nochim <- track2$nochim
sample_summary$dada2_prop <- track2$prop
sample_summary$STD_all <- rowSums(new_std_table)
sample_summary$NonSTD_all <- rowSums(seqtab_nochim) - rowSums(new_std_table)
sample_summary$STD_prop <- sample_summary$STD_all / rowSums(seqtab_nochim)
sample_summary$STD_coef <- sample_summary$STD_r2 <- NA
sample_summary$STD_coef[match(names(coef_summary), sample_summary$Sample_Name2)] <- coef_summary
sample_summary$STD_r2[match(names(r2_summary), sample_summary$Sample_Name2)] <- r2_summary
sample_summary$dna_copy_sum <- rowSums(seqtab_conv)

write.csv(sample_summary, sprintf("%s/summary_sheet_Prok.csv", output_folder06))

# Repalce sample_sheet
sample_sheet <- sample_summary

# Save important objects
saveRDS(sample_sheet, sprintf("%s/sample_sheet_Prok.obj", output_folder06))
saveRDS(seqtab_nochim, sprintf("%s/seqtab_nochim_Prok.obj", output_folder06))
saveRDS(seqtab_conv, sprintf("%s/seqtab_conv_Prok.obj", output_folder06))
saveRDS(new_std_table, sprintf("%s/new_std_table_Prok.obj", output_folder06))
saveRDS(tax_claident, sprintf("%s/tax_claident_Prok.obj", output_folder06))
saveRDS(tax_claident2, sprintf("%s/tax_claident2_Prok.obj", output_folder06))
saveRDS(taxa_w_std, sprintf("%s/taxa_w_std_Prok.obj", output_folder06))
saveRDS(taxa_wo_std, sprintf("%s/taxa_wo_std_Prok.obj", output_folder06))

# Save and output results
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder06, output_folder06))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder06, substr(Sys.time(), 1, 10)))

