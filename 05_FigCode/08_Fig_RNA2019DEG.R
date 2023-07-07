####
#### CER rice 2019
#### Rice RNA expression: Visualize specific genes
####

# Load workspace
result_folder01 <- "../04_2019_riceRNA"
load(sprintf("%s/05_DESeq2LocationOut/05_DESeq2LocationOut.RData", result_folder01))

# Create output folder
fig_output <- "00_FigRaw"

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

# Potential DEG data
data_folder <- result_folder01 %>% paste0("/06_SpecificGenesOut/")
deg_df <- read_csv(sprintf("%s/S2_DEGList.csv", data_folder))
deg_df$treatment <- factor(deg_df$treatment, levels = c("CT", "GN", "CK"))
deg_annot <- read_csv(sprintf("%s/S2_DEGannot.csv", data_folder))
## Select candicates for DEG examples (for Figure)
deg_example_cand <- deg_annot$Locus_ID %>% str_which("LOC_", negate = T) %>%
  deg_annot[.,] %>% arrange(padj)


# <----------------------------------------------------> #
#   Visualize each DEGs
# <----------------------------------------------------> #
# Figure 6: Example DEGs
genes_padj <- deg_example_cand %>% pull(Locus_ID)
genes_expr <- colSums(deg_df[,genes_padj]) %>% sort(decreasing = T) %>% names
genes_rank <- data.frame(locus_id = genes_padj, p_rank = 1:length(genes_padj),
                         expr_rank = 30 - rank(colSums(deg_df[,genes_padj]))) # Smaller is better
genes_rank$total_score <- genes_rank$p_rank + genes_rank$expr_rank
genes_rank$total_rank <- rank(genes_rank$total_score)
genes_rank <- genes_rank[order(genes_rank$total_rank),]
genes <- genes_rank$locus_id[c(1,4,5,6,8,15)]
# Select data
gene_i <- 1
deg_df2 <- deg_df %>% select(as.name("treatment"), genes[gene_i], as.name("plot")) %>% 
  rename_(value = as.name(genes[gene_i]))
## Visualize each DEG
set.seed(8181)
deg1 <-  deg_df2 %>% 
  ggplot(aes(x = treatment, y = value, group = treatment, color = treatment)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "white", color = "gray10") +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.8) +
  scale_y_log10() + scale_color_manual(values = c("gray70", "red3", "royalblue")) + theme(legend.position = "none") +
  ylab(genes[gene_i]) +
  NULL
deg2 <- deg1 %+% (deg_df2 %>% mutate(value = deg_df %>% pull(genes[2]))) + ylab(genes[2])
deg3 <- deg1 %+% (deg_df2 %>% mutate(value = deg_df %>% pull(genes[3]))) + ylab(genes[3])
deg4 <- deg1 %+% (deg_df2 %>% mutate(value = deg_df %>% pull(genes[4]))) + ylab(genes[4])
deg5 <- deg1 %+% (deg_df2 %>% mutate(value = deg_df %>% pull(genes[5]))) + ylab(genes[5])
deg6 <- deg1 %+% (deg_df2 %>% mutate(value = deg_df %>% pull(genes[6]))) + ylab(genes[6])


# <----------------------------------------------------> #
#                       Save figures
# <----------------------------------------------------> #
deg_fig_all <- list(deg1, deg2, deg3, deg4, deg5, deg6)
saveRDS(deg_fig_all, sprintf("%s/Fig_RNA2019_DEG.obj", fig_output))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/SessionInfo_05_DNA2019div_%s.txt", substr(Sys.time(), 1, 10)))




