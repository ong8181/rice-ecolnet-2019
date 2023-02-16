####
#### CER rice 2019 eDNA study
#### FigCode: Pattern visualization after data filtering
####

# Load workspace
result_folder01 <- "../03_2019_eDNA"
load(sprintf("%s/02_TimeSeriesCompile/05_TSfilter02Out/05_TSfilter02Out.RData", result_folder01))

# Load library and functions
library(tsnemicrobiota); packageVersion("tsnemicrobiota") # 0.1.0, 2021.1.14
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2022.4.20
library(lubridate); packageVersion("lubridate") #1.8.0, 2022.4.20
library(phyloseq); packageVersion("phyloseq") # 1.38.0, 2022.4.20
library(cowplot); packageVersion("cowplot") # 1.1.1, 2022.4.20
library(ggsci); packageVersion("ggsci") #2.9, 2020.1.14
theme_set(theme_cowplot())

# Load helper functions
source(sprintf("%s/03_VisualizePatterns/functions/F03_FigHelperFunctions.R", result_folder01))

# Create output folder
fig_output <- "00_FigRaw"
#dir.create(fig_output)

# Revise the treatment names (2023.2.7)
treatment_ori <- c("CT", "PN", "RM")
treatment_rev <- c("CT", "GN", "CK")
sample_data(ps_filt) <- sample_data(ps_filt) %>% data.frame %>% 
  mutate(treatment = factor(treatment, levels = treatment_ori, labels = treatment_rev)) %>% sample_data


# <----------------------------------------------------> #
#             Figures for general patterns
# <----------------------------------------------------> #
# Superkingdom visualization
ps_filt_m1 <- speedyseq::psmelt(ps_filt)
ps_filt_m1$date <- ymd(ps_filt_m1$date)
ps_filt_m2 <- ps_filt_m1 %>% group_by(date, superkingdom) %>% summarize(Abundance = sum(Abundance))
colnames(ps_filt_m2) <- c("date", "superkingdom", "dna_conc")
ps_filt_m2[ps_filt_m2$superkingdom == "", "superkingdom"] <- "Undetermined"

f1 <- ggplot(ps_filt_m2, aes(x = date, y = dna_conc/9, group = superkingdom, fill = superkingdom)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_startrek() + xlab(NULL) + ylab("DNA (copies/ml water)") +
  NULL

# Diversity visualization
ps_filt_div <- ps_filt
ps_filt_div <- transform_sample_counts(ps_filt_div, ceiling)
sample_data(ps_filt_div)$date <- ymd(sample_data(ps_filt_div)$date)
f2 <- plot_richness(ps_filt_div, x = "date", color = "treatment", measures = c("Observed"))
f2$layers <- f2$layers[-1]
f2 <- f2 + geom_line(aes(group = plot)) +
  scale_color_manual(values = c("gray70", "red3", "royalblue")) +
  geom_point(size = 0.5) + ylab("No. of ASV") + xlab(NULL) + ggtitle("ASV diversity")

ps_filt_div_undet <- subset_taxa(ps_filt_div, superkingdom == "")
ps_filt_div_pro <- subset_taxa(ps_filt_div, superkingdom == "Bacteria" | superkingdom == "Archaea")
ps_filt_div_fun <- subset_taxa(ps_filt_div, kingdom == "Fungi")
ps_filt_div_euk <- subset_taxa(ps_filt_div, superkingdom == "Eukaryota" & kingdom != "Fungi")

f2_1 <- plot_richness(ps_filt_div_undet, x = "date", color = "treatment", measures = c("Observed"))
f2_1$layers <- f2_1$layers[-1]
f2_1 <- f2_1 + geom_line(aes(group = plot)) +
  scale_color_manual(values = c("gray70", "red3", "royalblue")) +
  geom_point(size = 0.5) +
  ylab("No. of ASV") +
  xlab(NULL) + ggtitle("ASV diversity (Undetermined)")

f2_2 <- plot_richness(ps_filt_div_pro, x = "date", color = "treatment", measures = c("Observed"))
f2_2$layers <- f2_2$layers[-1]
f2_2 <- f2_2 + geom_line(aes(group = plot)) +
  scale_color_manual(values = c("gray70", "red3", "royalblue")) +
  geom_point(size = 0.5) +
  ylab("No. of ASV") +
  xlab(NULL) + ggtitle("ASV diversity (Prokaryote)")

f2_3 <- plot_richness(ps_filt_div_fun, x = "date", color = "treatment", measures = c("Observed"))
f2_3$layers <- f2_3$layers[-1]
f2_3 <- f2_3 + geom_line(aes(group = plot)) +
  scale_color_manual(values = c("gray70", "red3", "royalblue")) +
  geom_point(size = 0.5) +
  ylab("No. of ASV") +
  xlab(NULL) + ggtitle("ASV diversity (Fungi)")

f2_4 <- plot_richness(ps_filt_div_euk, x = "date", color = "treatment", measures = c("Observed"))
f2_4$layers <- f2_4$layers[-1]
f2_4 <- f2_4 + geom_line(aes(group = plot)) +
  scale_color_manual(values = c("gray70", "red3", "royalblue")) +
  geom_point(size = 0.5) +
  ylab("No. of ASV") +
  xlab(NULL) + ggtitle("ASV diversity (Other Eukaryote)")


# <----------------------------------------------------> #
#             Figures for target taxa
# <----------------------------------------------------> #
# --------------------- Extract Pythium eDNAs
target_id1 <- c(tax_table(ps_filt)[,"order"] == "Pythiales")
ps_target1 <- prune_taxa(target_id1, ps_filt)
ps_m1 <- speedyseq::psmelt(ps_target1)
ps_m1$date <- ymd(ps_m1$date)
ps_m2 <- ps_m1 %>% group_by(date, OTU, plot, miseq_run) %>% summarize(Abundance = sum(Abundance))
colnames(ps_m2) <- c("date", "ASV", "plot", "miseq_run", "dna_conc")

## Figures
s1 <- ggplot(ps_m2, aes(x = date, y = log10(dna_conc+1), group = ASV, fill = ASV)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_igv() + xlab(NULL) +
  ylab(expression(paste(Log[10], "(DNA + 1) (copies/ml water)"))) +
  facet_wrap(~plot) + geom_hline(yintercept = 0, linetype = 2) +
  xlim(ymd("2019-06-18"), ymd("2019-07-12")) + ggtitle("Pythiales") +
  NULL
s2 <- ggplot(ps_m2, aes(x = date, y = log10(dna_conc + 1))) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_igv() + xlab(NULL) +
  ylab(expression(paste(Log[10], "(DNA + 1) (copies/ml water)"))) +
  facet_wrap(~plot) + geom_hline(yintercept = 0, linetype = 2) +
  xlim(ymd("2019-06-18"), ymd("2019-07-12")) + ggtitle("Pythiales") +
  NULL
s3 <- ps_m2 %>% filter(ASV == "Fungi_Taxa00886" | ASV == "Fungi_Taxa01950") %>% 
  ggplot(aes(x = date, y = log10(dna_conc+1), group = ASV, fill = ASV)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_igv() + xlab(NULL) +
  ylab(expression(paste(Log[10], "(DNA + 1) (copies/ml water)"))) +
  facet_wrap(~plot) + geom_hline(yintercept = 0, linetype = 2) +
  ggtitle("Fungi_Taxa00886 + Fungi_Taxa01950") +
  NULL
## Boxplot + Jitterplot for Pythium
ps_pythium <- ps_m2 %>% group_by(date, plot) %>%
  summarise(total_dna_conc = sum(dna_conc)) %>% 
  filter(date > ymd("2019-06-17") & date < ymd("2019-07-13")) %>% 
  mutate(treatment = recode(plot,
                            "1" = "CT", "2" = "GN", "3" = "CK",
                            "4" = "CT", "5" = "GN", "6" = "CK",
                            "7" = "CT", "8" = "GN", "9" = "CK", .default = plot),
         before_after = case_when(date <= ymd("2019-06-24") ~ "before",
                                  date > ymd("2019-06-24") ~ "after")) #%>% 
  #mutate(before_after = fct_relevel(before_after, c("before", "after")))
ps_pythium$treatment <- factor(ps_pythium$treatment, levels = treatment_rev)
ps_pythium$before_after <- factor(ps_pythium$before_after, levels = c("before", "after"))
s1_1 <- ps_pythium %>%
  mutate(total_dna_conc = if_else(total_dna_conc > 0.5, total_dna_conc, 0.5)) %>% 
  ggplot(aes(x = before_after, y = total_dna_conc, color = treatment)) +
  geom_boxplot(outlier.color = NA, outlier.size = 0, outlier.shape = NA, color = "gray20") +
  geom_jitter(width = 0.25, height = 0, alpha = 0.8) +
  scale_color_manual(values = c("gray20", "red3", "royalblue")) +
  facet_wrap(~ treatment) + panel_border() +
  xlab("Sampling before or after treatment") +
  ylab("Total DNA (copies/ml water)") +
  scale_y_log10(limits = c(0.5, 20000)) +
  NULL
s1_2 <- ps_m2 %>% #filter(dna_conc > 0) %>% 
  ggplot(aes(x = date, y = log10(dna_conc+1), group = ASV, fill = ASV)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_igv() + xlab(NULL) +
  ylab(expression(paste(Log[10], "(DNA + 1) (copies/ml water)"))) +
  facet_wrap(~ miseq_run + plot, ncol = 9) + geom_hline(yintercept = 0, linetype = 2) +
  xlim(ymd("2019-06-18"), ymd("2019-07-12")) + ggtitle("Pythiales") +
  NULL


# ---------------------------------------------------------------
# --------------------- Extract Chironomus kiiensis eDNAs
target_id2 <- c(tax_table(ps_filt)[,"species"] == "Chironomus kiiensis")
ps_target2 <- prune_taxa(target_id2, ps_filt)
ps_m3 <- speedyseq::psmelt(ps_target2)
ps_m3$date <- ymd(ps_m3$date)
ps_m4 <- ps_m3 %>% group_by(date, OTU, plot) %>% summarize(Abundance = sum(Abundance))
colnames(ps_m4) <- c("date", "ASV", "plot", "dna_conc")

# Figures
s3 <- ggplot(ps_m4, aes(x = date, y = log10(dna_conc + 1), group = ASV, fill = ASV)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_igv() + xlab(NULL) +
  ylab(expression(paste(Log[10], "(DNA + 1) (copies/ml water)"))) +
  facet_wrap(~plot) + geom_hline(yintercept = 0, linetype = 2) +
  xlim(ymd("2019-06-18"), ymd("2019-07-12")) +
  ggtitle("Chironomus kiiensis") +
  NULL
s4 <- ggplot(ps_m4, aes(x = date, y = log10(dna_conc + 1))) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_igv() + xlab(NULL) +
  ylab(expression(paste(Log[10], "(DNA + 1) (copies/ml water)"))) +
  facet_wrap(~plot) + geom_hline(yintercept = 0, linetype = 2) +
  xlim(ymd("2019-06-18"), ymd("2019-07-12")) +
  ggtitle("Chironomus kiiensis") +
  NULL
## Boxplot + Jitterplot for Midge
ps_midge <- ps_m4 %>% group_by(date, plot) %>%
  summarise(total_dna_conc = sum(dna_conc)) %>% 
  filter(date > ymd("2019-06-17") & date < ymd("2019-07-13")) %>% 
  mutate(treatment = recode(plot,
                            "1" = "CT", "2" = "GN", "3" = "CK",
                            "4" = "CT", "5" = "GN", "6" = "CK",
                            "7" = "CT", "8" = "GN", "9" = "CK", .default = plot),
         before_after = case_when(date <= ymd("2019-06-24") ~ "before",
                                  date > ymd("2019-06-24") ~ "after"))
ps_midge$treatment <- factor(ps_midge$treatment, levels = treatment_rev)
ps_midge$before_after <- factor(ps_midge$before_after, levels = c("before", "after"))
s2_1 <- ps_midge %>%
  mutate(total_dna_conc = if_else(total_dna_conc > 0.5, total_dna_conc, 0.5)) %>% 
  ggplot(aes(x = before_after, y = total_dna_conc, color = treatment)) +
  geom_boxplot(outlier.color = NA, outlier.size = 0, outlier.shape = NA, color = "gray20") +
  geom_jitter(width = 0.25, height = 0, alpha = 0.8) +
  scale_color_manual(values = c("gray20", "red3", "royalblue")) +
  facet_wrap(~ treatment) + panel_border() +
  xlab("Sampling before or after treatment") +
  ylab("Total DNA (copies / ml water)") +
  scale_y_log10(limits = c(0.5, 20000)) +
  NULL


# <----------------------------------------------------> #
#              Figures for dimension reduction
# <----------------------------------------------------> #
# NMDS for overall composition
ps_filt_bray <- ordinate(ps_filt, "NMDS", "bray")
n1 <- plot_ordination(ps_filt, ps_filt_bray, color = "treatment", shape = "treatment") +
  geom_point(size = 2) + scale_color_startrek()

# NMDS for 2019/6/18 - 2019/7/12
sample_data(ps_filt)$date <- ymd(sample_data(ps_filt)$date)
ps_filt2 <- subset_samples(ps_filt, date >= ymd("2019-06-18") & date <= ymd("2019-07-12"))
ps_filt_bray2 <- ordinate(ps_filt2, "NMDS", "bray")
n2 <- plot_ordination(ps_filt2, ps_filt_bray2, color = "treatment", shape = "treatment") +
  geom_point(size = 2) + scale_color_startrek() + facet_wrap(~treatment)
n3 <- plot_ordination(ps_filt2, ps_filt_bray2, color = "date", shape = "treatment") +
  geom_point(size = 2) +
  scale_color_gradient(low = "royalblue", high = "red3") +
  facet_wrap(~plot)
n4 <- plot_ordination(ps_filt2, ps_filt_bray2, color = "plot", shape = "treatment") +
  geom_point(size = 2) + scale_color_igv() + facet_wrap(~plot)

# t-SNE
tsne_res <- tsne_phyloseq(ps_filt2, distance='bray', perplexity = 50, verbose = 0, rng_seed = 3901)
t1 <- plot_tsne_phyloseq(ps_filt2, tsne_res, color = 'date', shape = "plot") +
  geom_point(size=2) + facet_wrap(~ plot) + scale_shape_manual(values = rep(1,9))
t2 <- plot_tsne_phyloseq(ps_filt2, tsne_res, color = 'treatment', shape = "plot") +
  geom_point(size = 3) + scale_shape_manual(values = rep(16,9)) + scale_color_igv() +
  xlab("t-SNE 1") + ylab("t-SNE 2")


# Dimension reduction for 2019/6/24 - 2019/7/12
ps_filt3 <- subset_samples(ps_filt, date >= ymd("2019-06-29") & date <= ymd("2019-07-12"))
ps_filt_bray3 <- ordinate(ps_filt3, "NMDS", "bray")
n5 <- plot_ordination(ps_filt3, ps_filt_bray3, color = "date", shape = "treatment") +
  geom_point(size = 2) + scale_color_gradient(low = "royalblue", high = "red3") + facet_wrap(~plot)
n6 <- plot_ordination(ps_filt3, ps_filt_bray3, color = "treatment", shape = "treatment") +
  geom_point(size = 2) + scale_color_igv()

tsne_res2 <- tsne_phyloseq(ps_filt3, distance="bray", perplexity = 40, verbose = 0, rng_seed = 3901)
t3 <- plot_tsne_phyloseq(ps_filt3, tsne_res2, color = 'treatment', shape = "plot") +
  geom_point(size = 3) + scale_shape_manual(values = 1:9) + scale_color_igv() +
  xlab("t-SNE 1") + ylab("t-SNE 2") +
  NULL
t4 <- plot_tsne_phyloseq(ps_filt3, tsne_res2, color = 'treatment') +
  geom_point(size = 2, alpha = 0.5, shape = 16) + scale_color_manual(values = c("gray50", "red3", "royalblue")) +
  xlab("t-SNE 1") + ylab("t-SNE 2") +
  NULL


# <----------------------------------------------------> #
#                       Save figures
# <----------------------------------------------------> #
general_fig_all <- list(f1, f2, f2_1, f2_2, f2_3, f2_4)
target_fig_all <- list(s1, s2, s3, s4)
target_fig_jitter <- list(s1_1, s2_1)
nmds_fig_all <- list(n1, n2, n3, n4, n5, n6)
tsne_fig_all <- list(t1, t2, t3, t4)
saveRDS(general_fig_all, sprintf("%s/Fig_DNA2019flt_general.obj", fig_output))
saveRDS(target_fig_all, sprintf("%s/Fig_DNA2019flt_target.obj", fig_output))
saveRDS(target_fig_jitter, sprintf("%s/Fig_DNA2019flt_target_jitter.obj", fig_output))
saveRDS(tsne_fig_all, sprintf("%s/Fig_DNA2019flt_tSNE.obj", fig_output))
ggsave(sprintf("%s/Pythiales.pdf", fig_output), plot = s1_2,
       width = 24, height = 10)

 #### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/SessionInfo_06_DNA2019flt_%s.txt", substr(Sys.time(), 1, 10)))
