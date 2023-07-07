####
#### CER rice 2019 eDNA study
#### FigCode: Overall diversity
####

# Load workspace
result_folder01 <- "../03_2019_eDNA"
load(sprintf("%s/02_TimeSeriesCompile/02_TSfilterPrepOut/02_TSfilterPrepOut.RData", result_folder01))

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2021.1.14
library(lubridate); packageVersion("lubridate") #1.7.9.2, 2021.1.14
library(phyloseq); packageVersion("phyloseq") # 1.32.0, 2021.1.14
library(cowplot); packageVersion("cowplot") # 1.1.0, 2021.1.14
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
sample_data(ps_pro_sample1) <- sample_data(ps_pro_sample1) %>% data.frame %>% 
  mutate(treatment = factor(treatment, levels = treatment_ori, labels = treatment_rev)) %>% sample_data
sample_data(ps_fun_sample1) <- sample_data(ps_fun_sample1) %>% data.frame %>% 
  mutate(treatment = factor(treatment, levels = treatment_ori, labels = treatment_rev)) %>% sample_data
sample_data(ps_inv_sample1) <- sample_data(ps_inv_sample1) %>% data.frame %>% 
  mutate(treatment = factor(treatment, levels = treatment_ori, labels = treatment_rev)) %>% sample_data
sample_data(ps_euk_sample1) <- sample_data(ps_euk_sample1) %>% data.frame %>% 
  mutate(treatment = factor(treatment, levels = treatment_ori, labels = treatment_rev)) %>% sample_data


# <----------------------------------------------------> #
# Climate time series
# <----------------------------------------------------> #
# Load climate data 2019
clim_all <- readRDS("../02_2019_Rice/data/data_climate.obj")
clim_fig <- na.omit(clim_all)
clim_fig$date <- ymd(clim_fig$date)
clim_fig <- clim_fig %>% filter(date >= "2019-05-20" & date <= "2019-09-20")
r1 <- ggplot(clim_fig, aes(x = date, y = temp_mean)) +
  geom_line() + geom_point(size = 0.5) + scale_color_startrek() +
  geom_line(aes(x = date, y = temp_max), linetype = 2, linewidth = 0.2) +
  geom_line(aes(x = date, y = temp_min), linetype = 2, linewidth = 0.2) +
  xlab(NULL) + ylab(expression(paste("Air temperature (", degree, "C)"))) +
  NULL


# <----------------------------------------------------> #
#                 Generate figures v1
# <----------------------------------------------------> #
#---------- CMR-010 Prokaryote ----------#
# Compile phyloseq object
ps_pro_sample3 <- transform_sample_counts(ps_pro_sample1, ceiling)
# Exclude non-prokayote DNAs
pro_tax_cond <- tax_table(ps_pro_sample3)[,"superkingdom"] == "Bacteria" | tax_table(ps_pro_sample3)[,"superkingdom"] == "Archaea"
pro_tax_extract <- taxa_names(tax_table(ps_pro_sample3)[pro_tax_cond])
ps_pro_sample3 <- prune_taxa(pro_tax_extract, ps_pro_sample3)
# Figures
sample_data(ps_pro_sample3)$date <- ymd(sample_data(ps_pro_sample3)$date)
f1 <- plot_richness(ps_pro_sample3, x = "date", color = "treatment", measures = c("Observed"))
f1$layers <- f1$layers[-1]
f1 <- f1 + geom_line(aes(group = plot)) +
  scale_color_startrek() + geom_point(size = 0.5) +
  ylab("No. of ASV") +
  xlab(NULL) + ggtitle("Overall ASV diversity (Prokaryote)")
#----------------------------------------#

#---------- CMR-011 Fungi ----------#
# Compile phyloseq object
ps_fun_sample3 <- transform_sample_counts(ps_fun_sample1, ceiling)
# Exclude non-fungal DNAs
fun_tax_cond <- tax_table(ps_fun_sample3)[,"kingdom"] == "Fungi"
fun_tax_extract <- taxa_names(tax_table(ps_fun_sample3)[fun_tax_cond])
ps_fun_sample3 <- prune_taxa(fun_tax_extract, ps_fun_sample3)
# Figures
sample_data(ps_fun_sample3)$date <- ymd(sample_data(ps_fun_sample3)$date)
f2 <- plot_richness(ps_fun_sample3, x = "date", color = "treatment", measures = c("Observed"))
f2$layers <- f2$layers[-1]
f2 <- f2 + geom_line(aes(group = plot)) +
  scale_color_startrek() + geom_point(size = 0.5) +
  ylab("No. of ASV") +
  xlab(NULL) + ggtitle("Overall ASV diversity (Fungi)")
#----------------------------------------#

#---------- CMR-011 Invertebrate ----------#
# Compile phyloseq object
ps_inv_sample3 <- transform_sample_counts(ps_inv_sample1, ceiling)
# Exclude non-fungal DNAs
inv_tax_cond <- tax_table(ps_inv_sample3)[,"kingdom"] == "Metazoa"
inv_tax_extract <- taxa_names(tax_table(ps_inv_sample3)[inv_tax_cond])
ps_inv_sample3 <- prune_taxa(inv_tax_extract, ps_inv_sample3)
# Figures
sample_data(ps_inv_sample3)$date <- ymd(sample_data(ps_inv_sample3)$date)
f3 <- plot_richness(ps_inv_sample3, x = "date", color = "treatment", measures = c("Observed"))
f3$layers <- f3$layers[-1]
f3 <- f3 + geom_line(aes(group = plot)) +
  scale_color_startrek() + geom_point(size = 0.5) +
  ylab("No. of ASV") +
  xlab(NULL) + ggtitle("Overall ASV diversity (Invertebrate)")
#----------------------------------------#

#---------- CMR-009 Eukaryote ----------#
# Compile phyloseq object
ps_euk_sample3 <- transform_sample_counts(ps_euk_sample1, ceiling)
# Exclude non-fungal DNAs
euk_tax_cond1 <- tax_table(ps_euk_sample3)[,"superkingdom"] == "Eukaryota"
euk_tax_cond2 <- tax_table(ps_euk_sample3)[,"kingdom"] != "Fungi"
euk_tax_cond3 <- tax_table(ps_euk_sample3)[,"kingdom"] != "Metazoa"
euk_tax_extract <- taxa_names(tax_table(ps_euk_sample3)[euk_tax_cond1 & euk_tax_cond2 & euk_tax_cond3])
ps_euk_sample3 <- prune_taxa(euk_tax_extract, ps_euk_sample3)
# Figures
sample_data(ps_euk_sample3)$date <- ymd(sample_data(ps_euk_sample3)$date)
f4 <- plot_richness(ps_euk_sample3, x = "date", color = "treatment", measures = c("Observed"))
f4$layers <- f4$layers[-1]
f4 <- f4 + geom_line(aes(group = plot)) +
  scale_color_startrek() + geom_point(size = 0.5) +
  ylab("No. of ASV") +
  xlab(NULL) + ggtitle("Overall ASV diversity (Invertebrate)")
#----------------------------------------#



# <----------------------------------------------------> #
#                 Generate figures v2
# <----------------------------------------------------> #
#---------- CMR-010 Prokaryote ----------#
# Compile phyloseq object
# Include all (i.e., non-prokayote DNAs)
ps_pro_m1 <- speedyseq::psmelt(taxa_name_summarize(ps_pro_sample1, "phylum", top_n = 12))
ps_pro_m1$date <- ymd(ps_pro_m1$date)
ps_pro_m2 <- ps_pro_m1 %>% group_by(date, phylum, plot) %>% summarize(Abundance = sum(Abundance))
ps_pro_m3 <- ps_pro_m1 %>% group_by(date, rep_tax, plot) %>% summarize(Abundance = sum(Abundance))
colnames(ps_pro_m2) <- c("date", "phylum", "plot", "dna_conc")
colnames(ps_pro_m3) <- c("date", "rep_tax", "plot", "dna_conc")
# Figures
g1_1 <- ggplot(ps_pro_m2, aes(x = date, y = dna_conc, group = phylum, fill = phylum)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_igv() + xlab(NULL) + ylab("DNA (copies/ml water)") + facet_wrap(~plot)
g1_2 <- ggplot(ps_pro_m3, aes(x = date, y = dna_conc, group = rep_tax, fill = rep_tax)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_igv() + xlab(NULL) + ylab("DNA (copies/ml water)") + facet_wrap(~plot)
#----------------------------------------#

#---------- CMR-011 Fungi ----------#
# Compile phyloseq object
# Include all (i.e., non-fungal DNAs)
ps_fun_m1 <- speedyseq::psmelt(taxa_name_summarize(ps_fun_sample1, "phylum", top_n = 9))
ps_fun_m1$date <- ymd(ps_fun_m1$date)
ps_fun_m2 <- ps_fun_m1 %>% group_by(date, phylum, plot) %>% summarize(Abundance = sum(Abundance))
ps_fun_m3 <- ps_fun_m1 %>% group_by(date, rep_tax, plot) %>% summarize(Abundance = sum(Abundance))
colnames(ps_fun_m2) <- c("date", "phylum", "plot", "dna_conc")
colnames(ps_fun_m3) <- c("date", "rep_tax", "plot", "dna_conc")
# Figures
g1_3 <- ggplot(ps_fun_m2, aes(x = date, y = dna_conc, group = phylum, fill = phylum)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_igv() + xlab(NULL) + ylab("DNA (copies/ml water)") + facet_wrap(~plot)
g1_4 <- ggplot(ps_fun_m3, aes(x = date, y = dna_conc, group = rep_tax, fill = rep_tax)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_igv() + xlab(NULL) + ylab("DNA (copies/ml water)") + facet_wrap(~plot)
#----------------------------------------#

#---------- CMR-011 Invertebrate ----------#
# Compile phyloseq object
# Include all (i.e., non-fungal DNAs)
ps_inv_m1 <- speedyseq::psmelt(taxa_name_summarize(ps_inv_sample1, "order", top_n = 12))
ps_inv_m1$date <- ymd(ps_inv_m1$date)
ps_inv_m2 <- ps_inv_m1 %>% group_by(date, phylum, plot) %>% summarize(Abundance = sum(Abundance))
ps_inv_m3 <- ps_inv_m1 %>% group_by(date, rep_tax, plot) %>% summarize(Abundance = sum(Abundance))
colnames(ps_inv_m2) <- c("date", "order", "plot", "dna_conc")
colnames(ps_inv_m3) <- c("date", "rep_tax", "plot", "dna_conc")
# Figures
g1_5 <- ggplot(ps_inv_m2, aes(x = date, y = dna_conc, group = order, fill = order)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_igv() + xlab(NULL) + ylab("DNA (copies/ml water)") + facet_wrap(~plot)
g1_6 <- ggplot(ps_inv_m3, aes(x = date, y = dna_conc, group = rep_tax, fill = rep_tax)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_igv() + xlab(NULL) + ylab("DNA (copies/ml water)") + facet_wrap(~plot)
#----------------------------------------#

#---------- CMR-009 Eukaryote ----------#
# Compile phyloseq object
# Include all (i.e., non-fungal DNAs)
ps_euk_m1 <- speedyseq::psmelt(taxa_name_summarize(ps_euk_sample1, "phylum", top_n = 12))
ps_euk_m1$date <- ymd(ps_euk_m1$date)
ps_euk_m2 <- ps_euk_m1 %>% group_by(date, phylum, plot) %>% summarize(Abundance = sum(Abundance))
ps_euk_m3 <- ps_euk_m1 %>% group_by(date, rep_tax, plot) %>% summarize(Abundance = sum(Abundance))
colnames(ps_euk_m2) <- c("date", "phylum", "plot", "dna_conc")
colnames(ps_euk_m3) <- c("date", "rep_tax", "plot", "dna_conc")
# Figures
g1_7 <- ggplot(ps_euk_m2, aes(x = date, y = dna_conc, group = phylum, fill = phylum)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_igv() + xlab(NULL) + ylab("DNA (copies/ml water)") + facet_wrap(~plot)
g1_8 <- ggplot(ps_euk_m3, aes(x = as.Date(date), y = dna_conc/5, group = rep_tax, fill = rep_tax)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_igv() + xlab(NULL) + ylab("DNA (copies/ml water)") + facet_wrap(~plot)
#----------------------------------------#



# <----------------------------------------------------> #
#                       Save figures
# <----------------------------------------------------> #
saveRDS(r1, sprintf("%s/Fig_DNA2019_clim.obj", fig_output))

fig_all_1 <- list(f1, f2, f3, f4)
fig_all_2 <- list(g1_1, g1_2, g1_3, g1_4, g1_5, g1_6, g1_7, g1_8)
saveRDS(fig_all_1, sprintf("%s/Fig_DNA2019div_1.obj", fig_output))
saveRDS(fig_all_2, sprintf("%s/Fig_DNA2019bar_1.obj", fig_output))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/SessionInfo_05_DNA2019div_%s.txt", substr(Sys.time(), 1, 10)))

