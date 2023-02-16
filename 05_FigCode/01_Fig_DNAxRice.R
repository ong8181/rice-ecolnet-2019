####
#### CER eDNA study: DNAxRice
#### FigCode: Pattern visualization for all time series
####

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2022.4.19
library(gganimate); packageVersion("gganimate") # 1.0.7, 2021.3.23
library(lubridate); packageVersion("lubridate") #1.8.0, 2022.4.19
library(phyloseq); packageVersion("phyloseq") # 1.38.0, 2022.4.19
library(reshape2); packageVersion("reshape2") # 1.4.4, 2021.1.12
library(cowplot); packageVersion("cowplot") # 1.1.0, 2021.1.12
library(cols4all); packageVersion("cols4all") # 0.1, 2022.4.21
library(ggsci); packageVersion("ggsci") #2.9, 2018.8.10
theme_set(theme_cowplot())

# Create output folder
fig_output <- "00_FigRaw"
dir.create(fig_output)
dir.create("00_SessionInfo")

# Load workspace
result_folder01 <- "../01_2017_DNAxRice"
load(sprintf("%s/01_CompileAllDataOut/01_CompileAllDataOut.RData", result_folder01))


# <----------------------------------------------------> #
#                 Compile figures
# <----------------------------------------------------> #
# Combine all to make composite time series
#rice_all; clim_all; edna_all
# Rice time series
rice_fig <- na.omit(rice_all)
r1_1 <- ggplot(rice_fig, aes(x = date, y = gr, colour = plot, group = plot)) +
  geom_line() + geom_point(size = 0.5) + scale_color_startrek() +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab(NULL) + ylab("Rice growth rate (cm/day)") + ylim(-1, 7) +
  NULL

# Climate time series
clim_fig <- na.omit(clim_all)
r1_2 <- ggplot(clim_fig, aes(x = date, y = temp_mean, colour = plot, group = plot)) +
  geom_line() + geom_point(size = 0.5) + scale_color_startrek() +
  geom_line(aes(x = date, y = temp_max), linetype = 2, size = 0.2) +
  geom_line(aes(x = date, y = temp_min), linetype = 2, size = 0.2) +
  xlab(NULL) + ylab(expression(paste("Air temperature (", degree, "C)"))) +
  NULL

r1_3 <- ggplot(rice_fig, aes(x = date, y = height, colour = plot, group = plot)) +
  geom_line() + geom_point(size = 0.5) + scale_color_startrek() +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab(NULL) + ylab("Rice height (cm)") + ylim(0, 150) +
  NULL


# <----------------------------------------------------> #
#              Generate eDNA animations
# <----------------------------------------------------> #
# eDNA time series
# Extract taxa information
find_lowest_phylogeny <- function(x){
  lowest_id <- which(!(is.na(edna_tax[x,2:20]) | nchar(edna_tax[x,2:20]) < 1))
  if(length(lowest_id)<1){
    lowest_phylogeny_print <- sprintf("[%s:%s]", rownames(edna_tax)[x], "[Undetermined]")
  }else{
    lowest_phylogeny <- max(lowest_id)
    lowest_phylogeny_print <- sprintf("[%s:%s] %s",
                                      rownames(edna_tax)[x],
                                      colnames(edna_tax[x,lowest_phylogeny+1]),
                                      edna_tax[x,lowest_phylogeny+1])
  }
  return(lowest_phylogeny_print)
}

lowest_tax <- c(NULL)
for(i in 1:nrow(edna_tax)) lowest_tax[i] <- find_lowest_phylogeny(i)

# eDNA ggplot
edna_prep <- edna_all
colnames(edna_prep)[edna_var_coln] <- lowest_tax
top_cols1 <- 1:50 + 40 # Prokaryote
top_cols2 <- 501:550 + 40 # Fungi
top_cols3 <- 832:851 + 40 # Invertebrate
top_cols4 <- 978:(978+50) + 40 # Eukaryote
edna_melt1 <- na.omit(melt(edna_prep, id.vars = c("date", "plot"), measure.vars = top_cols1))
edna_melt2 <- na.omit(melt(edna_prep, id.vars = c("date", "plot"), measure.vars = top_cols2))
edna_melt3 <- na.omit(melt(edna_prep, id.vars = c("date", "plot"), measure.vars = top_cols3))
edna_melt4 <- na.omit(melt(edna_prep, id.vars = c("date", "plot"), measure.vars = top_cols4))
edna_melt1[edna_melt1$value < 1e+02,]$value <- 1e+02
edna_melt2[edna_melt2$value < 1,]$value <- 1
edna_melt3[edna_melt3$value < 1,]$value <- 1
edna_melt4[edna_melt4$value < 10,]$value <- 10

r1_4 <- ggplot(edna_melt1, aes(x = date, y = value + 0.5, group = plot, colour = plot)) +
  geom_point(size = 0.5) + geom_line() + scale_color_startrek() + scale_y_log10() +
  xlab(NULL) + ylab("Log(DNA copy + 0.5)") +
  labs(title = "{current_frame}") + transition_manual(variable) +
  NULL
r1_5 <- ggplot(edna_melt2, aes(x = date, y = value + 0.5, group = plot, colour = plot, frame = variable)) +
  geom_point(size = 0.5) + geom_line() + scale_color_startrek() + scale_y_log10() +
  xlab(NULL) + ylab("Log(DNA copy + 0.5)") +
  labs(title = "{current_frame}") + transition_manual(variable) +
  NULL
r1_6 <- ggplot(edna_melt3, aes(x = date, y = value + 0.5, group = plot, colour = plot, frame = variable)) +
  geom_point(size = 0.5) + geom_line() + scale_color_startrek() + scale_y_log10() +
  xlab(NULL) + ylab("Log(DNA copy + 0.5)") +
  labs(title = "{current_frame}") + transition_manual(variable) +
  NULL
r1_7 <- ggplot(edna_melt4, aes(x = date, y = value + 0.5, group = plot, colour = plot, frame = variable)) +
  geom_point(size = 0.5) + geom_line() + scale_color_startrek() + scale_y_log10() +
  xlab(NULL) + ylab("Log(DNA copy + 0.5)") +
  labs(title = "{current_frame}") + transition_manual(variable) +
  NULL


# <----------------------------------------------------> #
#                 Compile UIC figures
# <----------------------------------------------------> #
# Specify UIC result folder
result_folder02 <- "~/Work/Result/AllResult/2016_PRESTOwork/1_CER/6_RiceComAnalysis/02_UICRiceClimateDNAOut"
# Load UIC results
uic_rice_clim <- readRDS(sprintf("%s/uic_rice_clim.obj", result_folder02))
uic_clim_rice <- readRDS(sprintf("%s/uic_clim_rice.obj", result_folder02))
uic_rice_edna <- readRDS(sprintf("%s/uic_rice_edna.obj", result_folder02))
uic_edna_rice <- readRDS(sprintf("%s/uic_edna_rice.obj", result_folder02))
# Brief check by visualization
u1 <- ggplot(uic_rice_clim, aes(x = pval)) + geom_histogram() + ylim(0, 40) + xlab(expression(paste(italic(P), "-value"))) +
  ggtitle("UIC (Climate to Rice)") + geom_vline(xintercept = 0.05, linetype = 2, col = "red3")
u2 <- ggplot(uic_clim_rice, aes(x = pval)) + geom_histogram() + ylim(0, 40) + xlab(expression(paste(italic(P), "-value"))) +
  ggtitle("UIC (Rice to Climate)") + geom_vline(xintercept = 0.05, linetype = 2, col = "red3")
u3 <- ggplot(uic_rice_edna, aes(x = pval)) + geom_histogram() + ylim(0, 1500) + xlab(expression(paste(italic(P), "-value"))) +
  ggtitle("UIC (eDNA to Rice)") + geom_vline(xintercept = 0.05, linetype = 2, col = "red3")
u4 <- ggplot(uic_edna_rice, aes(x = pval)) + geom_histogram() + ylim(0, 1500) + xlab(expression(paste(italic(P), "-value"))) +
  ggtitle("UIC (Rice to eDNA)") + geom_vline(xintercept = 0.05, linetype = 2, col = "red3")

# Demonstration of UIC results
tmp_df <- rbind(uic_rice_clim[uic_rice_clim$cause_var == "temp_mean",c("cause_var", "te", "tp")],
                uic_clim_rice[uic_clim_rice$effect_var == "temp_mean",c("cause_var", "te", "tp")])
tmp_df$cause_var <- as.character(tmp_df$cause_var)
tmp_df$cause_var[tmp_df$cause_var == "temp_mean"] <- "Temperature to rice growth"
tmp_df$cause_var[tmp_df$cause_var == "gr"] <- "Rice growth to temperature"
tmp_df$cause_var <- factor(tmp_df$cause_var, levels = unique(tmp_df$cause_var)[1:2])
u5 <- ggplot(tmp_df, aes(x = tp, y = te)) + 
  geom_point() + geom_line() +
  facet_wrap(~ cause_var) + panel_border() +
  xlab("Time lag (day)") + ylab("Strength of causality\n(Information transfer)") +
  xlim(-14, 0) + geom_hline(yintercept = 0, linetype = 2)

# Specify UIC result compile folder
load(sprintf("%s/03_CompileUICresOut/03_CompileUICresOut.RData", result_folder01))


## Extract taxa with species-level identification
causal_spp1 <- intersect(max_rice_edna$cause_var,
                         rownames(edna_tax)[which(edna_tax[,"species"] != "")])
causal_spp2 <- sort(c(causal_spp1, "Fungi_Taxa00402"))
max_rice_edna_sub <- max_rice_edna %>%
  filter(cause_var %in% causal_spp2) %>%
  .[order(.$te, decreasing = F),]
max_rice_edna_sub$cause_var <- factor(max_rice_edna_sub$cause_var, levels = max_rice_edna_sub$cause_var)

u6 <- ggplot(max_rice_edna_sub, aes(y = factor(cause_var), x = te)) +
  geom_point(size = 1, alpha = 0.5) +
  theme(axis.text.y = element_text(size = 6)) +
  scale_x_log10() +
  ylab("Causal species to rice growth") +
  xlab("Influences to rice growth")

## Output taxa table
max_rice_edna_sub2 <- max_rice_edna_sub %>%
  bind_cols(edna_tax[as.character(.$cause_var),]) %>% arrange(superkingdom, phylum, te)
uncultured_id <- max_rice_edna_sub2$species %>% grep("uncultured", .)
max_rice_edna_sub2 <- max_rice_edna_sub2[-uncultured_id,]
max_rice_edna_sub2 <- read.csv("00_TaxData/max_rice_edna_sub2.csv") %>% 
  mutate(tax_name2 = paste0(tax_order, " ", tax_name)) %>% 
  arrange(tax_name2, -te) %>% 
  mutate(cause_var = factor(.$cause_var, levels = rev(.$cause_var))) %>% 
  mutate(tax_name = factor(.$tax_name, levels = unique(.$tax_name)))
u7 <- max_rice_edna_sub2 %>% 
  ggplot(aes(y = factor(cause_var), x = te, color = tax_name)) +
  geom_point(size = 1, alpha = 1) +
  theme(axis.text.y = element_text(size = 6)) +
  scale_color_manual(values = c4a("watlington", 12), name = "Taxa name") +
  scale_x_log10() +
  ylab("Causal species to rice growth") +
  xlab("Influences to rice growth")
#write.csv(max_rice_edna_sub2, "00_TaxData/causal_taxa_table.csv")



# <----------------------------------------------------> #
#                    Save figures
# <----------------------------------------------------> #
rice_plot_figs <- list(r1_1, r1_2, r1_3)
uic_figs <- list(u1, u2, u3, u4, u5, u6, u7)
saveRDS(rice_plot_figs, sprintf("%s/Fig_DNAxRice_RiceTS.obj", fig_output))
saveRDS(uic_figs, sprintf("%s/Fig_DNAxRice_UIChist.obj", fig_output))


# <----------------------------------------------------> #
#                 Save eDNA animations
# <----------------------------------------------------> #
# gganimate output
# save_file_name1 <- sprintf("%s/animation_pro_top50.gif", fig_output)
# save_file_name2 <- sprintf("%s/animation_fun_top50.gif", fig_output)
# save_file_name3 <- sprintf("%s/animation_inv_top50.gif", fig_output)
# save_file_name4 <- sprintf("%s/animation_euk_top50.gif", fig_output)
# anim_save(save_file_name1, r1_4, nframes = 50, fps = 5, height = 200, width = 500)
# anim_save(save_file_name2, r1_5, nframes = 50, fps = 5, height = 200, width = 500)
# anim_save(save_file_name3, r1_6, nframes = 50, fps = 5, height = 200, width = 500)
# anim_save(save_file_name4, r1_7, nframes = 50, fps = 5, height = 200, width = 500)

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/SessionInfo_01_DNAxRiceFig_%s.txt", substr(Sys.time(), 1, 10)))


