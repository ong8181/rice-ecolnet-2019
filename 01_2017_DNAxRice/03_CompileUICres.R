####
#### CERrice2017 All data analysis
#### No. 3 Compile UIC results
####

# Load workspace
load("01_CompileAllDataOut/01_CompileAllDataOut.RData")
output_folder02 <- "02_UICRiceClimateDNAOut"
output_folder03 <- "03_CompileUICresOut"
dir.create(output_folder03)

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2021.1.7
library(cowplot); packageVersion("cowplot") # 1.1.0, 2021.1.8
library(ggsci); packageVersion("ggsci") # 2.9, 2021.1.8
theme_set(theme_cowplot())

# Load UIC results
uic_rice_clim <- readRDS(sprintf("%s/uic_rice_clim.obj", output_folder02))
uic_clim_rice <- readRDS(sprintf("%s/uic_clim_rice.obj", output_folder02))
uic_rice_edna <- readRDS(sprintf("%s/uic_rice_edna.obj", output_folder02))
uic_edna_rice <- readRDS(sprintf("%s/uic_edna_rice.obj", output_folder02))

# Brief check by visualization
g1 <- ggplot(uic_rice_clim, aes(x = pval)) + geom_histogram() + ylim(0, 40) + xlab(expression(paste(italic(P), "-value"))) +
  ggtitle("UIC (Climate to Rice)") + geom_vline(xintercept = 0.05, linetype = 2, col = "red3")
g2 <- ggplot(uic_clim_rice, aes(x = pval)) + geom_histogram() + ylim(0, 40) + xlab(expression(paste(italic(P), "-value"))) +
  ggtitle("UIC (Rice to Climate)") + geom_vline(xintercept = 0.05, linetype = 2, col = "red3")
g3 <- ggplot(uic_rice_edna, aes(x = pval)) + geom_histogram() + ylim(0, 1500) + xlab(expression(paste(italic(P), "-value"))) +
  ggtitle("UIC (eDNA to Rice)") + geom_vline(xintercept = 0.05, linetype = 2, col = "red3")
g4 <- ggplot(uic_edna_rice, aes(x = pval)) + geom_histogram() + ylim(0, 1500) + xlab(expression(paste(italic(P), "-value"))) +
  ggtitle("UIC (Rice to eDNA)") + geom_vline(xintercept = 0.05, linetype = 2, col = "red3")
hist_all <- plot_grid(g1, g2, g3, g4, ncol = 2, align = "hv")
ggsave(sprintf("%s/HistPval.pdf", output_folder03), plot = hist_all,
       width = 10, height = 6)

# Extract causal variables
## Select maximum TE in significant tp values
## No.1: From Climate variables to Rice growth rate
uic_rice_clim$pair <- sprintf("effect_from_%s_to_%s", uic_rice_clim$cause_var,  uic_rice_clim$effect_var)
max_rice_clim <- uic_rice_clim[uic_rice_clim$pval < 0.05 & uic_rice_clim$tp < 1,] %>%
  group_by(pair) %>%
  summarise(max_te_id = which.max(te),
            effect_var = effect_var[max_te_id],
            cause_var = cause_var[max_te_id],
            E = E[max_te_id], tau = tau[max_te_id], tp = tp[max_te_id], nn = nn[max_te_id],
            n_lib = n_lib[max_te_id], n_pred = n_pred[max_te_id], rmse = rmse[max_te_id],
            te = max(te), pval = pval[max_te_id])

## No.2: From Rice growth rate to Climate variables
uic_clim_rice$pair <- sprintf("effect_from_%s_to_%s", uic_clim_rice$cause_var,  uic_clim_rice$effect_var)
max_clim_rice <- uic_clim_rice[uic_clim_rice$pval < 0.05 & uic_clim_rice$tp < 1,] %>%
  group_by(pair) %>%
  summarise(max_te_id = which.max(te),
            effect_var = effect_var[max_te_id],
            cause_var = cause_var[max_te_id],
            E = E[max_te_id], tau = tau[max_te_id], tp = tp[max_te_id], nn = nn[max_te_id],
            n_lib = n_lib[max_te_id], n_pred = n_pred[max_te_id], rmse = rmse[max_te_id],
            te = max(te), pval = pval[max_te_id])

## No.3: From Rice growth rate to eDNA
uic_rice_edna$pair <- sprintf("effect_from_%s_to_%s", uic_rice_edna$cause_var,  uic_rice_edna$effect_var)
max_rice_edna <- uic_rice_edna[uic_rice_edna$pval < 0.05 &
                                 uic_rice_edna$tp < 1 &
                                 TRUE,] %>%
  group_by(pair) %>%
  summarise(max_te_id = which.max(te),
            effect_var = effect_var[max_te_id],
            cause_var = cause_var[max_te_id],
            E = E[max_te_id], tau = tau[max_te_id], tp = tp[max_te_id], nn = nn[max_te_id],
            n_lib = n_lib[max_te_id], n_pred = n_pred[max_te_id], rmse = rmse[max_te_id],
            te = max(te), pval = pval[max_te_id])

## No.4: From eDNA to Rice growth rate
uic_edna_rice$pair <- sprintf("effect_from_%s_to_%s", uic_edna_rice$cause_var,  uic_edna_rice$effect_var)
max_edna_rice <- uic_edna_rice[uic_edna_rice$pval < 0.05 & uic_edna_rice$tp < 1,] %>%
  group_by(pair) %>%
  summarise(max_te_id = which.max(te),
            effect_var = effect_var[max_te_id],
            cause_var = cause_var[max_te_id],
            E = E[max_te_id], tau = tau[max_te_id], tp = tp[max_te_id], nn = nn[max_te_id],
            n_lib = n_lib[max_te_id], n_pred = n_pred[max_te_id], rmse = rmse[max_te_id],
            te = max(te), pval = pval[max_te_id])

# Extracted UIC objects
#max_rice_clim
#max_clim_rice
#max_rice_edna
#max_edna_rice

# Save workspace and objects
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder03, output_folder03))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder03, substr(Sys.time(), 1, 10)))

