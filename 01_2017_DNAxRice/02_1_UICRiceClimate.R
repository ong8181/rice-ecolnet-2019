####
#### CERrice2017 All data analysis
#### No. 2 Rice cross-map Climate (with lagged UIC)
####

# Load workspace
load("01_CompileAllDataOut/01_CompileAllDataOut.RData")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder02 <- "02_UICRiceClimateDNAOut"
dir.create(output_folder02)

# Load library and functions
library(rUIC); packageVersion("rUIC") # 0.1.5, 2020.12.17
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2021.1.7
library(cowplot); packageVersion("cowplot") # 1.1.0, 2021.1.7
library(ggsci); packageVersion("ggsci") # 2.9, 2021.1.7
theme_set(theme_cowplot())

# Determine the best embedding dimensions
TP_RANGE <- seq(2, -14, by = -1)
E_RANGE <- 0:14

# Prepare effect variable
x_std <- data.frame(x = as.numeric(scale(rice_all$gr)))
# Prepare output object
uic_rice_clim0 <- data.frame()
total_process <- length(clim_var_coln[1:5]) #* length(TP_RANGE)
process_i <- 1

for (clim_i in clim_var_coln[1:5]){
  # Record start time
  start_time <- proc.time()[3]
  
  # Set potential causal variable
  y_std <- data.frame(y = as.numeric(scale(clim_all[clim_i])))
  
  # Determine an optimal E
  simp_res <- rUIC::simplex(cbind(x_std, y_std), lib_var = "x", cond_var = "y",
                            lib = rice_lib, E = E_RANGE, tp = 1, tau = 1, 
                            Enull = "adaptive", n_boot = 2000, seed = 1234)
  Eopt <- with(simp_res, max(c(0, E[pval < 0.05])))
  
  # Calculate UIC
  for (tp_i in TP_RANGE){
    uic_res_tmp <- rUIC::uic(cbind(x_std, y_std), lib_var = "x", tar_var = "y",
                             E = Eopt + 1, tp = tp_i, lib = rice_lib,
                             n_boot = 2000, seed = 1234) %>%
      cbind(data.frame(effect_var = "gr", cause_var = colnames(clim_all[clim_i])), .)
    uic_rice_clim0 <- rbind(uic_rice_clim0, uic_res_tmp)
  }
  
  # Output process message
  elapsed_time <- round(proc.time()[3] - start_time, 2)
  message(paste("Process", process_i, "/", total_process, "finished:", elapsed_time, "sec elapsed"))
  process_i <- process_i + 1
}
# Delete temporal objects
rm(start_time); rm(elapsed_time)
rm(uic_res_tmp); rm(x_std); rm(y_std)
rm(process_i); rm(total_process)
rm(Eopt); rm(simp_res)

# Visualize results
uic_rice_clim0$signif <- NaN
uic_rice_clim0$signif[uic_rice_clim0$pval <= 0.05] <- "P < 0.05"
uic_rice_clim0$signif[uic_rice_clim0$pval > 0.05] <- "N.S."
uic_rice_clim0$signif <- factor(uic_rice_clim0$signif, levels = c("P < 0.05", "N.S."))
selected_clim_var <- c("temp_mean", "relhm_mean", "light_mean", "actvp_mean", "satdf_mean")
uic_rice_clim <- uic_rice_clim0[!is.na(match(uic_rice_clim0$cause_var, selected_clim_var)),]
uic_rice_clim$cause_var <- factor(uic_rice_clim$cause_var,
                                  levels = c("temp_mean", "light_mean", "actvp_mean", "relhm_mean", "satdf_mean"))

g1 <- uic_rice_clim %>%
  ggplot(aes(x = tp, y = te, facet = cause_var)) + geom_line() +
  geom_point(aes(x = tp, y = te, color = signif)) +
  scale_color_manual(values = c("red3", "black"), name = NULL) +
  geom_hline(yintercept = 0, linetype = 2) + facet_wrap(~ cause_var) +
  ylab("TE") + xlab("UIC time lag (tp)") + ylim(-0.015, 0.065)

# Save workspace and ojcects
ggsave(sprintf("%s/uic_rice_clim.pdf", output_folder02), plot = g1, width = 10, height = 6)
saveRDS(uic_rice_clim, sprintf("%s/uic_rice_clim.obj", output_folder02))

save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s_1.RData", output_folder02, output_folder02))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_1_SessionInfo_%s.txt", output_folder02, substr(Sys.time(), 1, 10)))

