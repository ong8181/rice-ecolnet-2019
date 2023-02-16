####
#### CERrice2017 All data analysis
#### No. 2.4 DNA cross-map Rice time series (with lagged CCM)
####

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder02 <- "02_UICRiceClimateDNAOut"

# Load library and functions
library(rUIC); packageVersion("rUIC") # 0.1.5, 2020.12.17
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2021.1.7

# Load workspace
load("01_CompileAllDataOut/01_CompileAllDataOut.RData")
uic_rice_clim <- readRDS(sprintf("%s/uic_rice_clim.obj", output_folder02))

# Determine the best embedding dimensions
TP_RANGE <- seq(2, -14, by = -1)
E_RANGE <- 0:14

# Prepare effect variable
y_std <- data.frame(y = as.numeric(scale(rice_all$gr)))

# Prepare conditional variables
temp_tp <- uic_rice_clim[uic_rice_clim$cause_var == "temp_mean", "te"] %>%
  which.max() %>% uic_rice_clim[uic_rice_clim$cause_var == "temp_mean",][., "tp"] %>% as.integer()
light_tp <- uic_rice_clim[uic_rice_clim$cause_var == "light_mean", "te"] %>%
  which.max() %>% uic_rice_clim[uic_rice_clim$cause_var == "light_mean",][., "tp"] %>% as.integer()
relhm_tp <- uic_rice_clim[uic_rice_clim$cause_var == "relhm_mean", "te"] %>%
  which.max() %>% uic_rice_clim[uic_rice_clim$cause_var == "relhm_mean",][., "tp"] %>% as.integer()
z1_std <- data.frame(z1 = dplyr::lag(as.numeric(scale(clim_all$temp_mean)), n = -temp_tp - 1))
z2_std <- data.frame(z2 = dplyr::lag(as.numeric(scale(clim_all$light_mean)), n = -light_tp - 1))
z3_std <- data.frame(z3 = dplyr::lag(as.numeric(scale(clim_all$relhm_mean)), n = -relhm_tp - 1))

# Prepare output object
uic_edna_rice0 <- data.frame()
total_process <- length(edna_var_coln)# * length(TP_RANGE)
process_i <- 1

# UIC main loop
for (edna_i in edna_var_coln){
  # Record start time
  start_time <- proc.time()[3]
  
  # Set potential causal variable
  x_std <- data.frame(x = as.numeric(scale(edna_all[edna_i])))
  
  # Determine an optimal E
  simp_res <- rUIC::simplex(cbind(x_std, y_std, z1_std, z2_std, z3_std), lib_var = "x",
                            cond_var = c("y", "z1", "z2", "z3"),
                            lib = edna_lib, E = E_RANGE, tp = 1, tau = 1, 
                            Enull = "adaptive", n_boot = 2000, seed = 12345)
  Eopt <- with(simp_res, max(c(0, E[pval < 0.05])))
  
  # Calculate UIC
  uic_res_tmp <- rUIC::uic(cbind(x_std, y_std, z1_std, z2_std, z3_std), lib_var = "x",
                           tar_var = "y", cond_var = c("z1", "z2", "z3"),
                           E = Eopt + 1, tp = TP_RANGE, lib = edna_lib,
                           n_boot = 2000, seed = 1234) %>%
    cbind(data.frame(effect_var =  rep(colnames(edna_all)[edna_i], length(TP_RANGE)),
                     cause_var = rep("gr", length(TP_RANGE))), .,
          data.frame(cond_var1 = rep("temp_mean", length(TP_RANGE)),
                     cond_var2 = rep("light_mean", length(TP_RANGE)),
                     cond_var3 = rep("relhm_mean", length(TP_RANGE))))
  
  # Combine results
  uic_edna_rice0 <- rbind(uic_edna_rice0, uic_res_tmp)
  
  # Output process message
  elapsed_time <- round(proc.time()[3] - start_time, 2)
  message(paste("Process", process_i, "/", total_process, "finished:", elapsed_time, "sec elapsed"))
  process_i <- process_i + 1
}

# Delete temporal objects
rm(start_time); rm(elapsed_time)
rm(uic_res_tmp); rm(x_std); rm(y_std)
rm(z1_std); rm(z2_std); rm(z3_std)
rm(process_i); rm(total_process)
rm(Eopt); rm(simp_res)

# Rename object
uic_edna_rice <- uic_edna_rice0
rm(uic_edna_rice0)

# Save workspace and objects
saveRDS(uic_edna_rice, sprintf("%s/uic_edna_rice.obj", output_folder02))
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s_4.RData", output_folder02, output_folder02))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_4_SessionInfo_%s.txt", output_folder02, substr(Sys.time(), 1, 10)))
