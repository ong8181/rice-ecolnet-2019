####
#### CERrice2017 All data analysis
#### No. 1 All data compile (Rice, Climate and eDNA data)
####

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder01 <- "01_CompileAllDataOut"
dir.create("00_SessionInfo")
dir.create(output_folder01)

# Load library and functions
library(phyloseq); packageVersion("phyloseq") # 1.22.3, 2018.12.20
library(lubridate); packageVersion("lubridate") # 1.7.4, 2018.12.20
library(rEDM); packageVersion("rEDM") # 0.7.5, 2019.4.3
library(rUIC); packageVersion("rUIC") # 0.1.5, 2020.12.16
source("functions/F01_HelperFunctions.R")
source("functions/BidirectSimplex_rEDM_0.7.4.R")

# Load main data (Rice, Climate and eDNA data)
rice_all0 <- readRDS("data_rice/rice_all.obj")
clim_all0 <- readRDS("data_climate/clm_d_mean.obj")
edna_all0 <- readRDS("data_edna/ps_comb_filt.obj")

# Load supplement data
clim_max <- readRDS("data_climate/clm_d_max.obj")
clim_min <- readRDS("data_climate/clm_d_min.obj")
clim_cer_mean <- readRDS("data_climate/cer_d_mean.obj")
clim_cer_max <- readRDS("data_climate/cer_d_max.obj")
clim_cer_min <- readRDS("data_climate/cer_d_min.obj")

# Adjust data format
clim_all0 <- cbind(clim_all0[,1:6], clim_all0[,17:21], clim_all0[,27:31], clim_all0[,7:11], clim_all0[,12:16])
clim_max <- cbind(clim_max[,1:6], clim_max[,17:21], clim_max[,12:16])
clim_min <- cbind(clim_min[,1:6], clim_min[,17:21], clim_min[,12:16])

# Set climate dataframe as a default
# Climate dataframe for each plot
date_for_df <- ymd(clim_all0[,"date"])
#all(colnames(clim_min) == colnames(clim_max))
p1_col <- substr(colnames(clim_all0), 1, 5) == "D0202"
p2_col <- substr(colnames(clim_all0), 1, 5) == "D0203"
p3_col <- substr(colnames(clim_all0), 1, 5) == "D0204"
p4_col <- substr(colnames(clim_all0), 1, 5) == "D0205"
p5_col <- substr(colnames(clim_all0), 1, 5) == "D0206"
p1_col2 <- substr(colnames(clim_max), 1, 5) == "D0202"
p2_col2 <- substr(colnames(clim_max), 1, 5) == "D0203"
p3_col2 <- substr(colnames(clim_max), 1, 5) == "D0204"
p4_col2 <- substr(colnames(clim_max), 1, 5) == "D0205"
p5_col2 <- substr(colnames(clim_max), 1, 5) == "D0206"

plot1_clim_df <- cbind(date_for_df, rep(1, nrow(clim_all0)), clim_all0[,p1_col], clim_max[,p1_col2], clim_min[,p1_col2])
plot2_clim_df <- cbind(date_for_df, rep(2, nrow(clim_all0)), clim_all0[,p2_col], clim_max[,p2_col2], clim_min[,p2_col2])
plot3_clim_df <- cbind(date_for_df, rep(3, nrow(clim_all0)), clim_all0[,p3_col], clim_max[,p3_col2], clim_min[,p3_col2])
plot4_clim_df <- cbind(date_for_df, rep(4, nrow(clim_all0)), clim_all0[,p4_col], clim_max[,p4_col2], clim_min[,p4_col2])
plot5_clim_df <- cbind(date_for_df, rep(5, nrow(clim_all0)), clim_all0[,p5_col], clim_max[,p5_col2], clim_min[,p5_col2])
clim_colnames <- c("date", "plot", "temp_mean", "relhm_mean", "satdf_mean", "actvp_mean", "light_mean", "temp_max", "relhm_max", "light_max", "temp_min", "relhm_min", "light_min")
colnames(plot1_clim_df) <- colnames(plot2_clim_df) <- colnames(plot3_clim_df) <- colnames(plot4_clim_df) <- colnames(plot5_clim_df) <- clim_colnames
plot1_clim_df$plot <- as.factor(plot1_clim_df$plot)
plot2_clim_df$plot <- as.factor(plot2_clim_df$plot)
plot3_clim_df$plot <- as.factor(plot3_clim_df$plot)
plot4_clim_df$plot <- as.factor(plot4_clim_df$plot)
plot5_clim_df$plot <- as.factor(plot5_clim_df$plot)

# Add cumulative climate variables
plot1_clim_df <- calculate_cum_vals(plot1_clim_df, "temp_mean", cum_range = 2:14, cum_func = "sum")
plot1_clim_df <- calculate_cum_vals(plot1_clim_df, "relhm_mean", cum_range = 2:14, cum_func = "mean")
plot1_clim_df <- calculate_cum_vals(plot1_clim_df, "satdf_mean", cum_range = 2:14, cum_func = "sum")
plot1_clim_df <- calculate_cum_vals(plot1_clim_df, "actvp_mean", cum_range = 2:14, cum_func = "sum")
plot1_clim_df <- calculate_cum_vals(plot1_clim_df, "light_mean", cum_range = 2:14, cum_func = "sum")

plot2_clim_df <- calculate_cum_vals(plot2_clim_df, "temp_mean", cum_range = 2:14, cum_func = "sum")
plot2_clim_df <- calculate_cum_vals(plot2_clim_df, "relhm_mean", cum_range = 2:14, cum_func = "mean")
plot2_clim_df <- calculate_cum_vals(plot2_clim_df, "satdf_mean", cum_range = 2:14, cum_func = "sum")
plot2_clim_df <- calculate_cum_vals(plot2_clim_df, "actvp_mean", cum_range = 2:14, cum_func = "sum")
plot2_clim_df <- calculate_cum_vals(plot2_clim_df, "light_mean", cum_range = 2:14, cum_func = "sum")

plot3_clim_df <- calculate_cum_vals(plot3_clim_df, "temp_mean", cum_range = 2:14, cum_func = "sum")
plot3_clim_df <- calculate_cum_vals(plot3_clim_df, "relhm_mean", cum_range = 2:14, cum_func = "mean")
plot3_clim_df <- calculate_cum_vals(plot3_clim_df, "satdf_mean", cum_range = 2:14, cum_func = "sum")
plot3_clim_df <- calculate_cum_vals(plot3_clim_df, "actvp_mean", cum_range = 2:14, cum_func = "sum")
plot3_clim_df <- calculate_cum_vals(plot3_clim_df, "light_mean", cum_range = 2:14, cum_func = "sum")

plot4_clim_df <- calculate_cum_vals(plot4_clim_df, "temp_mean", cum_range = 2:14, cum_func = "sum")
plot4_clim_df <- calculate_cum_vals(plot4_clim_df, "relhm_mean", cum_range = 2:14, cum_func = "mean")
plot4_clim_df <- calculate_cum_vals(plot4_clim_df, "satdf_mean", cum_range = 2:14, cum_func = "sum")
plot4_clim_df <- calculate_cum_vals(plot4_clim_df, "actvp_mean", cum_range = 2:14, cum_func = "sum")
plot4_clim_df <- calculate_cum_vals(plot4_clim_df, "light_mean", cum_range = 2:14, cum_func = "sum")

plot5_clim_df <- calculate_cum_vals(plot5_clim_df, "temp_mean", cum_range = 2:14, cum_func = "sum")
plot5_clim_df <- calculate_cum_vals(plot5_clim_df, "relhm_mean", cum_range = 2:14, cum_func = "mean")
plot5_clim_df <- calculate_cum_vals(plot5_clim_df, "satdf_mean", cum_range = 2:14, cum_func = "sum")
plot5_clim_df <- calculate_cum_vals(plot5_clim_df, "actvp_mean", cum_range = 2:14, cum_func = "sum")
plot5_clim_df <- calculate_cum_vals(plot5_clim_df, "light_mean", cum_range = 2:14, cum_func = "sum")

# Rice dataframe for each plot
plot1_rice_df <- adjust_rice_df_format(plot_number = 1)
plot2_rice_df <- adjust_rice_df_format(plot_number = 2)
plot3_rice_df <- adjust_rice_df_format(plot_number = 3)
plot4_rice_df <- adjust_rice_df_format(plot_number = 4)
plot5_rice_df <- adjust_rice_df_format(plot_number = 5)

# eDNA dataframe for each plot
edna_plot1 <- subset_samples(edna_all0, plot == 1)
edna_plot2 <- subset_samples(edna_all0, plot == 2)
edna_plot3 <- subset_samples(edna_all0, plot == 3)
edna_plot4 <- subset_samples(edna_all0, plot == 4)
edna_plot5 <- subset_samples(edna_all0, plot == 5)

plot1_edna_df <- adjust_edna_df_format(data.frame(sample_data(edna_plot1)), data.frame(otu_table(edna_plot1)), plot_number = 1)
plot2_edna_df <- adjust_edna_df_format(data.frame(sample_data(edna_plot2)), data.frame(otu_table(edna_plot2)), plot_number = 2)
plot3_edna_df <- adjust_edna_df_format(data.frame(sample_data(edna_plot3)), data.frame(otu_table(edna_plot3)), plot_number = 3)
plot4_edna_df <- adjust_edna_df_format(data.frame(sample_data(edna_plot4)), data.frame(otu_table(edna_plot4)), plot_number = 4)
plot5_edna_df <- adjust_edna_df_format(data.frame(sample_data(edna_plot5)), data.frame(otu_table(edna_plot5)), plot_number = 5)

# Check dimensions of dataframes
all(plot1_edna_df$date == plot1_clim_df$date)
all(plot1_edna_df$date == plot1_rice_df$date)
all(plot2_edna_df$date == plot2_clim_df$date)
all(plot2_edna_df$date == plot2_rice_df$date)
all(plot3_edna_df$date == plot3_clim_df$date)
all(plot3_edna_df$date == plot3_rice_df$date)
all(plot4_edna_df$date == plot4_clim_df$date)
all(plot4_edna_df$date == plot4_rice_df$date)
all(plot5_edna_df$date == plot5_clim_df$date)
all(plot5_edna_df$date == plot5_rice_df$date)

# Combine all to make composite time series
rice_all <- rbind(plot1_rice_df, plot2_rice_df, plot3_rice_df, plot4_rice_df, plot5_rice_df)
clim_all <- rbind(plot1_clim_df, plot2_clim_df, plot3_clim_df, plot4_clim_df, plot5_clim_df)
edna_all <- rbind(plot1_edna_df, plot2_edna_df, plot3_edna_df, plot4_edna_df, plot5_edna_df)
edna_all$plot <- as.factor(edna_all$plot)
edna_tax <- tax_table(edna_all0)

# Find column and row numbers to separate plots and data
colnames(rice_all); (rice_meta_coln <- 1:6); (rice_var_coln <- 7:ncol(rice_all))
colnames(clim_all); (clim_meta_coln <- 1:2); (clim_var_coln <- 3:ncol(clim_all))
colnames(edna_all); (edna_meta_coln <- 1:40); (edna_var_coln <- 41:ncol(edna_all))

rice_lib <- cbind(matrix(which(rice_all$date == "2017-05-01")+24, ncol = 1), matrix(which(rice_all$date == "2017-09-30")-8, ncol = 1))
clim_lib <- cbind(matrix(which(clim_all$date == "2017-05-01")+24, ncol = 1), matrix(which(clim_all$date == "2017-09-30")-8, ncol = 1))
edna_lib <- cbind(matrix(which(edna_all$date == "2017-05-01")+24, ncol = 1), matrix(which(edna_all$date == "2017-09-30")-8, ncol = 1))
all(rice_lib==clim_lib)
all(rice_lib==edna_lib)


# Determine the optimal embedding dimension using rUIC
E_RANGE <- 0:14

## eDNA time series
edna_std <- data.frame(x = as.numeric(scale(edna_all[,edna_var_coln[1]])))
Eedna_tmp <- rUIC::simplex(edna_std, lib_var = "x", lib = edna_lib,
                           cond_var = NULL, E = E_RANGE, tau = 1, tp = 1, Enull = "adaptive", n_boot = 2000, seed = 1234)
Eedna <- with(Eedna_tmp, max(c(0, E[pval < 0.05])))
for(i in edna_var_coln[-1]){
  start_time <- proc.time()[3]
  edna_std <- data.frame(x = as.numeric(scale(edna_all[,i])))
  Eedna_tmp <- rUIC::simplex(edna_std, lib_var = "x", lib = edna_lib,
                             cond_var = NULL, E = E_RANGE, tau = 1, tp = 1, Enull = "adaptive", n_boot = 2000, seed = 1234)
  Eedna[(i-edna_var_coln[1]+1)] <- with(Eedna_tmp, max(c(0, E[pval < 0.05])))
  cat(sprintf("eDNA Tax %04d / 1197 determined: %s sec\n", i - 40, round(proc.time()[3] - start_time, digits = 2)))
}
names(Eedna) <- colnames(edna_all)[edna_var_coln]

## Climate time series
clim_std <- data.frame(x = as.numeric(scale(clim_all[,clim_var_coln[1]])))
Eclim_tmp <- rUIC::simplex(clim_std, lib_var = "x", lib = clim_lib,
                           cond_var = NULL, E = E_RANGE, tau = 1, tp = 1, Enull = "adaptive", n_boot = 2000, seed = 1234)
Eclim <- with(Eclim_tmp, max(c(0, E[pval < 0.05])))
for(i in clim_var_coln[-1]){
  start_time <- proc.time()[3]
  clim_std <- data.frame(x = as.numeric(scale(clim_all[,i])))
  Eclim_tmp <- rUIC::simplex(clim_std, lib_var = "x", lib = clim_lib,
                             cond_var = NULL, E = E_RANGE, tau = 1, tp = 1, Enull = "adaptive", n_boot = 2000, seed = 1234)
  Eclim[(i-clim_var_coln[1]+1)] <- with(Eclim_tmp, max(c(0, E[pval < 0.05])))
  cat(sprintf("Cliamte variables %04d / 76 determined: %s sec\n", i - 2, round(proc.time()[3] - start_time, digits = 2)))
}
names(Eclim) <- colnames(clim_all)[clim_var_coln]

## Rice time series
rice_std <- data.frame(x = as.numeric(scale(rice_all[,rice_var_coln[1]])))
Erice_tmp <- rUIC::simplex(rice_std, lib_var = "x", lib = rice_lib,
                           cond_var = NULL, E = E_RANGE, tau = 1, tp = 1, Enull = "adaptive", n_boot = 2000, seed = 1234)
Erice <- with(Erice_tmp, max(c(0, E[pval < 0.05])))
for(i in rice_var_coln[-1]){
  start_time <- proc.time()[3]
  rice_std <- data.frame(x = as.numeric(scale(rice_all[,i])))
  Erice_tmp <- rUIC::simplex(rice_std, lib_var = "x", lib = rice_lib,
                             cond_var = NULL, E = E_RANGE, tau = 1, tp = 1, Enull = "adaptive", n_boot = 2000, seed = 1234)
  Erice[(i-rice_var_coln[1]+1)] <- with(Erice_tmp, max(c(0, E[pval < 0.05])))
  cat(sprintf("Rice variables %04d / 12 determined: %s sec\n", i - 6, round(proc.time()[3] - start_time, digits = 2)))
}
names(Erice) <- colnames(rice_all)[rice_var_coln]

# Save best E
saveRDS(Eedna, sprintf("%s/BestE_eDNA_ts.obj", output_folder01))
saveRDS(Eclim, sprintf("%s/BestE_Clim_ts.obj", output_folder01))
saveRDS(Erice, sprintf("%s/BestE_Rice_ts.obj", output_folder01))

# Save workspace and ojcects
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder01, output_folder01))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder01, substr(Sys.time(), 1, 10)))
