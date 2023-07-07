####
#### PRESTO 2019, Rice manipulation experiment
#### Visualization of CER rice yield data
####

# Load library
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2021.1.12
library(lubridate); packageVersion("lubridate") # 1.7.9.2, 2021.1.12
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.1.12
theme_set(theme_cowplot())

# Create output folder
fig_output <- "00_FigRaw"
#dir.create(fig_output)

# Load data
result_folder01 <- "../02_2019_Rice" 
d <- read.csv(sprintf("%s/data/data_rice_yield.csv", result_folder01))
# Revise the treatment names (2023.2.7)
treatment_ori <- c("CT", "PN", "RM")
treatment_rev <- c("CT", "GN", "CK")
d <- d %>% mutate(treatment = factor(treatment, levels = treatment_ori, labels = treatment_rev))


# <----------------------------------------------------> #
#                 Generate figures
# <----------------------------------------------------> #
# Visualize data (Basic pattern)
g1 <- ggplot(d, aes(x = treatment, y = max_rice_height, group = treatment, color = treatment)) +
  geom_boxplot(outlier.color = "white", outlier.shape = NA, width = 0.4) +
  geom_jitter(height = 0, width = 0.2) +
  scale_color_manual(values = c("gray10", "red3", "royalblue")) +
  ylab("Rice height (cm)") + xlab(NULL) +
  NULL

g2 <- ggplot(d, aes(x = treatment, y = spad, group = treatment, color = treatment)) +
  geom_boxplot(outlier.color = "white", outlier.shape = NA, width = 0.4) +
  geom_jitter(height = 0, width = 0.2) +
  scale_color_manual(values = c("gray10", "red3", "royalblue")) +
  ylab("SPAD") + xlab(NULL) +
  NULL

g3 <- ggplot(d, aes(x = treatment, y = total_dry_wt, group = treatment, color = treatment)) +
  geom_boxplot(outlier.color = "white", outlier.shape = NA, width = 0.4) +
  geom_jitter(height = 0, width = 0.2) +
  scale_color_manual(values = c("gray10", "red3", "royalblue")) +
  ylab("Total dry weight (g)") + xlab(NULL) +
  NULL

g4 <- ggplot(d, aes(x = treatment, y = n_head, group = treatment, color = treatment)) +
  geom_boxplot(outlier.color = "white", outlier.shape = NA, width = 0.4) +
  geom_jitter(height = 0, width = 0.2) +
  scale_color_manual(values = c("gray10", "red3", "royalblue")) +
  ylab("No. of head") + xlab(NULL) +
  NULL

# Visualization version 2
g5 <- ggplot(d, aes(x = n_head, y = total_dry_wt, group = treatment, color = treatment, shape = treatment, fill = treatment)) +
  geom_point(alpha = 0.6, size = 3) +
  scale_color_manual(values = c("gray10", "red3", "royalblue")) +
  scale_fill_manual(values = c("gray10", "red3", "royalblue")) +
  geom_smooth(method = "lm", se = FALSE) +
  xlab("No. of head") + ylab("Total dry weight (g)") +
  NULL

g6 <- ggplot(d, aes(x = treatment, y = total_dry_wt/n_head, group = treatment, color = treatment)) +
  geom_boxplot(outlier.color = "white", outlier.shape = NA, width = 0.4) +
  geom_jitter(height = 0, width = 0.2) + scale_color_manual(values = c("gray40", "red3", "royalblue")) +
  ylab("Dry weight per head (g)") + xlab(NULL) +
  NULL

# Visualization version 3 ####
g7 <- ggplot(d, aes(y = prop_fert_ster_wt, x = treatment, group = treatment, color = treatment)) +
  geom_boxplot(outlier.color = "white", outlier.shape = NA, width = 0.4) +
  geom_jitter(height = 0, width = 0.2) + scale_color_manual(values = c("gray40", "red3", "royalblue")) +
  ylab("Grain per stem wt") + xlab(NULL) +
  NULL

g8 <- ggplot(d, aes(y = total_head_wt/n_head, x = treatment, group = treatment, color = treatment)) +
  geom_boxplot(outlier.color = "white", outlier.shape = NA, width = 0.4) +
  geom_jitter(height = 0, width = 0.2) + scale_color_manual(values = c("gray40", "red3", "royalblue")) +
  ylab("Head weight (g/head)") + xlab(NULL) +
  NULL

g9 <- ggplot(d, aes(y = fert_grain_wt/stem_only_wt, x = treatment, group = treatment, color = treatment)) +
  geom_boxplot(outlier.color = "white", outlier.shape = NA, width = 0.4) +
  geom_jitter(height = 0, width = 0.2) + scale_color_manual(values = c("gray40", "red3", "royalblue")) +
  ylab("Grain weight (g/head)") + xlab(NULL) +
  NULL


# <----------------------------------------------------> #
#                Save rice plot figures
# <----------------------------------------------------> #
# Save figures
rice_yield2019_figs <- list(g1, g2, g3, g4, g5, g6, g7, g8, g9)
saveRDS(rice_yield2019_figs, sprintf("%s/Fig_RiceYield2019.obj", fig_output))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/SessionInfo_04_RiceYield2019_%s.txt", substr(Sys.time(), 1, 10)))


