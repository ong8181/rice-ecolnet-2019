####
#### PRESTO 2019, Rice manipulation experiment
#### Visualization of CER rice yield data
####

# Create output folder
output_folder02 <- "02_VisualizeYieldOut"
dir.create(output_folder02)

# Load library
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2021.1.12
library(lubridate); packageVersion("lubridate") # 1.7.9.2, 2021.1.12
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.1.12
library(ggrepel); packageVersion("ggrepel") # 0.9.0, 2021.1.12
library(ggsci); packageVersion("ggsci") # 2.9, 2021.1.12
#library(ggrepel); packageVersion("ggrepel")
theme_set(theme_cowplot())

# load data
d <- read.csv("data/data_rice_yield.csv")

# Visualize data (Basic pattern)
g1 <- ggplot(d, aes(x = treatment, y = max_rice_height, group = treatment, color = treatment)) +
  geom_boxplot(outlier.color = "white", outlier.shape = NA, width = 0.4) +
  geom_jitter(height = 0, width = 0.2) +
  scale_color_manual(values = c("gray40", "red3", "royalblue")) +
  ylab("Rice height (cm)") + xlab(NULL)

g2 <- ggplot(d, aes(x = treatment, y = spad, group = treatment, color = treatment)) +
  geom_boxplot(outlier.color = "white", outlier.shape = NA, width = 0.4) +
  geom_jitter(height = 0, width = 0.2) +
  scale_color_manual(values = c("gray40", "red3", "royalblue")) +
  ylab("SPAD") + xlab(NULL)

g3 <- ggplot(d, aes(x = treatment, y = total_dry_wt, group = treatment, color = treatment)) +
  geom_boxplot(outlier.color = "white", outlier.shape = NA, width = 0.4) +
  geom_jitter(height = 0, width = 0.2) +
  scale_color_manual(values = c("gray40", "red3", "royalblue")) +
  ylab("Total dry weight (g)") + xlab(NULL)

g4 <- ggplot(d, aes(x = treatment, y = n_head, group = treatment, color = treatment)) +
  geom_boxplot(outlier.color = "white", outlier.shape = NA, width = 0.4) +
  geom_jitter(height = 0, width = 0.2) +
  scale_color_manual(values = c("gray40", "red3", "royalblue")) +
  ylab("No. of head") + xlab(NULL)

g_all1 <- plot_grid(g1 + theme(legend.position = "none"),
                    g2 + theme(legend.position = "none"),
                    g3 + theme(legend.position = "none"),
                    g4 + theme(legend.position = "none"),
                    ncol = 2, labels = "auto")


# Visualization version 2
g5 <- ggplot(d, aes(x = n_head, y = total_dry_wt, group = treatment, color = treatment, shape = treatment, fill = treatment)) +
  geom_point(alpha = 0.6, size = 3) +
  scale_color_manual(values = c("gray20", "red3", "royalblue")) +
  scale_fill_manual(values = c("gray20", "red3", "royalblue")) +
  geom_smooth(method = "lm", se = FALSE) +
  xlab("No. of head") + ylab("Total dry weight (g)")

g6 <- ggplot(d, aes(x = treatment, y = total_dry_wt/n_head, group = treatment, color = treatment)) +
  geom_boxplot(outlier.color = "white", outlier.shape = NA, width = 0.4) +
  geom_jitter(height = 0, width = 0.2) + scale_color_manual(values = c("gray40", "red3", "royalblue")) +
  ylab("Dry weight per head (g)") + xlab(NULL)

g_all2 <- plot_grid(g5 + theme(legend.position = "none"),
                    g6 + theme(legend.position = "none"), ncol = 2, align = "hv")


# Visualization version 3 ####
g7 <- ggplot(d, aes(y = prop_fert_ster_wt, x = treatment, group = treatment, color = treatment)) +
  geom_boxplot(outlier.color = "white", outlier.shape = NA, width = 0.4) +
  geom_jitter(height = 0, width = 0.2) + scale_color_manual(values = c("gray40", "red3", "royalblue")) +
  ylab("Grain per stem wt") + xlab(NULL)

g8 <- ggplot(d, aes(y = total_head_wt/n_head, x = treatment, group = treatment, color = treatment)) +
  geom_boxplot(outlier.color = "white", outlier.shape = NA, width = 0.4) +
  geom_jitter(height = 0, width = 0.2) + scale_color_manual(values = c("gray40", "red3", "royalblue")) +
  ylab("Head weight (g/head)") + xlab(NULL)

g9 <- ggplot(d, aes(y = fert_grain_wt/stem_only_wt, x = treatment, group = treatment, color = treatment)) +
  geom_boxplot(outlier.color = "white", outlier.shape = NA, width = 0.4) +
  geom_jitter(height = 0, width = 0.2) + scale_color_manual(values = c("gray40", "red3", "royalblue")) +
  ylab("Grain weight (g/head)") + xlab(NULL)

g_all3 <- plot_grid(g7 + theme(legend.position = "none"),
                    g8 + theme(legend.position = "none"),
                    g9 + theme(legend.position = "none"),
                    ncol = 2, align = "hv")

# Save figures
ggsave(sprintf("%s/RiceYield_Fig1.pdf", output_folder02), g_all1, width = 8, height = 8)
ggsave(sprintf("%s/RiceYield_Fig2.pdf", output_folder02), g_all2, width = 8, height = 5)
ggsave(sprintf("%s/RiceYield_Fig3.pdf", output_folder02), g_all3, width = 8, height = 8)
