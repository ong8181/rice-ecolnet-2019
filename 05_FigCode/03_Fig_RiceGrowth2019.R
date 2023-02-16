####
#### PRESTO 2019, Rice manipulation experiment
#### Visualization of CER rice growth data
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
d_gr <- read_csv(sprintf("%s/data/data_rice_growth.csv", result_folder01))
d_hd <- read_csv(sprintf("%s/data/data_rice_heading.csv", result_folder01))
d_yd <- read_csv(sprintf("%s/data/data_rice_yield.csv", result_folder01))

# Revise the treatment names (2023.2.7)
treatment_ori <- c("CT", "PN", "RM")
treatment_rev <- c("CT", "GN", "CK")
d_gr <- d_gr %>% mutate(treatment = factor(treatment, levels = treatment_ori, labels = treatment_rev))
d_hd <- d_hd %>% mutate(treatment = factor(treatment, levels = treatment_ori, labels = treatment_rev))
d_yd <- d_yd %>% mutate(treatment = factor(treatment, levels = treatment_ori, labels = treatment_rev))


# <----------------------------------------------------> #
#         Compile data before generating figures
# <----------------------------------------------------> #
d_gr$plot <- as.factor(d_gr$plot)
d_gr$treatment <- as.factor(d_gr$treatment)
d_hd$plot <- as.factor(d_hd$plot)
d_hd$treatment <- as.factor(d_hd$treatment)

# Add treatment date ("Before" and "After")
## The first treatment was 6/24
## The third (last) treatment was 6/28
treatment_date <- ymd_hms("2019-06-30 12:00:00") # One day after the date when the treatment completed
cumulative_end_date <- ymd("2019-07-10")
d_gr$before_after <- d_hd$before_after <- NaN
d_gr$before_after[which(d_gr$date_time < treatment_date)] <- "Before"
d_gr$before_after[which(d_gr$date_time > treatment_date)] <- "After"
d_hd$before_after[which(d_hd$date_time < treatment_date)] <- "Before"
d_hd$before_after[which(d_hd$date_time > treatment_date)] <- "After"

d_gr$before_after <- factor(d_gr$before_after, levels = c("Before", "After"))
d_hd$before_after <- factor(d_hd$before_after, levels = c("Before", "After"))

# Make the table longer
d_gr_long1 <- d_gr %>% select(date, plot, treatment, height1, height2, height3, before_after) %>%
  pivot_longer(-c(date, plot, treatment, before_after), names_to = "ind_id")
d_gr_long2 <- d_gr %>% select(date, plot, treatment, gr1, gr2, gr3, before_after) %>%
  pivot_longer(-c(date, plot, treatment, before_after), names_to = "ind_id")
d_gr_long3 <- d_gr %>% select(date, plot, treatment, spad1, spad2, spad3, before_after) %>%
  pivot_longer(-c(date, plot, treatment, before_after), names_to = "ind_id")
## Add individual IDs
d_gr_long1$ind <- paste0("p", d_gr_long1$plot, "_ind", str_sub(d_gr_long1$ind_id, start = -1))
d_gr_long2$ind <- paste0("p", d_gr_long2$plot, "_ind", str_sub(d_gr_long2$ind_id, start = -1))
d_gr_long3$ind <- paste0("p", d_gr_long3$plot, "_ind", str_sub(d_gr_long3$ind_id, start = -1))
## Calculate plot-mean values
d_ht_mean <- d_gr_long1 %>% group_by(date, treatment) %>% summarize(plot_mean = mean(value))
d_gr_mean <- d_gr_long2 %>% group_by(date, treatment) %>% summarize(plot_mean = mean(value))
d_sp_mean <- d_gr_long3 %>% group_by(date, treatment) %>% summarize(plot_mean = mean(value))

# Adding up growth rate after treatment
d_cum_gr <- d_gr_long2[d_gr_long2$date <= cumulative_end_date,] %>% group_by(before_after, ind) %>%
  summarize(value = sum(value),
            treatment = unique(treatment))


# <----------------------------------------------------> #
#                    Generate figures
# <----------------------------------------------------> #
p1 <- ggplot(d_gr_long1, aes(x = date, y = value, color = treatment, group = treatment)) +
  geom_vline(xintercept = ymd("2019-06-24"), size = 1, color = "gray10", alpha = 0.5) +
  geom_vline(xintercept = ymd("2019-06-26"), size = 1, color = "gray10", alpha = 0.5) +
  geom_vline(xintercept = ymd("2019-06-28"), size = 1, color = "gray10", alpha = 0.5) +
  geom_point(alpha = 0.5, size = 0.8) + scale_color_manual(values = c("gray70", "red3", "royalblue")) +
  geom_line(data = d_ht_mean, aes(x = date, y = plot_mean, color = treatment, group = treatment), size = 1) +
  theme(axis.text.x = element_text(size = 12)) + xlab(NULL) + ylab("Rice height (cm)") +
  ggtitle("Rice height at CER 2019") + geom_hline(yintercept = 0, linetype = 2) +
  NULL

p2 <- ggplot(d_gr_long2, aes(x = date, y = value, color = treatment, group = treatment)) +
  geom_vline(xintercept = ymd("2019-06-24"), size = 1, color = "gray10", alpha = 0.5) +
  geom_vline(xintercept = ymd("2019-06-26"), size = 1, color = "gray10", alpha = 0.5) +
  geom_vline(xintercept = ymd("2019-06-28"), size = 1, color = "gray10", alpha = 0.5) +
  geom_point(alpha = 0.5, size = 0.8) + scale_color_manual(values = c("gray70", "red3", "royalblue")) +
  geom_line(data = d_gr_mean, aes(x = date, y = plot_mean, color = treatment, group = treatment), size = 1) +
  theme(axis.text.x = element_text(size = 12)) + xlab(NULL) + ylab("Rice growth rate (cm/day)") +
  ggtitle("Rice growth rate at CER 2019") + geom_hline(yintercept = 0, linetype = 2) + 
  NULL

p3 <- ggplot(d_gr_long3, aes(x = date, y = value, color = treatment, group = treatment)) +
  geom_vline(xintercept = ymd("2019-06-24"), size = 1, color = "gray10", alpha = 0.5) +
  geom_vline(xintercept = ymd("2019-06-26"), size = 1, color = "gray10", alpha = 0.5) +
  geom_vline(xintercept = ymd("2019-06-28"), size = 1, color = "gray10", alpha = 0.5) +
  geom_point(alpha = 0.5, size = 0.8) + scale_color_manual(values = c("gray70", "red3", "royalblue")) +
  geom_line(data = d_sp_mean, aes(x = date, y = plot_mean, color = treatment, group = treatment), size = 1) +
  theme(axis.text.x = element_text(size = 12)) + xlab(NULL) + ylab("SPAD") +
  ggtitle("Rice SPAD at CER 2019") + geom_hline(yintercept = 0, linetype = 2) + 
  NULL

p4 <- ggplot(d_gr_long2, aes(x = treatment, y = value, facet = before_after, color = treatment)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, color = "gray10") +
  facet_wrap(~ before_after) + panel_border() +
  geom_jitter(alpha = 0.5, width = 0.1, height = 0) +
  ylab("Growth rate (cm/day)") + xlab("Treatment") +
  scale_color_manual(values = c("gray70", "red3", "royalblue")) +
  NULL

p5 <- ggplot(d_cum_gr, aes(x = treatment, y = value, facet = before_after, color = treatment)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, color = "gray10") +
  facet_wrap(~ before_after) + panel_border() +
  geom_jitter(alpha = 1, width = 0.1, height = 0) +
  ylab("Cumulative growth\nafter the treatment (cm)") + xlab("Treatment") +
  scale_color_manual(values = c("gray70", "red3", "royalblue")) +
  NULL


# <----------------------------------------------------> #
#                Save rice plot figures
# <----------------------------------------------------> #
rice_gr2019_figs <- list(p1, p2, p3, p4, p5)
saveRDS(rice_gr2019_figs, sprintf("%s/Fig_RiceGrowth2019.obj", fig_output))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/SessionInfo_03_RiceGrowth2019_%s.txt", substr(Sys.time(), 1, 10)))

