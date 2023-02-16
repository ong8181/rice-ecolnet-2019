####
#### CERrice2017 All data analysis
#### Visualizing the dynamics of two target taxa
####

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2021.1.12
library(cowplot); packageVersion("cowplot") # 1.1.0, 2021.1.12
library(ggsci); packageVersion("ggsci") # 2.9, 2021.1.12
theme_set(theme_cowplot())

# Load workspace
result_folder01 <- "../01_2017_DNAxRice"
load(sprintf("%s/03_CompileUICresOut/03_CompileUICresOut.RData", result_folder01))


# Create output folder
fig_output <- "00_FigRaw"
#dir.create(fig_output)


# <----------------------------------------------------> #
#                  Generate figures
# <----------------------------------------------------> #
# Set target taxa
tax1 <- "Fungi_Taxa00402"
tax2 <- "Inv_Taxa00042"
tax3 <- "Inv_Taxa00057"
tax4 <- "Inv_Taxa00139"
tax5 <- "Inv_Taxa00145"
tax6 <- "Inv_Taxa00176"

# Visualize eDNA dynamics
edna_fig <- na.omit(edna_all)
g1 <- ggplot(edna_fig, aes(x = date, y = Fungi_Taxa00402, colour = plot, group = plot)) +
  geom_line() + geom_point(size = 0.5) + scale_color_startrek() +
  geom_hline(yintercept = 0, linetype = 2) + xlab(NULL) +
  ylab("eDNA conc. (copies/ml water)") +
  ggtitle(expression(paste("putative ", italic("Globisporangium")))) +
  NULL

g2 <- ggplot(edna_fig, aes(x = date, y = Inv_Taxa00042, colour = plot, group = plot)) +
  geom_line() + geom_point(size = 0.5) + scale_color_startrek() +
  geom_hline(yintercept = 0, linetype = 2) + xlab(NULL) +
  ylab("Dynamics of Chironomus kiiensis 1\n(eDNA copies/ml water)") + #scale_y_log10() +
  NULL

g3 <- ggplot(edna_fig, aes(x = date, y = Inv_Taxa00057, colour = plot, group = plot)) +
  geom_line() + geom_point(size = 0.5) + scale_color_startrek() +
  geom_hline(yintercept = 0, linetype = 2) + xlab(NULL) +
  ylab("Dynamics of Chironomus kiiensis 2\n(eDNA copies/ml water)") + #scale_y_log10() +
  NULL

g4 <- ggplot(edna_fig, aes(x = date, y = Inv_Taxa00139, colour = plot, group = plot)) +
  geom_line() + geom_point(size = 0.5) + scale_color_startrek() +
  geom_hline(yintercept = 0, linetype = 2) + xlab(NULL) +
  ylab("Dynamics of Chironomus kiiensis 3\n(eDNA copies/ml water)") + #scale_y_log10() +
  NULL

g5 <- ggplot(edna_fig, aes(x = date, y = Inv_Taxa00145, colour = plot, group = plot)) +
  geom_line() + geom_point(size = 0.5) + scale_color_startrek() +
  geom_hline(yintercept = 0, linetype = 2) + xlab(NULL) +
  ylab("Dynamics of Chironomus kiiensis 4\n(eDNA copies/ml water)") + #scale_y_log10() +
  NULL

g6 <- ggplot(edna_fig, aes(x = date, y = Inv_Taxa00176, colour = plot, group = plot)) +
  geom_line() + geom_point(size = 0.5) + scale_color_startrek() +
  geom_hline(yintercept = 0, linetype = 2) + xlab(NULL) +
  ylab("Dynamics of Chironomus kiiensis 5\n(eDNA copies/ml water)") + #scale_y_log10() +
  NULL

g7 <- ggplot(edna_fig, aes(x = date,
                           y = Inv_Taxa00042 + Inv_Taxa00057 + Inv_Taxa00139 + Inv_Taxa00145 + Inv_Taxa00176,
                           colour = plot, group = plot)) +
  geom_line() + geom_point(size = 0.5) + scale_color_startrek() +
  geom_hline(yintercept = 0, linetype = 2) + xlab(NULL) +
  ylab("eDNA conc. (copies/ml water)") +
  ggtitle(expression(italic("Chironomus kiiensis"))) +
  #theme(title = element_text(size = 10)) + #scale_y_log10() +
  NULL


# <----------------------------------------------------> #
#                     Save figures
# <----------------------------------------------------> #
all_figs <- list(g1, g2, g3, g4, g5, g6, g7)
saveRDS(all_figs, sprintf("%s/Fig_DNAxRice_TargetTS.obj", fig_output))


# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/SessionInfo_02_DNAxRiceTarget_%s.txt", substr(Sys.time(), 1, 10)))

