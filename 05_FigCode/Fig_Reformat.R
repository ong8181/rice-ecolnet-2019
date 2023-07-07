####
#### CER2019 Rice manipulation study
#### FigCode: Format figures
####

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2022.4.19
library(lubridate); packageVersion("lubridate") #1.8.0, 2022.4.19
library(cowplot); packageVersion("cowplot") # 1.1.1, 2022.4.19
library(ggimage); packageVersion("ggimage") # 0.3.0, 2022.4.19
library(magick); packageVersion("magick") # 2.7.3, 2022.4.19
library(ggsci); packageVersion("ggsci") # 2.9, 2018.8.10
library(ggpubr); packageVersion("ggpubr") # 0.4.0, 2023.2.9
library(scales); packageVersion("scales") # 1.1.1, 2021.3.24
library(patchwork); packageVersion("patchwork") # 1.1.1, 2022.4.19
theme_set(theme_cowplot())

# For statistical tests
library(lme4); packageVersion("lme4") # 1.1.31, 2023.01.11

# Specify figure folder
setwd("05_FigCode/")
fig_output <- "00_FigRaw"
formatfig_output <- "00_ReformatFigs"
dir.create(formatfig_output)


# <----------------------------------------------------> #
#               Load figure objects
# <----------------------------------------------------> #
# Load images
## Define a function to import image
import_image <- function (image_url) {
  img_imported <- image_url %>% magick::image_read()
  return(cowplot::ggdraw() + cowplot::draw_image(img_imported))
}
## Rice plot 2017
Fig_2017_RicePlot <- import_image("00_FieldImage_2017/PlotImages_v2.png")
## CER field 2017
Fig_2017_Field <- import_image("00_FieldImage_2017/PlotImages_v3.png")
## Monitoring method
Fig_eDNAMonitoring <- import_image("00_FieldImage_2017/ExpDesign_v3.jpg")
## Manipulative experiment 2019
Fig_2019_RicePlot <- import_image("00_FieldImage_2019/RiceExp2019_Setup.png")
Fig_2019_Pot <- import_image("00_FieldImage_2019/RiceExp2019_Pot.png")
Fig_2019_Midge <- import_image("00_FieldImage_2019/RiceExp2019_Midge.png")
Fig_2019_Pythium1 <- import_image("00_FieldImage_2019/RiceExp2019_Pythium1.png")
Fig_2019_Pythium2 <- import_image("00_FieldImage_2019/RiceExp2019_Pythium2.png")
# <----------------------------------------------------> #
# eDNA x Rice in 2017
Fig_DNAts <- readRDS("00_FigRaw/Fig_ushio2022.obj")
Fig_DNAxRice_RiceTS <- readRDS(sprintf("%s/Fig_DNAxRice_RiceTS.obj", fig_output))
Fig_DNAxRice_TargetS <- readRDS(sprintf("%s/Fig_DNAxRice_TargetTS.obj", fig_output))
Fig_DNAxRice_UIC <- readRDS(sprintf("%s/Fig_DNAxRice_UIChist.obj", fig_output))
# <----------------------------------------------------> #
# Rice growth data in 2019
Fig_RiceGrowth2019 <- readRDS(sprintf("%s/Fig_RiceGrowth2019.obj", fig_output))
# <----------------------------------------------------> #
# eDNA in 2019
#Fig_DNA2019bar <- readRDS(sprintf("%s/Fig_DNA2019bar_1.obj", fig_output))
#Fig_DNA2019div <- readRDS(sprintf("%s/Fig_DNA2019div_1.obj", fig_output))
#Fig_DNA2019flt_general <- readRDS(sprintf("%s/Fig_DNA2019flt_general.obj", fig_output))
Fig_DNA2019flt_target <- readRDS(sprintf("%s/Fig_DNA2019flt_target.obj", fig_output))
Fig_DNA2019flt_jitter <- readRDS(sprintf("%s/Fig_DNA2019flt_target_jitter.obj", fig_output))
Fig_DNA2019flt_tSNE <- readRDS(sprintf("%s/Fig_DNA2019flt_tSNE.obj", fig_output))
# <----------------------------------------------------> #
# RNA expression in 2019
#Fig_RNAhistogram <- readRDS(sprintf("%s/Fig_RNAhistogram.obj", fig_output))
Fig_RNAmaplot <- readRDS(sprintf("%s/Fig_RNAmaplot.obj", fig_output))
Fig_RNAmaplot_PN <- readRDS(sprintf("%s/Fig_RNAmaplot_PN.obj", fig_output))
Fig_RNAmaplot_RM <- readRDS(sprintf("%s/Fig_RNAmaplot_RM.obj", fig_output))
# <----------------------------------------------------> #
# DEG data
Fig_DEGs <- readRDS(sprintf("%s/Fig_RNA2019_DEG.obj", fig_output))
# <----------------------------------------------------> #
## Import from Ushio (2022) Proceedings B
Fig_Network_img <- image_read("00_FigRaw/Fig_ushio2022_network.jpg")
Fig_NetworkLegend_img <- image_read("00_FigRaw/Fig_ushio2022_legend.jpg")

# Revision (eDNA 2017 & 2019)
Fig_eDNAbar <- readRDS(sprintf("%s/Fig_DNA_Barplot.obj", fig_output))
# Revision (Climate 2019)
Fig_clim2019 <- readRDS(sprintf("%s/Fig_DNA2019_clim.obj", fig_output))


# <----------------------------------------------------> #
#                 Reformat figs
# <----------------------------------------------------> #
# Figure 1: Experimental design and data in 2017
plot_legend <- get_legend(Fig_DNAts[[2]] + theme(legend.position = "top"))
fg01 <- Fig_DNAxRice_RiceTS[[1]] + ylab("Rice growth rate\n(cm/day)")
fg02 <- Fig_DNAxRice_RiceTS[[2]] + ylab(expression(paste("Temperature (", degree, "C)")))
fg03 <- Fig_DNAts[[1]] + ylab("DNA\n(copies/ml water)") + theme(legend.position = "right")
fg04 <- Fig_DNAts[[2]]$data %>%
  ggplot(aes(x = date, y = value, color = plot)) +
  geom_point(size = 0.5) +
  geom_line() + ylab("No of ASV") +
  scale_color_startrek()
Fig_Rice2017 <- Fig_2017_RicePlot / (fg01 + fg02 + fg03 + fg04 +
  plot_layout(ncol = 2, byrow = FALSE)) +
  plot_annotation(tag_levels = c("a"))

# Figure 2: UIC results
Fig_UICres <- (((Fig_DNAxRice_UIC[[5]]  + labs(tag = "a")) +
                 (Fig_DNAxRice_TargetS[[1]] + theme(legend.position = "none") + labs(tag = "c") + ylab("eDNA conc.\n(copies/ml water)")) +
                 (Fig_DNAxRice_TargetS[[7]] + theme(legend.position = "bottom") + labs(tag = "d") + ylab("eDNA conc.\n(copies/ml water)")) +
               plot_layout(ncol = 1, heights = c(1,1,1.1))) |
                 (Fig_DNAxRice_UIC[[7]] + theme_light() +
                    geom_point(size = 2) +
                    theme(axis.text.y = element_text(size = 6),
                          legend.text = element_text(size = 12)) +
                    labs(tag = "b")))


# Figure 3: Ecological communities in 2019
Fig_DNA2019flt_jitter1 <- Fig_DNA2019flt_jitter[[1]] + labs(tag = "b") + theme(legend.position = "none") + ylab("Total eDNA (copies/ml water)")
Fig_DNA2019flt_jitter2 <- Fig_DNA2019flt_jitter[[2]] + labs(tag = "c") + theme(legend.position = "none") + ylab("Total eDNA (copies/ml water)")
## Add statistical clarity
text1 <- data.frame(x = c(1.5, 1.5, 1.5), y = c(20000, 20000, 20000),
                    lab = c("italic(P) == 0.032", "italic(P) < 2.0 %*% 10^{-16}", "italic(P) == 0.906"),
                    treatment = factor(c("CT", "GN", "CK"), levels = c("CT", "GN", "CK")))
text2 <- data.frame(x = c(1.5, 1.5, 1.5), y = c(20000, 20000, 20000),
                    lab = c("italic(P) == 0.045", "italic(P) == 6.6 %*% 10^{-6}", "italic(P) == 0.182"),
                    treatment = factor(c("CT", "GN", "CK"), levels = c("CT", "GN", "CK")))
Fig_DNA2019flt_jitter1 <- Fig_DNA2019flt_jitter1 + geom_text(data = text1, aes(x = x, y = y, label = lab), color = "black", parse = T, size = 4)
Fig_DNA2019flt_jitter2 <- Fig_DNA2019flt_jitter2 + geom_text(data = text2, aes(x = x, y = y, label = lab), color = "black", parse = T, size = 4)
Fig_DNA2019flt_jitter1 <- Fig_DNA2019flt_jitter1 + theme(plot.title = element_text(size = 12))
Fig_DNA2019flt_jitter2 <- Fig_DNA2019flt_jitter2 + theme(plot.title = element_text(size = 12))
Fig_DNA2019_1 <-(Fig_DNA2019flt_jitter1) / (Fig_DNA2019flt_jitter2)
Fig_DNA2019_2 <- Fig_DNA2019flt_tSNE[[4]] + labs(tag = "d") + theme(legend.position = "right") + stat_ellipse(type = "t", linetype = 1, level = 0.95)
Fig_DNA2019 <- (Fig_DNA2019_1 | Fig_DNA2019_2) + plot_layout(ncol = 2, widths = c(1, 1.5))
## Statistical tests
## Treatment specific (Pythium)
## Pythium all
Fig_DNA2019flt_jitter[[1]]$data %>%
  glmer(total_dna_conc ~ before_after*treatment + (1|plot), data = ., family = Gamma(link = "log")) %>% summary()
## Pythium each
Fig_DNA2019flt_jitter[[1]]$data %>% subset(treatment == "CT") %>%  # p = 0.032 (+)
  glmer(total_dna_conc ~ before_after + (1|plot), data = ., family = Gamma(link = "log")) %>% summary()
Fig_DNA2019flt_jitter[[1]]$data %>% subset(treatment == "GN") %>%  # p < 2.0e-16 (+)
  glmer(total_dna_conc ~ before_after + (1|plot), data = ., family = Gamma(link = "log")) %>% summary()
Fig_DNA2019flt_jitter[[1]]$data %>% subset(treatment == "CK") %>%  # p = 0.9062
  glmer(total_dna_conc ~ before_after + (1|plot), data = ., family = Gamma(link = "log")) %>% summary()
## Treatment specific (Midge)
## Midge all
Fig_DNA2019flt_jitter[[2]]$data %>% 
  glmer(total_dna_conc ~ before_after*treatment + (1|plot), data = ., family = Gamma(link = "log")) %>% summary()
## Midge each
Fig_DNA2019flt_jitter[[2]]$data %>% subset(treatment == "CT") %>%  # p = 0.0452 (-)
  glmer(total_dna_conc ~ before_after + (1|plot), data = ., family = Gamma(link = "log")) %>% summary()
Fig_DNA2019flt_jitter[[2]]$data %>% subset(treatment == "GN") %>%  # p = 6.64e-06 (+)
  glmer(total_dna_conc ~ before_after + (1|plot), data = ., family = Gamma(link = "log")) %>% summary()
Fig_DNA2019flt_jitter[[2]]$data %>% subset(treatment == "CK") %>%  # p = 0.182 (-)
  glmer(total_dna_conc ~ before_after + (1|plot), data = ., family = Gamma(link = "log")) %>% summary()
## Assemble figures
Fig_clim2019$layers[[1]] <- geom_line(linewidth = 0.3)
Fig_clim2019$layers[[2]] <- geom_point(size = 0.3)
Fig_Rice2019 <- Fig_2019_RicePlot + labs(tag = "a") +
  inset_element(Fig_2019_Pot, 0.72, 0.15, 0.93, 0.40) +
  inset_element(Fig_clim2019 + theme_classic(), 0.02, 0.65, 0.32, 0.91)
Fig_RiceDNA2019 <- (Fig_Rice2019) / (Fig_DNA2019)


# Figure 4: Rice growth rate in 2019
## Add statistical test results
text3 <- data.frame(x = c(2, 2), y = c(6.5, 6.5), lab = c("N.S.", "N.S."),
                    before_after = factor(c("Before", "After"), levels = levels(Fig_RiceGrowth2019[[5]]$data$before_after)))
text4 <- data.frame(x = c(2), y = c(34), lab = c("N.S."), before_after = factor("Before", levels = c("Before", "After")))
text5 <- data.frame(x = c(1.5, 2), y = c(23.7, 24.7), lab = c("italic(P) == 0.016", "italic(P) == 0.021"),
                    before_after = factor("After", levels = c("Before", "After")))
lines1 <- data.frame(x = c(1.1, 1.1), xend = c(1.9, 2.9), y = c(23.5, 24.5), yend = c(23.5, 24.5),
                    before_after = factor("After", levels = c("Before", "After")))
Fig_RiceGrowth2019_0_1 <- Fig_RiceGrowth2019[[4]] + ylim(0,7) +
  geom_text(data = text3, aes(x = x, y = y, label = lab), color = "black", parse = T, size = 4)
Fig_RiceGrowth2019_0_2 <- Fig_RiceGrowth2019[[5]]$data %>% 
  ggplot(aes(x = treatment, y = value, facet = before_after, color = treatment)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, color = "gray10") +
  facet_wrap(~ before_after, scales = "free_y") + panel_border() +
  geom_jitter(alpha = 1, width = 0.1, height = 0) +
  geom_text(data = text4, aes(x = x, y = y, label = lab), color = "black", parse = T, size = 4) +
  geom_text(data = text5, aes(x = x, y = y, label = lab), color = "black", parse = T, size = 4) +
  geom_segment(data = lines1, aes(x = x, xend = xend, y = y, yend = yend), color = "black") +
  ylab("Cumulative growth after the treatment (cm)") + xlab("Treatment") +
  scale_color_manual(values = c("gray70", "red3", "royalblue")) +
  theme(legend.position = "none")
## Assemble figures
Fig_RiceGrowth <- plot_grid(Fig_RiceGrowth2019_0_1, Fig_RiceGrowth2019_0_2,
                            labels = c("a", "b"),
                            nrow = 1, align = "hv", axis = "lrbt")
## Statistical tests
## Growth rate per day (overall)
Fig_RiceGrowth2019[[4]]$data %>% lmerTest::lmer(value ~ before_after*treatment + (1|ind/plot), data = .) %>% summary
## Growth rate per day (Before/after specific)
Fig_RiceGrowth2019[[4]]$data %>% subset(before_after == "Before") %>%  # p = N.S.
  lmer(value ~ treatment + (1|ind/plot), data = .) %>% multcomp::glht(linfct = multcomp::mcp(treatment = "Tukey")) %>% summary
Fig_RiceGrowth2019[[4]]$data %>% subset(before_after == "After") %>%  # p = N.S.
  lmer(value ~ treatment + (1|ind/plot), data = .) %>% multcomp::glht(linfct = multcomp::mcp(treatment = "Tukey")) %>% summary
## Cumulative growth (overall)
Fig_RiceGrowth2019[[5]]$data %>% mutate(plot = str_sub(ind, end = 2)) %>%
  lmerTest::lmer(value ~ before_after*treatment + (1|plot), data = .) %>% summary
  ## Cumulative growth (Before/after specific)
Fig_RiceGrowth2019[[5]]$data %>% mutate(plot = str_sub(ind, end = 2)) %>% 
  subset(before_after == "Before") %>%  # p = N.S.
  lmer(value ~ treatment + (1|plot), data = .) %>% multcomp::glht(linfct = multcomp::mcp(treatment = "Tukey")) %>% summary
Fig_RiceGrowth2019[[5]]$data %>% mutate(plot = str_sub(ind, end = 2)) %>% 
  subset(before_after == "After") %>%  # PN-CT, p = 0.0156; RM-CT, p = 0.0208
  lmer(value ~ treatment + (1|plot), data = .) %>% multcomp::glht(linfct = multcomp::mcp(treatment = "Tukey")) %>% summary


# Figure 5: RNA expression in 2019
Fig_RNA2019 <- plot_grid(Fig_RNAmaplot[[3]] + theme(legend.position = "none"),
                         Fig_RNAmaplot[[4]] + theme(legend.position = "none"),
                         ncol = 2)
#Fig_DNA2019flt_general[[1]] # DNA copies/ml time-series
#Fig_DNA2019flt_general[[2]] # ASV time-series


# Figure 6: Potential DEGs
set.seed(8181)
deg1 <- Fig_DEGs[[1]]
deg2 <- Fig_DEGs[[2]]
deg3 <- Fig_DEGs[[3]]
deg4 <- Fig_DEGs[[4]]
deg5 <- Fig_DEGs[[5]]
deg6 <- Fig_DEGs[[6]]
Fig_DEGs <- (deg1 + deg4 + deg6 + deg2 + deg3 + deg5 + plot_layout(ncol = 3, byrow = T) + plot_annotation(tag_levels = "a"))
## Independent NB GLMM
### DEG1
deg1$data %>% glmer.nb(value ~ treatment + (1|plot), data = .) %>%
  multcomp::glht(linfct = multcomp::mcp(treatment = "Tukey")) %>% summary
# GN-CT < 1e-04, GN-CK < 1e-04
### DEG2
deg2$data %>% glmer.nb(value ~ treatment + (1|plot), data = .) %>%
  multcomp::glht(linfct = multcomp::mcp(treatment = "Tukey")) %>% summary
# GN-CT < 1e-04, GN-CK < 1e-04
### DEG3
deg3$data %>% glmer.nb(value ~ treatment + (1|plot), data = .) %>%
  multcomp::glht(linfct = multcomp::mcp(treatment = "Tukey")) %>% summary
# GN-CT < 1e-04, GN-CK < 5e-04
### DEG4
deg4$data %>% glmer.nb(value ~ treatment + (1|plot), data = .) %>%
  multcomp::glht(linfct = multcomp::mcp(treatment = "Tukey")) %>% summary
# GN-CT < 1e-04, GN-CK < 1e-02
### DEG5
deg5$data %>% glmer.nb(value ~ treatment + (1|plot), data = .) %>%
  multcomp::glht(linfct = multcomp::mcp(treatment = "Tukey")) %>% summary
# GN-CT < 1e-04, GN-CK < 1e-04
### DEG6
deg6$data %>% glmer.nb(value ~ treatment + (1|plot), data = .) %>%
  multcomp::glht(linfct = multcomp::mcp(treatment = "Tukey")) %>% summary
# GN-CT < 1e-04, GN-CK < 1e-04
### For all figures, PN is statistically clearly different from the other two treatments


# <----------------------------------------------------> #
#  Format supplementary figures
# <----------------------------------------------------> #
## Supplementary figures
### Ecological community network
Fig_Network2 <- ggdraw() + draw_image(Fig_Network_img)
Fig_NetworkLegend2 <- ggdraw() + draw_image(Fig_NetworkLegend_img)
Fig_EcolNetwork <- plot_grid(Fig_Network2, Fig_NetworkLegend2,
                             ncol = 1, rel_heights = c(1, 0.11))
### Field setting and eDNA analysis workflow
Fig_Workflow <- ((Fig_2017_Field + labs(tag = "a")) +
                 (Fig_EcolNetwork + labs(tag = "c")) +
                  plot_layout(width = c(1.5, 1)))/
                 (Fig_eDNAMonitoring + labs(tag = "b")) +
                  plot_layout(height = c(1, 1.8))
### Pythium and midge images
Fig_2019_Organism <- (Fig_2019_Pythium2 + Fig_2019_Pythium1) / (Fig_2019_Midge + plot_spacer()) +
  plot_annotation(tag_levels = "a")
### DNA community in 2019
Fig_DNACom2019 <- plot_grid(Fig_DNA2019flt_target[[2]] + panel_border(),
                            Fig_DNA2019flt_target[[4]] + panel_border(),
                            ncol = 2, labels = "auto")

# For Supplementary Figure for rice growth trajectory in 2019
Fig_RiceGrowth_SI <- plot_grid(Fig_RiceGrowth2019[[1]], #+ xlim(c(ymd("2019-06-18"), ymd("2019-07-13"))),
                               Fig_RiceGrowth2019[[2]], #+ xlim(c(ymd("2019-06-18"), ymd("2019-07-13"))),
                               Fig_RiceGrowth2019[[3]], #+ xlim(c(ymd("2019-06-18"), ymd("2019-07-13"))),
                               labels = c("a", "b", "c"),
                               ncol = 1, align = "hv")

### RNA expression in 2019
Fig_RNA2019_loc <- plot_grid(Fig_RNAmaplot_PN[[4]] + theme(legend.position = "none"),
                             Fig_RNAmaplot_PN[[5]] + theme(legend.position = "none"),
                             Fig_RNAmaplot_PN[[6]] + theme(legend.position = "none"),
                             Fig_RNAmaplot_RM[[4]] + theme(legend.position = "none"),
                             Fig_RNAmaplot_RM[[5]] + theme(legend.position = "none"),
                             Fig_RNAmaplot_RM[[6]] + theme(legend.position = "none"),
                             ncol = 3, rel_widths = c(1,1,1), labels = "auto")

### eDNA patterns
Fig_DNAbar_2017 <- plot_grid(Fig_eDNAbar[[3]] + panel_border() +
                               theme(axis.text.x = element_text(vjust = 0.5)),
                             Fig_eDNAbar[[1]] + panel_border() +
                               theme(axis.text.x = element_text(vjust = 0.5)),
                             ncol = 1, align = "hv", axis = "lrbt", labels = "auto")
Fig_DNAbar_2019 <- plot_grid(Fig_eDNAbar[[4]] + panel_border() +
                               geom_segment(aes(x = "2019-06-18", y = 3.4e+7, xend = "2019-07-12", yend = 3.4e+7), color = "royalblue", linewidth = 1) +
                               theme(axis.text.x = element_text(vjust = 0.5)),
                             Fig_eDNAbar[[2]] + panel_border() +
                               theme(axis.text.x = element_text(vjust = 0.5)),
                             ncol = 1, align = "hv", axis = "lrbt", labels = "auto")


# <----------------------------------------------------> #
#               Save main figures
# <----------------------------------------------------> #
## Main figures
### Figure 1
ggsave(sprintf("%s/Figure_01.pdf", formatfig_output),
       plot = Fig_Rice2017,
       width = 14, height = 14)
### Figure 2
ggsave(sprintf("%s/Figure_02.pdf", formatfig_output),
       plot = Fig_UICres,
       width = 14, height = 10)
### Figure 3
ggsave(sprintf("%s/Figure_03.pdf", formatfig_output),
       plot = Fig_RiceDNA2019,
       width = 12, height = 14)
### Figure 4
ggsave(sprintf("%s/Figure_04.pdf", formatfig_output),
       plot = Fig_RiceGrowth,
       width = 14, height = 6)
### Figure 5
ggsave(sprintf("%s/Figure_05.jpg", formatfig_output),
       plot = Fig_RNA2019,
       dpi = 300, width = 10, height = 6)
### Figure 6
ggsave(sprintf("%s/Figure_06.pdf", formatfig_output),
       plot = Fig_DEGs,
       width = 9, height = 7)


# <----------------------------------------------------> #
#               Save Supplementary figures
# <----------------------------------------------------> #
## Save supplementary figures
### Figure S01
ggsave(filename = sprintf("%s/Figure_S01.jpg", formatfig_output),
       plot = Fig_Workflow,
       dpi = 300, width = 14, height = 10)
### Figure S02
ggsave(sprintf("%s/Figure_S02.pdf", formatfig_output),
       plot = Fig_2019_Organism,
       width = 9, height = 5)
### Figure S03
# ggsave(sprintf("%s/Figure_S03.pdf", formatfig_output),
#        plot = Fig_DNACom2019,
#        width = 12, height = 8)
ggsave(sprintf("%s/Figure_S03.pdf", formatfig_output),
       plot = Fig_RiceGrowth_SI,
       width = 12, height = 10)
### Figure S04
ggsave(sprintf("%s/Figure_S04.jpg", formatfig_output),
       plot = Fig_RNA2019_loc,
       dpi = 300, width = 12, height = 8)

### Figure R01-02
ggsave(sprintf("%s/Figure_R01.pdf", formatfig_output),
       plot = Fig_DNAbar_2017,
       width = 15, height = 15)
ggsave(sprintf("%s/Figure_R02.pdf", formatfig_output),
       plot = Fig_DNAbar_2019,
       width = 15, height = 20)


# <----------------------------------------------------> #
#             Save session information
# <----------------------------------------------------> #
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/SessionInfo_x1_FormatFigs_%s.txt", substr(Sys.time(), 1, 10)))


