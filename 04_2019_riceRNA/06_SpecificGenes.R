####
#### CER rice 2019
#### Rice RNA expression: Visualize specific genes
####

# Load workspace
load("05_DESeq2LocationOut/05_DESeq2LocationOut.RData")

# create directory
output_folder06 <- "06_SpecificGenesOut"
dir.create(output_folder06)

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)

# Load library
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2019.12.9
library(lubridate); packageVersion("lubridate") # 1.7.9.2, 2021.1.18
library(phyloseq); packageVersion("phyloseq") # 1.34.0, 2021.1.18
library(cowplot); packageVersion("cowplot") # 1.1.0, 2021.1.18
library(ggsci); packageVersion("ggsci") # 2.9, 2019.12.9
theme_set(theme_cowplot())
## DEG packages
library(DESeq2); packageVersion("DESeq2") # 1.30.0, 2021.1.18
library(ggpubr); packageVersion("ggpubr") # 0.4.0, 2021.1.18

# Revise the treatment names (2023.2.7)
treatment_ori <- c("CT", "PN", "RM")
treatment_rev <- c("CT", "GN", "CK")
deseq_rna_s1$treatment <- factor(deseq_rna_s1$treatment, levels = treatment_ori, labels = treatment_rev)
deseq_rna_s2$treatment <- factor(deseq_rna_s2$treatment, levels = treatment_ori, labels = treatment_rev)
deseq_rna_s3$treatment <- factor(deseq_rna_s3$treatment, levels = treatment_ori, labels = treatment_rev)
deseq_rna_s4$treatment <- factor(deseq_rna_s4$treatment, levels = treatment_ori, labels = treatment_rev)


# <----------------------------------------------------> #
# Check significant DEGs
# <----------------------------------------------------> #
length(which(deseq_res_s1_1$padj < 0.05))
length(which(deseq_res_s1_2$padj < 0.05))
length(which(deseq_res_s1l1_1$padj < 0.05))
length(which(deseq_res_s1l1_2$padj < 0.05))
length(which(deseq_res_s1l2_1$padj < 0.05))
length(which(deseq_res_s1l2_2$padj < 0.05))
length(which(deseq_res_s1l3_1$padj < 0.05))
length(which(deseq_res_s1l3_2$padj < 0.05))

length(which(deseq_res_s2_1$padj < 0.05))
length(which(deseq_res_s2_2$padj < 0.05))
length(which(deseq_res_s2l1_1$padj < 0.05))
length(which(deseq_res_s2l1_2$padj < 0.05))
length(which(deseq_res_s2l2_1$padj < 0.05))
length(which(deseq_res_s2l2_2$padj < 0.05))
length(which(deseq_res_s2l3_1$padj < 0.05))
length(which(deseq_res_s2l3_2$padj < 0.05))

length(which(deseq_res_s3_1$padj < 0.05))
length(which(deseq_res_s3_2$padj < 0.05))
length(which(deseq_res_s3l1_1$padj < 0.05))
length(which(deseq_res_s3l1_2$padj < 0.05))
length(which(deseq_res_s3l2_1$padj < 0.05))
length(which(deseq_res_s3l2_2$padj < 0.05))
length(which(deseq_res_s3l3_1$padj < 0.05))
length(which(deseq_res_s3l3_2$padj < 0.05))

length(which(deseq_res_s4_1$padj < 0.05))
length(which(deseq_res_s4_2$padj < 0.05))
length(which(deseq_res_s4l1_1$padj < 0.05))
length(which(deseq_res_s4l1_2$padj < 0.05))
length(which(deseq_res_s4l2_1$padj < 0.05))
length(which(deseq_res_s4l2_2$padj < 0.05))
length(which(deseq_res_s4l3_1$padj < 0.05))
length(which(deseq_res_s4l3_2$padj < 0.05))

# Generate subset phyloseq objects
dir.create(sprintf("%s/S1_DEG_FigEachGene", output_folder06))
dir.create(sprintf("%s/S2_DEG_FigEachGene", output_folder06))
dir.create(sprintf("%s/S3_DEG_FigEachGene", output_folder06))
dir.create(sprintf("%s/S4_DEG_FigEachGene", output_folder06))

## Sampling event 1
tmp_gene1 <- rownames(as.data.frame(deseq_res_s1_1[which(deseq_res_s1_1$padj < 0.05),]))
tmp_gene2 <- rownames(as.data.frame(deseq_res_s1_2[which(deseq_res_s1_2$padj < 0.05),]))
tmp_gene3 <- rownames(as.data.frame(deseq_res_s1l1_1[which(deseq_res_s1l1_1$padj < 0.05),]))
tmp_gene4 <- rownames(as.data.frame(deseq_res_s1l1_2[which(deseq_res_s1l1_2$padj < 0.05),]))
tmp_gene5 <- rownames(as.data.frame(deseq_res_s1l2_1[which(deseq_res_s1l2_1$padj < 0.05),]))
tmp_gene6 <- rownames(as.data.frame(deseq_res_s1l2_2[which(deseq_res_s1l2_2$padj < 0.05),]))
tmp_gene7 <- rownames(as.data.frame(deseq_res_s1l3_1[which(deseq_res_s1l3_1$padj < 0.05),]))
tmp_gene8 <- rownames(as.data.frame(deseq_res_s1l3_2[which(deseq_res_s1l3_2$padj < 0.05),]))
tmp_gene_s1 <- tmp_gene <- unique(c(tmp_gene1, tmp_gene2, tmp_gene3, tmp_gene4,
                                    tmp_gene5, tmp_gene6, tmp_gene7, tmp_gene8))
gene_s1_noloc <- unique(c(tmp_gene1, tmp_gene2))

for(i in tmp_gene){
        # Visualize results
        g0 <- data.frame(gene_exp = counts(deseq_rna_s1, normalize = T)[i,],
                         treatment = deseq_rna_s1$treatment,
                         plot = deseq_rna_s1$plot) %>%
                ggplot(aes(x = treatment, y = gene_exp + 0.5, color = plot, shape = plot, group = treatment)) +
                geom_boxplot(outlier.shape = NA, outlier.color = "white", width = 0.3, color = "gray80") +
                geom_jitter(width = 0.2, height = 0, size = 2) +
                scale_color_igv() + scale_shape_manual(values = 15:23) +
                scale_y_log10() +
                ylab(paste(i, "\n(DESeq2-normalized)")) + ggtitle("S1")
        ggsave(sprintf("%s/S1_DEG_FigEachGene/S1_%s.pdf", output_folder06, i),
               plot = g0, width = 5, height =4)
}
rm(tmp_gene); rm(g0)
rm(tmp_gene1); rm(tmp_gene2)
rm(tmp_gene3); rm(tmp_gene4)
rm(tmp_gene5); rm(tmp_gene6)
rm(tmp_gene7); rm(tmp_gene8)

## Sampling event 2
tmp_gene1 <- rownames(as.data.frame(deseq_res_s2_1[which(deseq_res_s2_1$padj < 0.05),]))
tmp_gene2 <- rownames(as.data.frame(deseq_res_s2_2[which(deseq_res_s2_2$padj < 0.05),]))
tmp_gene3 <- rownames(as.data.frame(deseq_res_s2l1_1[which(deseq_res_s2l1_1$padj < 0.05),]))
tmp_gene4 <- rownames(as.data.frame(deseq_res_s2l1_2[which(deseq_res_s2l1_2$padj < 0.05),]))
tmp_gene5 <- rownames(as.data.frame(deseq_res_s2l2_1[which(deseq_res_s2l2_1$padj < 0.05),]))
tmp_gene6 <- rownames(as.data.frame(deseq_res_s2l2_2[which(deseq_res_s2l2_2$padj < 0.05),]))
tmp_gene7 <- rownames(as.data.frame(deseq_res_s2l3_1[which(deseq_res_s2l3_1$padj < 0.05),]))
tmp_gene8 <- rownames(as.data.frame(deseq_res_s2l3_2[which(deseq_res_s2l3_2$padj < 0.05),]))
tmp_gene_s2 <- tmp_gene <- unique(c(tmp_gene1, tmp_gene2, tmp_gene3, tmp_gene4,
                                    tmp_gene5, tmp_gene6, tmp_gene7, tmp_gene8))
gene_s2_noloc <- unique(c(tmp_gene1, tmp_gene2))
gene_s2_loc <- data.frame(genes = unique(c(tmp_gene3, tmp_gene4, tmp_gene5, tmp_gene6, tmp_gene7, tmp_gene8)),
                          treatment = c(rep("GN", length(tmp_gene3)), rep("CK", length(tmp_gene4)),
                                        rep("GN", length(tmp_gene5)), rep("CK", length(tmp_gene6)),
                                        rep("GN", length(tmp_gene7)), rep("CK", length(tmp_gene8))),
                          loc =  c(rep("1", length(c(tmp_gene3, tmp_gene4))),
                                   rep("2", length(c(tmp_gene5, tmp_gene6))),
                                   rep("3", length(c(tmp_gene7, tmp_gene8)))))

for(i in tmp_gene){
        # Visualize results
        g0 <- data.frame(gene_exp = counts(deseq_rna_s2, normalize = T)[i,],
                         treatment = deseq_rna_s2$treatment,
                         plot = deseq_rna_s2$plot) %>%
                ggplot(aes(x = treatment, y = gene_exp + 0.5, color = plot, shape = plot, group = treatment)) +
                geom_boxplot(outlier.shape = NA, outlier.color = "white", width = 0.3, color = "gray80") +
                geom_jitter(width = 0.2, height = 0, size = 2) +
                scale_color_igv() + scale_shape_manual(values = 15:23) +
                scale_y_log10() +
                ylab(paste(i, "\n(DESeq2-normalized)")) + ggtitle("S2")
        ggsave(sprintf("%s/S2_DEG_FigEachGene/S2_%s.pdf", output_folder06, i),
               plot = g0, width = 5, height =4)
}
rm(tmp_gene); rm(g0)
rm(tmp_gene1); rm(tmp_gene2)
rm(tmp_gene3); rm(tmp_gene4)
rm(tmp_gene5); rm(tmp_gene6)
rm(tmp_gene7); rm(tmp_gene8)

## Sampling event 3
tmp_gene1 <- rownames(as.data.frame(deseq_res_s3_1[which(deseq_res_s3_1$padj < 0.05),]))
tmp_gene2 <- rownames(as.data.frame(deseq_res_s3_2[which(deseq_res_s3_2$padj < 0.05),]))
tmp_gene3 <- rownames(as.data.frame(deseq_res_s3l1_1[which(deseq_res_s3l1_1$padj < 0.05),]))
tmp_gene4 <- rownames(as.data.frame(deseq_res_s3l1_2[which(deseq_res_s3l1_2$padj < 0.05),]))
tmp_gene5 <- rownames(as.data.frame(deseq_res_s3l2_1[which(deseq_res_s3l2_1$padj < 0.05),]))
tmp_gene6 <- rownames(as.data.frame(deseq_res_s3l2_2[which(deseq_res_s3l2_2$padj < 0.05),]))
tmp_gene7 <- rownames(as.data.frame(deseq_res_s3l3_1[which(deseq_res_s3l3_1$padj < 0.05),]))
tmp_gene8 <- rownames(as.data.frame(deseq_res_s3l3_2[which(deseq_res_s3l3_2$padj < 0.05),]))
tmp_gene_s3 <- tmp_gene <- unique(c(tmp_gene1, tmp_gene2, tmp_gene3, tmp_gene4,
                                    tmp_gene5, tmp_gene6, tmp_gene7, tmp_gene8))
gene_s3_noloc <- unique(c(tmp_gene1, tmp_gene2))

for(i in tmp_gene){
        # Visualize results
        g0 <- data.frame(gene_exp = counts(deseq_rna_s3, normalize = T)[i,],
                         treatment = deseq_rna_s3$treatment,
                         plot = deseq_rna_s3$plot) %>%
                ggplot(aes(x = treatment, y = gene_exp + 0.5, color = plot, group = treatment, shape = plot)) +
                geom_boxplot(outlier.shape = NA, outlier.color = "white", width = 0.3, color = "gray80") +
                geom_jitter(width = 0.2, height = 0, size = 2) +
                scale_color_igv() + scale_shape_manual(values = 15:23) +
                scale_y_log10() +
                ylab(paste(i, "\n(DESeq2-normalized)")) + ggtitle("S3")
        ggsave(sprintf("%s/S3_DEG_FigEachGene/S3_%s.pdf", output_folder06, i),
               plot = g0, width = 5, height =4)
}
rm(tmp_gene); rm(g0)
rm(tmp_gene1); rm(tmp_gene2)
rm(tmp_gene3); rm(tmp_gene4)
rm(tmp_gene5); rm(tmp_gene6)
rm(tmp_gene7); rm(tmp_gene8)

## Sampling event 4
tmp_gene1 <- rownames(as.data.frame(deseq_res_s4_1[which(deseq_res_s4_1$padj < 0.05),]))
tmp_gene2 <- rownames(as.data.frame(deseq_res_s4_2[which(deseq_res_s4_2$padj < 0.05),]))
tmp_gene3 <- rownames(as.data.frame(deseq_res_s4l1_1[which(deseq_res_s4l1_1$padj < 0.05),]))
tmp_gene4 <- rownames(as.data.frame(deseq_res_s4l1_2[which(deseq_res_s4l1_2$padj < 0.05),]))
tmp_gene5 <- rownames(as.data.frame(deseq_res_s4l2_1[which(deseq_res_s4l2_1$padj < 0.05),]))
tmp_gene6 <- rownames(as.data.frame(deseq_res_s4l2_2[which(deseq_res_s4l2_2$padj < 0.05),]))
tmp_gene7 <- rownames(as.data.frame(deseq_res_s4l3_1[which(deseq_res_s4l3_1$padj < 0.05),]))
tmp_gene8 <- rownames(as.data.frame(deseq_res_s4l3_2[which(deseq_res_s4l3_2$padj < 0.05),]))
tmp_gene_s4 <- tmp_gene <- unique(c(tmp_gene1, tmp_gene2, tmp_gene3, tmp_gene4,
                                    tmp_gene5, tmp_gene6, tmp_gene7, tmp_gene8))
gene_s4_noloc <- unique(c(tmp_gene1, tmp_gene2))

for(i in tmp_gene){
        # Visualize results
        g0 <- data.frame(gene_exp = counts(deseq_rna_s4, normalize = T)[i,],
                         treatment = deseq_rna_s4$treatment,
                         plot = deseq_rna_s4$plot) %>%
                ggplot(aes(x = treatment, y = gene_exp + 0.5, color = plot, group = treatment, shape = plot)) +
                geom_boxplot(outlier.shape = NA, outlier.color = "white", width = 0.3, color = "gray80") +
                geom_jitter(width = 0.2, height = 0, size = 2) +
                scale_color_igv() + scale_shape_manual(values = 15:23) +
                scale_y_log10() +
                ylab(paste(i, "\n(DESeq2-normalized)")) + ggtitle("S4")
        ggsave(sprintf("%s/S4_DEG_FigEachGene/S4_%s.pdf", output_folder06, i),
               plot = g0, width = 5, height =4)
}
rm(tmp_gene); rm(g0)
rm(tmp_gene1); rm(tmp_gene2)
rm(tmp_gene3); rm(tmp_gene4)
rm(tmp_gene5); rm(tmp_gene6)
rm(tmp_gene7); rm(tmp_gene8)

# Check potential DEG names for overall analysis & location-specific analysis
tmp_gene_s1
tmp_gene_s2
tmp_gene_s3
tmp_gene_s4

#### save potential DEGs
s2_df <- cbind(data.frame(sample_name = deseq_rna_s2$Sample_Name2,
                          treatment = deseq_rna_s2$treatment,
                          plot = deseq_rna_s2$plot,
                          ind_id = deseq_rna_s2$individual_id),
               t(data.frame(counts(deseq_rna_s2, normalize = T)[gene_s2_noloc,])))
s4_df <- cbind(data.frame(sample_name = deseq_rna_s4$Sample_Name2,
                          treatment = deseq_rna_s4$treatment,
                          plot = deseq_rna_s4$plot,
                          ind_id = deseq_rna_s4$individual_id),
               t(data.frame(counts(deseq_rna_s4, normalize = T)[gene_s4_noloc,])))

# Save DEG list
#write.csv(s1_df, sprintf("%s/S1_DEGList.csv", output_folder06), row.names = F)
write.csv(s2_df, sprintf("%s/S2_DEGList.csv", output_folder06), row.names = F)
#write.csv(s3_df, sprintf("%s/S3_DEGList.csv", output_folder06), row.names = F)
write.csv(s4_df, sprintf("%s/S4_DEGList.csv", output_folder06), row.names = F)

# Add adj-p to the annotation
annot_s2 <- des[gene_s2_noloc,]
annot_s2$padj <- c(deseq_res_s2_1[which(deseq_res_s2_1$padj < 0.05),"padj"],
                   deseq_res_s2_2[which(deseq_res_s2_2$padj < 0.05),"padj"])
annot_s2$log2FoldChange <- c(deseq_res_s2_1[which(deseq_res_s2_1$padj < 0.05),"log2FoldChange"],
                             deseq_res_s2_2[which(deseq_res_s2_2$padj < 0.05),"log2FoldChange"])
annot_s4 <- des[gene_s4_noloc,]
annot_s4$padj <- c(deseq_res_s4_1[which(deseq_res_s4_1$padj < 0.05),"padj"],
                   deseq_res_s4_2[which(deseq_res_s4_2$padj < 0.05),"padj"])
annot_s4$log2FoldChange <- c(deseq_res_s4_1[which(deseq_res_s4_1$padj < 0.05),"log2FoldChange"],
                             deseq_res_s4_2[which(deseq_res_s4_2$padj < 0.05),"log2FoldChange"])
annot_s2_loc <- cbind(gene_s2_loc, des[gene_s2_loc$genes,])
annot_s2_loc$padj <- c(deseq_res_s2l1_2[which(deseq_res_s2l1_2$padj < 0.05),"padj"],
                       deseq_res_s2l3_1[which(deseq_res_s2l3_1$padj < 0.05),"padj"],
                       deseq_res_s2l3_2[which(deseq_res_s2l3_2$padj < 0.05),"padj"])
annot_s2_loc$log2FoldChange <- c(deseq_res_s2l1_2[which(deseq_res_s2l1_2$padj < 0.05),"log2FoldChange"],
                                 deseq_res_s2l3_1[which(deseq_res_s2l3_1$padj < 0.05),"log2FoldChange"],
                                 deseq_res_s2l3_2[which(deseq_res_s2l3_2$padj < 0.05),"log2FoldChange"])

# Save annotation
write.csv(annot_s2, sprintf("%s/S2_DEGannot.csv", output_folder06), row.names = F)
write.csv(annot_s4, sprintf("%s/S4_DEGannot.csv", output_folder06), row.names = F)
write.csv(annot_s2_loc, sprintf("%s/S2_DEGannot_Loc.csv", output_folder06), row.names = F)

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder06, substr(Sys.time(), 1, 10)))
