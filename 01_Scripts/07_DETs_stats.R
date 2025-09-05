# Info -------------------------------------------------------------------------

# Temperature effect on gene expression
# Differential gene expression analysis
# Sebastes faciatus in common garden
# Using temperature as continuous variable
# And separating individuals from long vs short term temperature exposure
# 
# CL - Janvier 2023
# 

# Library ----------------------------------------------------------------------
rm(list = ls())
gc()

library(tidyverse)

library(vegan)
library(DESeq2)

library(eulerr)
library(ComplexHeatmap)
library(ggplotify)
library(ggpubr)
library(ggVennDiagram)
library(ggrepel)

library(topGO)
library(trinotateR)

library(here)

library(WGCNA)
library(RColorBrewer)

library(rrvgo)

# Create new dir for result -----------------------------------------------------

new_dir = "07_DETs_results_clean"
dir.create(file.path(here::here(), "02_Results", new_dir),
           showWarnings = FALSE)

# Library with overrepresented sequences ---------------------------------------

multiqc_R1 <- read.table(file.path(here::here(), "02_Results", "01_FastQC",
                                 "02_TrimGalore", "MultiQC_report",
                                 "fastqc_overrepresented_sequencesi_plot_R1.tsv"),
                         header = T)
names(multiqc_R1) <- c("ID", "top", "remain")

df <- multiqc_R1 %>% reshape2::melt() 

fig.R1 <- ggplot(df) +
  geom_bar(aes(x = ID, y=value, fill = variable),
           stat = "summary") +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  xlab("") + ylab("Percentage of Total Sequences (%)") +
  scale_fill_discrete(labels = c("Top over-represented sequences",
                                 "Sum of remaining over-represented sequences")) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.text.x=element_text(size = 8)) +
  ggtitle("FastQC: Overrepresented sequences for R1")
  
fig.R1

multiqc_R2 <- read.table(file.path(here::here(), "02_Results", "01_FastQC",
                                   "02_TrimGalore", "MultiQC_report",
                                   "fastqc_overrepresented_sequencesi_plot_R2.tsv"),
                         header = T)
names(multiqc_R2) <- c("ID", "top", "remain")

df <- multiqc_R2 %>% reshape2::melt() 

fig.R2 <- ggplot(df) +
  geom_bar(aes(x = ID, y=value, fill = variable),
           stat = "summary") +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  xlab("") + ylab("Percentage of Total Sequences (%)") +
  scale_fill_discrete(labels = c("Top over-represented sequences",
                                 "Sum of remaining over-represented sequences")) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.text.x=element_text(size = 8)) +
  ggtitle("FastQC: Overrepresented sequences for R2")

fig.R2

fig.fastqc <- ggarrange(fig.R1, fig.R2,
                       nrow = 2, labels = LETTERS, common.legend = TRUE,
                       legend = "bottom") +
  theme(plot.background = element_rect(fill = "white"))

fig.fastqc

ggsave(file.path(here::here(), "02_Results", new_dir, "FastQC_overrepresented_sequences.png"), 
       fig.fastqc,
       width = 7, height = 8, scale = 2)

# Data -------------------------------------------------------------------------

metaData.init <- read.csv("00_Data/00_Data_Info/metaData_all.csv")
metaData.init$ID <- paste(metaData.init$Seq_batch, metaData.init$ID_sample, sep = ".")

# remove the two libraries whith high amount of rRNA
metaData.init <- subset(metaData.init, !ID %in% c("1828.ARN_22_00004", "1828.ARN_22_00034"))

# Count table
counts.init <- read.table("00_Data/04_FeatureCounts/Count_transcript.txt", header = T, row.names = 1)[,-(1:5)]

colnames(counts.init) <- colnames(counts.init) %>% 
  str_replace(., "X.media.genobiwan.Extra_Storage.Projets.REDTANK_Transcripto.00_Data.02_Hisat2.",
              "") %>%
  str_replace(., "_hisat2.bam", "")

counts.init <- subset(counts.init, select = metaData.init$ID)

# Test for sequencing batch effect ---------------------------------------------

## Subset samples sequenced in both batches ----

ID.duplicated <- (metaData.init %>% filter(duplicated(.[["ID_sample"]])))$ID_sample

metaData.batches <- subset(metaData.init, ID_sample %in% ID.duplicated)

counts.batches <- subset(counts.init, select = metaData.batches$ID)
counts.batches <- round(counts.batches, 0)

metaData.batches <- metaData.batches[match(colnames(counts.batches), metaData.batches$ID),]
metaData.batches$Seq_batch <- as.factor(metaData.batches$Seq_batch)

# individual order
identical(colnames(counts.batches), metaData.batches$ID)

# Design
design <- as.formula(~ Seq_batch + ID_sample)

ddsObj.batch <- DESeqDataSetFromMatrix(
  countData = counts.batches,
  colData = metaData.batches,
  design = design
)

# Normalized data
count.vst.batch <- vst(ddsObj.batch, blind = T)

## Principal Component Analyis - PCA ----
pca.df.batch <-  plotPCA(count.vst.batch,
                         intgroup = c("Seq_batch", "ID_sample"),
                         ntop = nrow(ddsObj.batch), returnData = TRUE)

percentVar <- round(100 * attr(pca.df.batch, "percentVar"), 2)

colvec <- RColorBrewer::brewer.pal(n = 5, name = 'Set1')

batch.pca.plot <- ggplot(pca.df.batch, aes(PC1, PC2, 
                                           color = ID_sample, 
                                           shape = Seq_batch)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey") +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw() +
  labs(shape = "Sequencing batch") +
  labs(colour = "Samples") +
  scale_colour_manual(values = colvec) +
  scale_shape_manual(values = c(19, 17), labels = c("# 1812", "# 1828")) +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank()
  ) 

batch.pca.plot
ggsave(file.path(here::here(), "02_Results", new_dir, "PCA_batch_effect.png"), 
       batch.pca.plot,
       width = 3, height = 2, scale = 1.5)

## RDA ----

dist.batch <- vegdist(t(assay(count.vst.batch)), method = "euclidean")

rda.batch <- capscale(dist.batch ~ Seq_batch + Condition(ID_sample), metaData.batches)
anova.cca(rda.batch)
RsquareAdj(rda.batch)

rda.sample <- capscale(dist.batch ~ ID_sample + Condition(Seq_batch), metaData.batches)
anova.cca(rda.sample)
RsquareAdj(rda.sample)


res.batch.df <- data.frame(Effect = c("Batch", "Sample"),
                           rbind(data.frame(anova.cca(rda.batch))[1,],
                                 data.frame(anova.cca(rda.sample))[1,]),
                           R2 = c(RsquareAdj(rda.batch)$r.squared,
                                  RsquareAdj(rda.sample)$r.squared),
                           adj.R2 = c(RsquareAdj(rda.batch)$adj.r.squared,
                                      RsquareAdj(rda.sample)$adj.r.squared),
                           Residual = c(data.frame(anova.cca(rda.batch))[2,2],
                                        data.frame(anova.cca(rda.sample))[2,2]))

# write.csv(res.batch.df,
#           file = file.path(here::here(), "02_Results", new_dir, "Results_RDA_Batch_effect.csv"),
#           row.names = F)

# DETs across temperature -----------------------------------------------------

## subsample samples ----
# individual order
identical(colnames(counts.batches), metaData.batches$ID)

metaData <- subset(metaData.init, ! ID_sample %in% metaData.batches$ID_sample)
metaData <- rbind(metaData, subset(metaData.batches, Seq_batch == 1812))

metaData$acclim <- as.factor(metaData$Temperature_acclim)
metaData$sampling <- as.factor(metaData$Temperature_sampling)
metaData$Delta_short <- factor(metaData$Temperature_sampling - metaData$Temperature_acclim)
metaData$Delta_long <- factor(metaData$Temperature_acclim - 5.0)

counts.all <- round(subset(counts.init, select = metaData$ID),0)

meta.all <- metaData
counts.all <- counts.all

# write.csv(meta.all, 
#           file = file.path(here::here(), "00_Data", "metaData_subset_20nov2024.csv"))
# write.csv(counts.all, 
#           file = file.path(here::here(), "00_Data", "countsRNA_subset_20nov2024.csv"))


# Design
design <- as.formula(~ 1)
ddsObj <- DESeqDataSetFromMatrix(
  countData = counts.all,
  colData = meta.all,
  design = design
)

counts.vst.all <- vst(ddsObj, blind = T)

rm("metaData")

# For short and long-term temperature exposure

meta.short <- subset(meta.all, Temperature_acclim %in% c(5, 7.5))
meta.long <- subset(meta.all, Delta_short == 0)


# DETs across temperature for long term acclimation-----------------------------

metaData <- meta.long
counts <- subset(counts.all, select = metaData$ID)

# individual order
identical(colnames(counts), metaData$ID)

# Design
design <- as.formula(~ 1)
ddsObj <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metaData,
  design = design
)

# Normalized data
count.vst <- subset(assay(counts.vst.all), select = metaData$ID)

# Distance 

count.dist <- vegdist(t(count.vst), method = "euclidean")

# RDA

# Delta_acclim effect
rda.delta <- capscale(count.dist ~ metaData$Delta_acclim + Condition(metaData$Temperature_acclim))
anova.cca(rda.delta)
RsquareAdj(rda.delta)

rda.long.c <- capscale(count.dist ~ metaData$Temperature_acclim)
anova.cca(rda.long.c)
RsquareAdj(rda.long.c)
# 
# rda.long.f <- capscale(count.dist ~ metaData$acclim)
# anova.cca(rda.long.f)
# RsquareAdj(rda.long.f)


res.long.df <- data.frame(Effect = c("Acclim_10vs3_months", "Acclimation"),
                           rbind(data.frame(anova.cca(rda.delta))[1,],
                                 data.frame(anova.cca(rda.long.c))[1,]),
                           R2 = c(RsquareAdj(rda.delta)$r.squared,
                                  RsquareAdj(rda.long.c)$r.squared),
                           adj.R2 = c(RsquareAdj(rda.delta)$adj.r.squared,
                                      RsquareAdj(rda.long.c)$adj.r.squared),
                           Residual = c(data.frame(anova.cca(rda.delta))[2,2],
                                        data.frame(anova.cca(rda.long.c))[2,2]))

# write.csv(res.long.df,
#           file = file.path(here::here(), "02_Results", new_dir, "Results_RDA_Long-term.csv"),
#           row.names = F)
# 


## full RDA  and plot

## with or withou sex effect

rda.1 <- capscale(count.dist ~ metaData$Temperature_acclim + metaData$Delta_acclim)
RsquareAdj(rda.1)
anova.cca(rda.1, by = "margin")

rda.2 <- capscale(count.dist ~ metaData$Temperature_acclim + metaData$Delta_acclim + substring(metaData$Sex_genet, 1, 1))
RsquareAdj(rda.2)
anova.cca(rda.2, by = "margin")

RsquareAdj(capscale(count.dist ~ metaData$acclim + Condition(metaData$Delta_acclim + substring(metaData$Sex_genet, 1, 1))))
RsquareAdj(capscale(count.dist ~ metaData$Delta_acclim + Condition(metaData$acclim + substring(metaData$Sex_genet, 1, 1))))
RsquareAdj(capscale(count.dist ~ substring(metaData$Sex_genet, 1, 1) + Condition(metaData$acclim + metaData$Delta_acclim)))

# PCA plot

## Principal Component Analyis - PCA ----
pca.long <- rda(t(count.vst))

pca.df.long <- data.frame(scores(pca.long, display = "sites"),
                          acclim = metaData$acclim)

percentVar <- round(summary(pca.long)$cont$importance[2,1:2] * 100, 2)

pca.long.plot <- ggplot(pca.df.long, aes(PC1, PC2, 
                                           color = acclim)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey") +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(name = "Long-term\ntemperature (°C)",
                     values = c("#0571B0", "#92C5DE", "#F4A582", "#CA0020")) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank()
  ) +
  theme(plot.margin = margin(.5, .5, .5, .5, "cm")) 

pca.long.plot

# RDA Ordination

rda.df.long <- data.frame(scores(rda.1, display = "sites"),
                          acclim = metaData$acclim)

percentVar <- round(summary(rda.1)$cont$importance[2,1:2] * 100, 2)

rda.long.plot <- ggplot(rda.df.long, aes(CAP1, CAP2, 
                                         color = acclim)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey") +
  geom_point(size = 4) +
  xlab(paste0("RDA1: ", percentVar[1], "% variance")) +
  ylab(paste0("RDA2: ", percentVar[2], "% variance")) +
  scale_color_manual(name = "Long-term\ntemperature (°C)",
                     values = c("#0571B0", "#92C5DE", "#F4A582", "#CA0020")) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank()
  ) +
  theme(plot.margin = margin(.5, .5, .5, .5, "cm")) 

rda.long.plot


## Result for all - wald test

ddsObj <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metaData,
  design = design
)

# Make 5C as reference temperature
ddsObj$acclim <- relevel(ddsObj$acclim, ref = "5")

design(ddsObj) <- formula(~ acclim)
dds.all<- DESeq(ddsObj)
resultsNames(dds.all)

padj.cutoff <- 0.05
log2FC <- 1

# Long-term exposure
res.a.down <- results(dds.all, name = "acclim_2.5_vs_5")
res.a.down.sig <- subset(res.a.down, padj < padj.cutoff & abs(log2FoldChange) >= log2FC)
summary(res.a.down.sig)

res.a.up <- results(dds.all, name = "acclim_7.5_vs_5")
res.a.up.sig <- subset(res.a.up, padj < padj.cutoff & abs(log2FoldChange) >= log2FC)
summary(res.a.up.sig)

res.a.up2 <- results(dds.all, name = "acclim_10_vs_5")
res.a.up2.sig <- subset(res.a.up2, padj < padj.cutoff & abs(log2FoldChange) >= log2FC)
summary(res.a.up2.sig)

DET.acclim <- rbind(res.a.down %>% data.frame %>% mutate(Compare = "acclim_2.5_vs_5"),
                    res.a.up %>% data.frame %>% mutate(Compare = "acclim_7.5_vs_5"),
                    res.a.up2 %>%  data.frame %>% mutate(Compare = "acclim_10_vs_5"))
# 
# write.csv(DET.acclim,
#           file.path(here(), "02_Results", new_dir,
#                     "DETs_long_term.csv"))

## Venn plot
s1 <- list("5.0°C to 2.5°C" = row.names(res.a.down.sig),
           "5.0°C to 7.5°C" = row.names(res.a.up.sig),
           "5.0°C to 10.0°C" = row.names(res.a.up2.sig)
)

venn.plot.long.euler <- plot(euler(s1, shape = "circle"),
                             quantities = list(T, cex = .75),
                             legend = list(T, cex = 0.75, side = "right"),
                             fill = c("#0571B0", "#F4A582", "#CA0020"),
                             alpha = 0.75, edges = F) %>%
  ggplotify::as.ggplot(.) +
  theme(plot.margin = margin(.2, .2, .2, .2, "cm")) 

venn.plot.long.euler  


## Group genes expression ---------------------
deg_list.a <- unique(c(rownames(res.a.down.sig),
                       rownames(res.a.up.sig),
                       rownames(res.a.up2.sig)))

##deg_list.a <- unique(rownames(res.a.cont.sig))
deg_vst.a <- count.vst[deg_list.a, ]

meta.temp <- metaData
rownames(meta.temp) <- paste(metaData$Seq_batch, metaData$ID_sample, sep = ".")
meta.temp$acclim <- as.factor(meta.temp$Temperature_acclim)
meta.temp$sampling <- as.factor(meta.temp$Temperature_sampling)

clusters_acclim <- DEGreport::degPatterns(deg_vst.a, 
                                          metadata = meta.temp, 
                                          time = "acclim", col = NULL,
                                          minc = 15
) 

df.longterm <- clusters_acclim[["normalized"]] %>% 
  group_by(cluster) %>% 
  mutate(sum_cluster = length(unique(genes))) %>% 
  ungroup() %>% arrange(desc(sum_cluster)) %>%
  mutate(sum_cluster2 = factor(sum_cluster, levels = unique(sum_cluster))) %>%
  mutate(cluster2 = as.integer(sum_cluster2)) %>%
  mutate(cluster_new = paste("Group L", formatC(cluster2, width=2, flag="0"), " (", sum_cluster, ")", sep = "")) %>%
  data.frame()

cluster.acclim.plot <- ggplot(df.longterm, aes(x = Temperature_acclim, y = value)) +
  geom_violin(aes(group = acclim), alpha = 0.5, fill = "grey") +
  geom_smooth(#aes(group = acclim),
    col = "red", 
    linetype="dashed", method = "gam", formula = y ~ poly(x, 2)) +
  scale_x_continuous(name = "Long-term temperature (°C)", breaks = c(2.5, 5, 7.5, 10)) +
  xlab("Long-term temperature (°C)") +
  ylab("Z-score of transcript abundance") +
  theme_bw() +
  facet_wrap(~ cluster_new, ncol = 5) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  theme(plot.margin = margin(.5, .5, .5, .5, "cm")) 

cluster.acclim.plot

### Export results long-term ---------------------------------------------------

figure.long.up <- ggarrange(rda.long.plot, venn.plot.long.euler,
                            ncol = 2, labels = c("A.", "B."),
                            widths = c(1.3,1))
figure.long <- ggarrange(figure.long.up, cluster.acclim.plot,
                         nrow = 2, labels = c("", "C.")) +
  theme(plot.margin = margin(.2, .2, .2, .2, "cm"), 
        plot.background = element_rect(fill = "white")) 

figure.long

# ggsave(file.path(here::here(), "02_Results", new_dir, "Fig_LongTerm.png"),
#        figure.long,
#        width = 6, height = 4.5, scale = 1.5)

ggsave(file.path(here::here(), "02_Results", new_dir, "Fig_LongTerm_new_RDA.png"),
       figure.long,
       width = 6.5, height = 4.5, scale = 1.8)

ggsave(file.path(here::here(), "02_Results", new_dir, "Fig2_LongTerm_new_RDA.pdf"),
       figure.long,
       width = 16.5, height = 11.5, units = "cm", scale = 1.8)

# write.csv(df.longterm,
#           file.path(here::here(), "02_Results", new_dir, "Longterm_Gene_cluster.csv"))


# DETs across temperature for short-term temperature change-----------------------------

metaData <- meta.short
metaData$sampling <- factor(paste(metaData$sampling))
counts <- subset(counts.all, select = metaData$ID)

levels(metaData$Delta_short) <- c("down", "no", "up")

# individual order
identical(colnames(counts), metaData$ID)

# subset normalize data
count.vst <- subset(assay(counts.vst.all), select = metaData$ID)

# Distance 

count.dist <- vegdist(t(count.vst), method = "euclidean")

# RDA

rda.short.c <- capscale(count.dist ~ Temperature_sampling * Temperature_acclim, metaData)
anova.cca(rda.short.c, by = "margin")
RsquareAdj(rda.short.c)

# rda.short.f <- capscale(count.dist ~ acclim * sampling, metaData) #metaData$sampling + Condition(metaData$acclim))
# anova.cca(rda.short.f, by = "margin", scope = "acclim")
# RsquareAdj(rda.short.f)

res.short.c <- data.frame(Effect = c("Sampling", "Acclimation", "Sampling:Accimation"),
                          rbind(anova.cca(rda.short.c, by = "margin", scope = "Temperature_sampling")[1,],
                                anova.cca(rda.short.c, by = "margin", scope = "Temperature_acclim")[1,],
                                anova.cca(rda.short.c, by = "margin")[1,]),
                          R2 = c(RsquareAdj(capscale(count.dist ~ Temperature_sampling + Condition(Temperature_acclim + Temperature_sampling:Temperature_acclim), metaData))$r.squared,
                                 RsquareAdj(capscale(count.dist ~ Temperature_acclim + Condition(Temperature_sampling + Temperature_sampling:Temperature_acclim), metaData))$r.squared,
                                 RsquareAdj(capscale(count.dist ~ Temperature_sampling:Temperature_acclim + Condition(Temperature_sampling + Temperature_acclim), metaData))$r.squared),
                          adjusted_R2 = c(RsquareAdj(capscale(count.dist ~ Temperature_sampling + Condition(Temperature_acclim + Temperature_sampling:Temperature_acclim), metaData))$adj.r.squared,
                                          RsquareAdj(capscale(count.dist ~ Temperature_acclim + Condition(Temperature_sampling + Temperature_sampling:Temperature_acclim), metaData))$adj.r.squared,
                                          RsquareAdj(capscale(count.dist ~ Temperature_sampling:Temperature_acclim + Condition(Temperature_sampling + Temperature_acclim), metaData))$adj.r.squared),
                          Factor = rep("Temperature as continous variable", 3))

# write.csv(res.short.c,
#           file = file.path(here::here(), "02_Results", new_dir, "Results_RDA_Short-term.csv"),
#           row.names = F)
# 
# res.short.f <- data.frame(Effect = c("Sampling", "Acclimation", "Sampling:Accimation"),
#                           rbind(anova.cca(rda.short.f, by = "margin", scope = "sampling")[1,],
#                                 anova.cca(rda.short.f, by = "margin", scope = "acclim")[1,],
#                                 anova.cca(rda.short.f, by = "margin")[1,]),
#                           R2 = c(RsquareAdj(capscale(count.dist ~ sampling + Condition(acclim + sampling:acclim), metaData))$r.squared,
#                                  RsquareAdj(capscale(count.dist ~ acclim + Condition(sampling + sampling:acclim), metaData))$r.squared,
#                                  RsquareAdj(capscale(count.dist ~ sampling:acclim + Condition(sampling + acclim), metaData))$r.squared),
#                           adjusted_R2 = c(RsquareAdj(capscale(count.dist ~ sampling + Condition(acclim + sampling:acclim), metaData))$adj.r.squared,
#                                           RsquareAdj(capscale(count.dist ~ acclim + Condition(sampling + sampling:acclim), metaData))$adj.r.squared,
#                                           RsquareAdj(capscale(count.dist ~ sampling:acclim + Condition(sampling + acclim), metaData))$adj.r.squared),
#                           Factor = rep("Temperature as factor", 3))


# using sex or not

rda.s.1 <- capscale(count.dist ~ Temperature_sampling * Temperature_acclim, metaData)
  
rda.s.2 <- capscale(count.dist ~ metaData$Temperature_acclim + metaData$Temperature_sampling + 
                      metaData$Temperature_acclim:metaData$Temperature_sampling + 
                      substring(metaData$Sex_genet, 1, 1))
anova.cca(rda.s.2, by = "margin")

capscale(count.dist ~ metaData$Temperature_acclim +
           Condition(metaData$Temperature_sampling + 
                      metaData$Temperature_acclim:metaData$Temperature_sampling + 
                      substring(metaData$Sex_genet, 1, 1))) %>%
# RsquareAdj()
  anova.cca

capscale(count.dist ~ metaData$Temperature_sampling +
           Condition(metaData$Temperature_acclim + 
                       metaData$Temperature_acclim:metaData$Temperature_sampling + 
                       substring(metaData$Sex_genet, 1, 1))) %>%
# RsquareAdj()
  anova.cca()

capscale(count.dist ~ metaData$Temperature_acclim:metaData$Temperature_sampling + 
           Condition(metaData$Temperature_acclim + 
                       metaData$Temperature_sampling +
                       substring(metaData$Sex_genet, 1, 1))) %>%
  # RsquareAdj()
  anova.cca()

capscale(count.dist ~ substring(metaData$Sex_genet, 1, 1) +
           Condition(metaData$Temperature_acclim + 
                       metaData$Temperature_sampling +
                       metaData$Temperature_acclim:metaData$Temperature_sampling)) %>%
  # RsquareAdj()
  anova.cca()

# PCA plot

## Principal Component Analyis - PCA ----
pca.short <- rda(t(count.vst))

pca.df.short <- data.frame(scores(pca.short, display = "sites"),
                          acclim = metaData$Temperature_acclim,
                          sampling = metaData$Temperature_sampling)

percentVar <- round(summary(pca.short)$cont$importance[2,1:2] * 100, 2)

pca.short.plot <- ggplot(pca.df.short, aes(PC1, PC2, 
                                         color = factor(sampling),
                                         shape = factor(acclim))) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey") +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(name = "Sampling\ntemperature (°C)",
                     values = c("#0571B0", "#92C5DE", "#F4A582", "#CA0020")) +
  scale_shape_manual(name = "Acclimation\ntemperature (°C)",
                     values = c(8, 17)) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank()
  ) +
  theme(plot.margin = margin(.5, .5, .5, .5, "cm")) 

pca.short.plot


## RDA plot

## Principal Component Analyis - PCA ----
rda.short <- capscale(count.dist ~ Temperature_sampling * Temperature_acclim, metaData)

rda.df.short <- data.frame(scores(rda.short, display = "sites"),
                           acclim = metaData$Temperature_acclim,
                           sampling = metaData$Temperature_sampling)

percentVar <- round(summary(rda.short)$cont$importance[2,1:2] * 100, 2)

rda.short.plot <- ggplot(rda.df.short, aes(CAP1, CAP2, 
                                           color = factor(sampling),
                                           shape = factor(acclim))) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey") +
  geom_point(size = 4) +
  xlab(paste0("RDA1: ", percentVar[1], "% variance")) +
  ylab(paste0("RDA2: ", percentVar[2], "% variance")) +
  scale_color_manual(name = "Short-term\ntemperature (°C)",
                     values = c("#0571B0", "#92C5DE", "#F4A582", "#CA0020")) +
  scale_shape_manual(name = "Long-term\ntemperature (°C)",
                     values = c(8, 17)) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank()
  ) +
  theme(plot.margin = margin(.5, .5, .5, .5, "cm")) 

rda.short.plot



## Result for all - wald test

design <- formula(~ acclim + acclim:Delta_short)

ddsObj <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metaData,
  design = design
)

ddsObj$Delta_short <- relevel(ddsObj$Delta_short, ref = "no")

dds.all <- DESeq(ddsObj)
resultsNames(dds.all)

# Results DETs

res.s5.down <- results(dds.all, name = "acclim5.Delta_shortdown")
res.s5.down.sig <- subset(res.s5.down, padj < padj.cutoff & abs(log2FoldChange) >= log2FC)
summary(res.s5.down.sig)

res.s5.up <- results(dds.all, name = "acclim5.Delta_shortup")
res.s5.up.sig <- subset(res.s5.up, padj < padj.cutoff & abs(log2FoldChange) >= log2FC)
summary(res.s5.up.sig)

res.s7.5.down <- results(dds.all, name = "acclim7.5.Delta_shortdown")
res.s7.5.down.sig <- subset(res.s7.5.down, padj < padj.cutoff & abs(log2FoldChange) >= log2FC)
summary(res.s7.5.down.sig)

res.s7.5.up <- results(dds.all, name = "acclim7.5.Delta_shortup")
res.s7.5.up.sig <- subset(res.s7.5.up, padj < padj.cutoff & abs(log2FoldChange) >= log2FC)
summary(res.s7.5.up.sig)

DET.shortterm <- rbind(res.s5.down %>% data.frame %>% mutate(Compare = "acclim5.Delta_shortdown"),
                       res.s5.up %>% data.frame %>% mutate(Compare = "acclim5.Delta_shortup"),
                       res.s7.5.down %>% data.frame %>% mutate(Compare = "acclim7.5.Delta_shortdown"),
                       res.s7.5.up %>% data.frame %>% mutate(Compare = "acclim7.5.Delta_shortup")
                       )
write.csv(DET.shortterm,
          file.path(here(), "02_Results", new_dir, "DETs_short_term.csv"))

## Venn plot
s1 <- list("5.0°C to 2.5°C" = row.names(res.s5.down.sig),
           "5.0°C to 7.5°C" = row.names(res.s5.up.sig),
           "7.5°C to 5.0°C" = row.names(res.s7.5.down.sig),
           "7.5°C to 10.0°C" = rownames(res.s7.5.up.sig)
)

venn.plot.short.euler <- plot(euler(s1, shape = "circle"),
                             quantities = list(T, cex = 0.75),
                             legend = list(T, cex = 0.75, side = "right"),
                             fill = c("#0571B0", "#F4A582", "blue4", "#CA0020"),
                             alpha = 0.5, edges = F) 
venn.plot.short.euler  


## Group genes expression ---------------------

deg_list.s <- unique(c(row.names(res.s5.down.sig),
                       row.names(res.s5.up.sig),
                       row.names(res.s7.5.down.sig),
                       rownames(res.s7.5.up.sig)))

deg_vst.s <- count.vst[deg_list.s, ]

meta.temp <- metaData
rownames(meta.temp) <- paste(metaData$Seq_batch, metaData$ID_sample, sep = ".")
meta.temp$acclim <- as.factor(meta.temp$Temperature_acclim)
meta.temp$sampling <- as.factor(meta.temp$Temperature_sampling)

clusters_short <- DEGreport::degPatterns(deg_vst.s, 
                                         metadata = meta.temp, 
                                         time = "sampling", col = "acclim",
                                         minc = 15
) 

df.shortterm <- clusters_short[["normalized"]] %>% 
  group_by(cluster) %>% 
  mutate(sum_cluster = length(unique(genes))) %>% 
  ungroup() %>% arrange(desc(sum_cluster)) %>%
  mutate(sum_cluster2 = factor(sum_cluster, levels = unique(sum_cluster))) %>%
  mutate(cluster2 = as.integer(sum_cluster2)) %>%
  mutate(cluster_new = paste("Group S", formatC(cluster2, width=2, flag="0"), " (", sum_cluster, ")", sep = "")) %>%
  data.frame()

cluster.short.plot <- ggplot(df.shortterm, aes(x = sampling, y = value, color = acclim, fill = acclim)) +
  geom_violin(alpha = 0.5, position = position_dodge(.5)) +
  geom_smooth(aes(group = acclim),
              linetype="dashed", method = "gam", formula = y ~ poly(x, 2)) +
  scale_fill_manual("Long-term (acclimation)\ntemperature",
                    values = c("#0571B0", "#F4A582"),
                    labels = c("5°C", "7.5°C" )) +
  scale_color_manual("Long-term (acclimation)\ntemperature",
                     values = c("#0571B0", "#F4A582"),
                     labels = c("5°C", "7.5°C" )) +
  xlab("Sampling temperature (°C)") +
  ylab("Z-score of transcript abundance") +
  theme_bw() +
  facet_wrap(~ cluster_new, ncol = 5) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  theme(plot.margin = margin(.5, .5, .5, .5, "cm"),
        legend.position = c(1, 0.05), legend.justification = c(1, 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10))

cluster.short.plot

### Export results short-term ---------------------------------------------------

figure.short.up <- ggarrange(rda.short.plot, venn.plot.short.euler,
                             ncol = 2, labels = c("A.", "B."),
                             widths = c(1.3,1))
figure.short <- ggarrange(figure.short.up, cluster.short.plot,
                          nrow = 2, labels = c("", "C.")) +
  theme(plot.margin = margin(.2, .2, .2, .2, "cm"), 
        plot.background = element_rect(fill = "white")) 

figure.short

# ggsave(file.path(here::here(), "02_Results", new_dir, "Fig_ShortTerm.png"),
#        figure.short,
#        width = 6, height = 4.5, scale = 1.5)
# 
# write.csv(df.shortterm,
#           file.path(here::here(), "02_Results", new_dir, "Shortterm_Gene_cluster.csv"))

ggsave(file.path(here::here(), "02_Results", new_dir, "Fig_ShortTerm_new_RDA.png"),
       figure.short,
       width = 6.5, height = 4.5, scale = 1.8)

ggsave(file.path(here::here(), "02_Results", new_dir, "Fig3_ShortTerm_new_RDA.pdf"),
       figure.short,
       width = 16.5, height = 11.5, units = "cm", scale = 1.8)

# Diff long vs short -----------------------------------------------------------

# list.DETs <- list("LT (5.0°C to 2.5°C)" = rownames(res.a.down.sig),
#                   "LT (5.0°C to 7.5°C)" = rownames(res.a.up.sig),
#                   "LT (5.0°C to 10.0°C)" = rownames(res.a.up2.sig),
#                   "ST (5.0°C to 2.5°C)" = row.names(res.s5.down.sig),
#                   "ST (5.0°C to 7.5°C)" = row.names(res.s5.up.sig),
#                   "ST (7.5°C to 5.0°C)" = row.names(res.s7.5.down.sig),
#                   "ST (7.5°C to 10.0°C)" = rownames(res.s7.5.up.sig))

list.DETs <- list("Long-term response" = unique(c(rownames(res.a.down.sig),
                                         rownames(res.a.up.sig) %>% unique,
                                         rownames(res.a.up2.sig) %>% unique)),
                  "Short-term response (fish acclimated at 5°C)" = unique(c(row.names(res.s5.down.sig),
                                               row.names(res.s5.up.sig))),
                  "Short-term response (fish acclimated at 7.5°C)" = unique(c(row.names(res.s7.5.down.sig),
                                               row.names(res.s7.5.up.sig))))

df.dets <- data.frame(Transcrip_name = c(list.DETs$`Long-term response`,
                                         list.DETs$`Short-term response (fish acclimated at 5°C)`,
                                         list.DETs$`Short-term response (fish acclimated at 7.5°C)`),
                      Compare = c(rep("LongTerm", length(list.DETs$`Long-term response`)),
                                  rep("ShortTerm_5", length(list.DETs$`Short-term response (fish acclimated at 5°C)`)),
                                  rep("ShortTerm_7.5", length(list.DETs$`Short-term response (fish acclimated at 7.5°C)`))
                                  ))

write.csv(df.dets,
          file.path(here(), "02_Results", new_dir, "List_DETs.csv"))

venn.all <- plot(euler(list.DETs, shape = "circle"),
                 quantities = list(T, cex = 1),
                 legend = list(T, cex = 0.75, side = "right"),
                # fill =  hcl.colors(7, "RdYlBu"),
                 fill = c("#CA0020", "#0571B0", "blue4"),
                 alpha = 0.75, edges = F) 
venn.all



# GO enrichment - Using TopGO --------------------------------------------------

annot.file <- read.csv("00_Data/03b_Annotations/annot_df.csv", header = T, row.names = 1)

# Number of transcript successfully annotated
annot.file$transcript_ID %>% unique() %>% length()
annot.file$transcript_ID %>% unique() %>% length() /35483 # total n transcript

# Number of GO terms
annot.file$GO %>% unique() %>% length()



gene2GO <- tapply(annot.file$GO, annot.file$transcript_ID, function(x)x)

# Lists of transcripts
all.genes <- row.names(counts.all)
long.all <- unique(c(rownames(res.a.down.sig),
                     rownames(res.a.up.sig),
                     rownames(res.a.up2.sig)))
short.5 <- unique(c(rownames(res.s5.down.sig),
                      rownames(res.s5.up.sig)))
short.75 <- unique(c(rownames(res.s7.5.down.sig),
                    rownames(res.s7.5.up.sig)))

long.excl <- setdiff(long.all, c(short.5, short.75))
short.excl.5 <- setdiff(short.5, long.all)
short.excl.75 <- setdiff(short.75, long.all)

short.long.comm <- intersect(c(short.5, short.75), long.all)

## Enrichment analyses ----
go2gene.total <- NULL

### For exclusively long term response  ----

res.GO.enrich.long<- mat.or.vec(0,0)

# Define vector that is 1 if gene belong to a given meta-modules and 0 otherwise
tmp <- ifelse(all.genes %in% long.excl, 1, 0)
geneList <- tmp

# geneList needs names that match those for GO terms,
names(geneList) <- sapply(all.genes, function(x)x[1])

ontology <- c("BP", "CC", "MF")

for (i in 1:3) {  
  # Create topGOdata object:
  GOdata <- new("topGOdata",
                ontology = ontology[i],
                allGenes = geneList,
                geneSelectionFun = function(x)(x == 1),
                annot = annFUN.gene2GO, gene2GO = gene2GO)
  
  # Run Fisher’s Exact Test:
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  tab <- GenTable(GOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score),
                  numChar = 1000)
  tab$padj <- p.adjust(tab$raw.p.value, method = "fdr")
  
  # Retrieve transcript id related to a GO.ID
  allGO = genesInTerm(GOdata)
  go2gene.total <- c(go2gene.total, allGO)
  
  ## - create table results
  rm(temp)
  temp <- tab
  #  temp <- subset(tab, padj <= 0.1)
  temp$ontology <- ontology[i]
  
  ## res.GO.enrich.pop <- mat.or.vec(0,0)
  res.GO.enrich.long <- rbind(res.GO.enrich.long, temp)
}

### For exclusively short term response at 5C ----

res.GO.enrich.short.5 <- mat.or.vec(0,0)

# Define vector that is 1 if gene belong to a given meta-modules and 0 otherwise
tmp <- ifelse(all.genes %in% short.excl.5, 1, 0)
geneList <- tmp

# geneList needs names that match those for GO terms,
names(geneList) <- unlist(all.genes, function(x)x[1])

ontology <- c("BP", "CC", "MF")

for (i in 1:3) {  
  # Create topGOdata object:
  GOdata <- new("topGOdata",
                ontology = ontology[i],
                allGenes = geneList,
                geneSelectionFun = function(x)(x == 1),
                annot = annFUN.gene2GO, gene2GO = gene2GO)
  
  # Run Fisher’s Exact Test:
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  tab <- GenTable(GOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score),
                  numChar = 1000)
  tab$padj <- p.adjust(tab$raw.p.value, method = "fdr")
  
  # Retrieve transcript id related to a GO.ID
  allGO = genesInTerm(GOdata)
  go2gene.total <- c(go2gene.total, allGO)
  
  ## - create table results
  rm(temp)
  temp <- tab
  #  temp <- subset(tab, padj <= 0.1)
  temp$ontology <- ontology[i]
  
  ## res.GO.enrich.pop <- mat.or.vec(0,0)
  res.GO.enrich.short.5 <- rbind(res.GO.enrich.short.5, temp)
}

### For exclusively short term response at 7.5C ----

res.GO.enrich.short.75 <- mat.or.vec(0,0)

# Define vector that is 1 if gene belong to a given meta-modules and 0 otherwise
tmp <- ifelse(all.genes %in% short.excl.75, 1, 0)
geneList <- tmp

# geneList needs names that match those for GO terms,
names(geneList) <- unlist(all.genes, function(x)x[1])

ontology <- c("BP", "CC", "MF")

for (i in 1:3) {  
  # Create topGOdata object:
  GOdata <- new("topGOdata",
                ontology = ontology[i],
                allGenes = geneList,
                geneSelectionFun = function(x)(x == 1),
                annot = annFUN.gene2GO, gene2GO = gene2GO)
  
  # Run Fisher’s Exact Test:
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  tab <- GenTable(GOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score),
                  numChar = 1000)
  tab$padj <- p.adjust(tab$raw.p.value, method = "fdr")
  
  # Retrieve transcript id related to a GO.ID
  allGO = genesInTerm(GOdata)
  go2gene.total <- c(go2gene.total, allGO)
  
  ## - create table results
  rm(temp)
  temp <- tab
  #  temp <- subset(tab, padj <= 0.1)
  temp$ontology <- ontology[i]
  
  ## res.GO.enrich.pop <- mat.or.vec(0,0)
  res.GO.enrich.short.75 <- rbind(res.GO.enrich.short.75, temp)
}

### For common response  ----

res.GO.enrich.comm <- mat.or.vec(0,0)

# Define vector that is 1 if gene belong to a given meta-modules and 0 otherwise
tmp <- ifelse(all.genes %in% short.long.comm, 1, 0)
geneList <- tmp

# geneList needs names that match those for GO terms,
names(geneList) <- unlist(all.genes, function(x)x[1])

ontology <- c("BP", "CC", "MF")

for (i in 1:3) {  
  # Create topGOdata object:
  GOdata <- new("topGOdata",
                ontology = ontology[i],
                allGenes = geneList,
                geneSelectionFun = function(x)(x == 1),
                annot = annFUN.gene2GO, gene2GO = gene2GO)
  
  # Run Fisher’s Exact Test:
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  tab <- GenTable(GOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score),
                  numChar = 1000)
  tab$padj <- p.adjust(tab$raw.p.value, method = "fdr")
  
  # Retrieve transcript id related to a GO.ID
  allGO = genesInTerm(GOdata)
  go2gene.total <- c(go2gene.total, allGO)
  
  ## - create table results
  rm(temp)
  temp <- tab
  #  temp <- subset(tab, padj <= 0.1)
  temp$ontology <- ontology[i]
  
  ## res.GO.enrich.pop <- mat.or.vec(0,0)
  res.GO.enrich.comm <- rbind(res.GO.enrich.comm, temp)
}

write.csv(res.GO.enrich.long,
          file.path(here::here(), "02_Results", new_dir, "Go_Enrich_LongTermExclu.csv"))

write.csv(res.GO.enrich.short.5,
          file.path(here::here(), "02_Results", new_dir, "Go_Enrich_ShortTermExclu_5C.csv"))

write.csv(res.GO.enrich.short.75,
          file.path(here::here(), "02_Results", new_dir, "Go_Enrich_ShortTermExclu_7.5C.csv"))

write.csv(res.GO.enrich.comm,
          file.path(here::here(), "02_Results", new_dir, "Go_Enrich_Comm_LongShort.csv"))



## Detail for longterm cluster -------------------------------------------------


cluster.lt <- read.csv(file.path(here::here(), "02_Results", new_dir,
                                 "Longterm_Gene_cluster.csv"),
                       row.names = 1)

### 1) Cluster L01--------------------------------------------------------------
grp.l01 <- subset(cluster.lt, cluster2 == 1)$genes %>% unique()

res.GO.enrich.l01 <- mat.or.vec(0,0)

# Define vector that is 1 if gene belong to a given meta-modules and 0 otherwise
tmp <- ifelse(all.genes %in% grp.l01, 1, 0)
geneList <- tmp

# geneList needs names that match those for GO terms,
names(geneList) <- unlist(all.genes, function(x)x[1])

ontology <- c("BP", "CC", "MF")

for (i in 1:3) {  
  # Create topGOdata object:
  GOdata <- new("topGOdata",
                ontology = ontology[i],
                allGenes = geneList,
                geneSelectionFun = function(x)(x == 1),
                annot = annFUN.gene2GO, gene2GO = gene2GO)
  
  # Run Fisher’s Exact Test:
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  tab <- GenTable(GOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score),
                  numChar = 1000)
  tab$padj <- p.adjust(tab$raw.p.value, method = "fdr")
  
  # Retrieve transcript id related to a GO.ID
  allGO = genesInTerm(GOdata)
  go2gene.total <- c(go2gene.total, allGO)
  
  ## - create table results
  rm(temp)
  temp <- tab
  #  temp <- subset(tab, padj <= 0.1)
  temp$ontology <- ontology[i]
  
  ## res.GO.enrich.pop <- mat.or.vec(0,0)
  res.GO.enrich.l01 <- rbind(res.GO.enrich.l01, temp)
}

write.csv(res.GO.enrich.l01,
          file.path(here::here(), "02_Results", new_dir, "Go_Enrich_long_L01.csv"))


### 2) Cluster L02--------------------------------------------------------------
grp.l02 <- subset(cluster.lt, cluster2 == 2)$genes %>% unique()

res.GO.enrich.l02 <- mat.or.vec(0,0)

# Define vector that is 1 if gene belong to a given meta-modules and 0 otherwise
tmp <- ifelse(all.genes %in% grp.l02, 1, 0)
geneList <- tmp

# geneList needs names that match those for GO terms,
names(geneList) <- unlist(all.genes, function(x)x[1])

ontology <- c("BP", "CC", "MF")

for (i in 1:3) {  
  # Create topGOdata object:
  GOdata <- new("topGOdata",
                ontology = ontology[i],
                allGenes = geneList,
                geneSelectionFun = function(x)(x == 1),
                annot = annFUN.gene2GO, gene2GO = gene2GO)
  
  # Run Fisher’s Exact Test:
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  tab <- GenTable(GOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score),
                  numChar = 1000)
  tab$padj <- p.adjust(tab$raw.p.value, method = "fdr")
  
  # Retrieve transcript id related to a GO.ID
  allGO = genesInTerm(GOdata)
  go2gene.total <- c(go2gene.total, allGO)
  
  ## - create table results
  rm(temp)
  temp <- tab
  #  temp <- subset(tab, padj <= 0.1)
  temp$ontology <- ontology[i]
  
  ## res.GO.enrich.pop <- mat.or.vec(0,0)
  res.GO.enrich.l02 <- rbind(res.GO.enrich.l02, temp)
}

write.csv(res.GO.enrich.l02,
          file.path(here::here(), "02_Results", new_dir, "Go_Enrich_long_L02.csv"))

# Reduce + Visualize Gene Ontology ---------------------------------------------
library(rrvgo)

# subset significantly enriched GO terms and paste to REVIGO

padj.cutoff.go <- 0.001

go.long <- subset(res.GO.enrich.long, padj <= padj.cutoff.go)
go.short5 <- subset(res.GO.enrich.short.5, padj <= padj.cutoff.go)
go.short75 <- subset(res.GO.enrich.short.75, padj <= padj.cutoff.go)
go.com <- subset(res.GO.enrich.comm, padj <= padj.cutoff.go)

go.l01 <- subset(res.GO.enrich.l01, padj <= padj.cutoff.go)
go.l02 <- subset(res.GO.enrich.l02, padj <= padj.cutoff.go)


sim.threshold <- 0.9 ## similarity threshold - allowed large similarity

## Long-term response to temperature -------------------------------------------

GO.df <- go.long

# Molecular function
simMatrix <- calculateSimMatrix(GO.df$GO.ID,
                                orgdb="org.Dr.eg.db", ## using Danio rerio
                                ont="MF",
                                method="Wang")
scores <- setNames(-log10(GO.df$padj), GO.df$GO.ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=sim.threshold,
                                orgdb="org.Dr.eg.db")

write.csv(reducedTerms,
          file.path(here::here(), "02_Results", new_dir, "ReducedTerms_LongTerm_MF.csv"))

LongTerm.MF.plot <- scatterPlot(simMatrix, reducedTerms) +
  ggtitle("Molecular Function") + 
  theme(plot.title = element_text(size=12)) 

pdf(file = file.path(here::here(), "02_Results", new_dir, "LongTerm_MF_treemap.pdf"),
    width = 4, height = 6)
treemapPlot(reducedTerms,
            title = "Long-term response - Molecular function")
dev.off()


# Biological processes
simMatrix <- calculateSimMatrix(GO.df$GO.ID,
                                orgdb="org.Dr.eg.db", ## using Danio rerio
                                ont="BP",
                                method="Wang")
scores <- setNames(-log10(GO.df$padj), GO.df$GO.ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=sim.threshold,
                                orgdb="org.Dr.eg.db")
write.csv(reducedTerms,
          file.path(here::here(), "02_Results", new_dir, "ReducedTerms_LongTerm_BP.csv"))


LongTerm.BP.plot <- scatterPlot(simMatrix, reducedTerms) +
  ggtitle("Biological Process") + 
  theme(plot.title = element_text(size=12)) 

pdf(file = file.path(here::here(), "02_Results", new_dir, "LongTerm_BP_treemap.pdf"),
    width = 4, height = 6)
treemapPlot(reducedTerms,
            title = "Long-term response - Biological Process")
dev.off()

# Cellular component 
simMatrix <- calculateSimMatrix(GO.df$GO.ID,
                                orgdb="org.Dr.eg.db", ## using Danio rerio
                                ont="CC",
                                method="Wang")
scores <- setNames(-log10(GO.df$padj), GO.df$GO.ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=sim.threshold,
                                orgdb="org.Dr.eg.db")

write.csv(reducedTerms,
          file.path(here::here(), "02_Results", new_dir, "ReducedTerms_LongTerm_CC.csv"))

LongTerm.CC.plot <- scatterPlot(simMatrix, reducedTerms) +
  ggtitle("Cellular component") + 
  theme(plot.title = element_text(size=12)) 

pdf(file = file.path(here::here(), "02_Results", new_dir, "LongTerm_CC_treemap.pdf"),
    width = 4, height = 6)
treemapPlot(reducedTerms,
            title = "Long-term response - Cellular component")
dev.off()

## Short term response to temperature ------------------------------------------

### acclim 5C ------------------------------------------------------------------

rm(list = c("scores", "reducedTerms", "GO.df"))

GO.df <- go.short5

# Molecular function
simMatrix <- calculateSimMatrix(GO.df$GO.ID,
                                orgdb="org.Dr.eg.db", ## using Danio rerio
                                ont="MF",
                                method="Wang")
scores <- setNames(-log10(GO.df$padj), GO.df$GO.ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=sim.threshold,
                                orgdb="org.Dr.eg.db")

write.csv(reducedTerms,
          file.path(here::here(), "02_Results", new_dir, "ReducedTerms_ShortTerm_5_MF.csv"))

ShortTerm.5.MF.plot <- scatterPlot(simMatrix, reducedTerms) +
  ggtitle("Molecular Function") + 
  theme(plot.title = element_text(size=12)) 

pdf(file = file.path(here::here(), "02_Results", new_dir, "ShortTerm_5_MF_treemap.pdf"),
    width = 4, height = 6)
treemapPlot(reducedTerms,
            title = "5C Short-term response - Molecular Function")
dev.off()

# Biological processes
simMatrix <- calculateSimMatrix(GO.df$GO.ID,
                                orgdb="org.Dr.eg.db", ## using Danio rerio
                                ont="BP",
                                method="Wang")
scores <- setNames(-log10(GO.df$padj), GO.df$GO.ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=sim.threshold,
                                orgdb="org.Dr.eg.db")
write.csv(reducedTerms,
          file.path(here::here(), "02_Results", new_dir, "ReducedTerms_ShortTerm_5_BP.csv"))


ShortTerm.5.BP.plot <- scatterPlot(simMatrix, reducedTerms) +
  ggtitle("Biological Process") + 
  theme(plot.title = element_text(size=12)) 

pdf(file = file.path(here::here(), "02_Results", new_dir, "ShortTerm_5_BP_treemap.pdf"),
    width = 4, height = 6)
treemapPlot(reducedTerms,
            title = "Short-term response - Biological Process")
dev.off()

# Cellular component 
simMatrix <- calculateSimMatrix(GO.df$GO.ID,
                                orgdb="org.Dr.eg.db", ## using Danio rerio
                                ont="CC",
                                method="Wang")

## No terms were found in orgdb for CC 


### acclim 7.5C ------------------------------------------------------------------

rm(list = c("scores", "reducedTerms", "GO.df"))

GO.df <- go.short75

# Molecular function
simMatrix <- calculateSimMatrix(GO.df$GO.ID,
                                orgdb="org.Dr.eg.db", ## using Danio rerio
                                ont="MF",
                                method="Wang")
scores <- setNames(-log10(GO.df$padj), GO.df$GO.ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=sim.threshold,
                                orgdb="org.Dr.eg.db")

write.csv(reducedTerms,
          file.path(here::here(), "02_Results", new_dir, "ReducedTerms_ShortTerm_75_MF.csv"))

ShortTerm.75.MF.plot <- scatterPlot(simMatrix, reducedTerms) +
  ggtitle("Molecular Function") + 
  theme(plot.title = element_text(size=12)) 

pdf(file = file.path(here::here(), "02_Results", new_dir, "ShortTerm_75_MF_treemap.pdf"),
    width = 4, height = 6)
treemapPlot(reducedTerms,
            title = "7.5C Short-term response - Molecular Function")
dev.off()

# Biological processes
simMatrix <- calculateSimMatrix(GO.df$GO.ID,
                                orgdb="org.Dr.eg.db", ## using Danio rerio
                                ont="BP",
                                method="Wang")
scores <- setNames(-log10(GO.df$padj), GO.df$GO.ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=sim.threshold,
                                orgdb="org.Dr.eg.db")
write.csv(reducedTerms,
          file.path(here::here(), "02_Results", new_dir, "ReducedTerms_ShortTerm_75_BP.csv"))


ShortTerm.75.BP.plot <- scatterPlot(simMatrix, reducedTerms) +
  ggtitle("Biological Process") + 
  theme(plot.title = element_text(size=12)) 

pdf(file = file.path(here::here(), "02_Results", new_dir, "ShortTerm_75_BP_treemap.pdf"),
    width = 4, height = 6)
treemapPlot(reducedTerms,
            title = "7.5-Short-term response - Biological Process")
dev.off()

# Cellular component 
simMatrix <- calculateSimMatrix(GO.df$GO.ID,
                                orgdb="org.Dr.eg.db", ## using Danio rerio
                                ont="CC",
                                method="Wang")

scores <- setNames(-log10(GO.df$padj), GO.df$GO.ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=sim.threshold,
                                orgdb="org.Dr.eg.db")

write.csv(reducedTerms,
          file.path(here::here(), "02_Results", new_dir, "ReducedTerms_ShortTerm_75_CC.csv"))

ShortTerm.75.CC.plot <- scatterPlot(simMatrix, reducedTerms) +
  ggtitle("Cellular component") + 
  theme(plot.title = element_text(size=12)) 

pdf(file = file.path(here::here(), "02_Results", new_dir, "ShortTerm_75_CC_treemap.pdf"),
    width = 4, height = 6)
treemapPlot(reducedTerms,
            title = "Cellular component")
dev.off()



## Common response -------------------------------------------------------------
GO.df <- go.com

# Molecular function
simMatrix <- calculateSimMatrix(GO.df$GO.ID,
                                orgdb="org.Dr.eg.db", ## using Danio rerio
                                ont="MF",
                                method="Wang")
scores <- setNames(-log10(GO.df$padj), GO.df$GO.ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=sim.threshold,
                                orgdb="org.Dr.eg.db")

write.csv(reducedTerms,
          file.path(here::here(), "02_Results", new_dir, "ReducedTerms_Common_MF.csv"))

Common.MF.plot <- scatterPlot(simMatrix, reducedTerms) +
  ggtitle("Molecular Function") + 
  theme(plot.title = element_text(size=12)) 

pdf(file = file.path(here::here(), "02_Results", new_dir, "Comm_MF_treemap.pdf"),
    width = 4, height = 6)
treemapPlot(reducedTerms,
            title = "Common response - Molecular Function")
dev.off()


# Biological processes
simMatrix <- calculateSimMatrix(GO.df$GO.ID,
                                orgdb="org.Dr.eg.db", ## using Danio rerio
                                ont="BP",
                                method="Wang")
scores <- setNames(-log10(GO.df$padj), GO.df$GO.ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=sim.threshold,
                                orgdb="org.Dr.eg.db")

write.csv(reducedTerms,
          file.path(here::here(), "02_Results", new_dir, "ReducedTerms_Common_BP.csv"))

Common.BP.plot <- scatterPlot(simMatrix, reducedTerms) +
  ggtitle("Biological Process") + 
  theme(plot.title = element_text(size=12)) 


pdf(file = file.path(here::here(), "02_Results", new_dir, "Comm_BP_treemap.pdf"),
    width = 4, height = 6)
treemapPlot(reducedTerms,
            title = "Common response - Biological Process")
dev.off()


# Cellular component 
simMatrix <- calculateSimMatrix(GO.df$GO.ID,
                                orgdb="org.Dr.eg.db", ## using Danio rerio
                                ont="CC",
                                method="Wang")

## No terms were found in orgdb for CC



### for L01 and L02 ------------------------------------------------------------
## Long-term response to temperature L01

GO.df <- go.l01

# Molecular function
simMatrix <- calculateSimMatrix(GO.df$GO.ID,
                                orgdb="org.Dr.eg.db", ## using Danio rerio
                                ont="MF",
                                method="Wang")
scores <- setNames(-log10(GO.df$padj), GO.df$GO.ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=sim.threshold,
                                orgdb="org.Dr.eg.db")

write.csv(reducedTerms,
          file.path(here::here(), "02_Results", new_dir, "ReducedTerms_L01_MF.csv"))

L01.MF.plot <- scatterPlot(simMatrix, reducedTerms) +
  ggtitle("Molecular Function") + 
  theme(plot.title = element_text(size=12)) 

pdf(file = file.path(here::here(), "02_Results", new_dir, "L01_MF_treemap.pdf"),
    width = 4, height = 6)
treemapPlot(reducedTerms,
            title = "Group L01 response - Molecular function")
dev.off()


# Biological processes
simMatrix <- calculateSimMatrix(GO.df$GO.ID,
                                orgdb="org.Dr.eg.db", ## using Danio rerio
                                ont="BP",
                                method="Wang")
scores <- setNames(-log10(GO.df$padj), GO.df$GO.ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=sim.threshold,
                                orgdb="org.Dr.eg.db")
write.csv(reducedTerms,
          file.path(here::here(), "02_Results", new_dir, "ReducedTerms_L01_BP.csv"))


L01.BP.plot <- scatterPlot(simMatrix, reducedTerms) +
  ggtitle("Biological Process") + 
  theme(plot.title = element_text(size=12)) 

pdf(file = file.path(here::here(), "02_Results", new_dir, "L01_BP_treemap.pdf"),
    width = 4, height = 6)
treemapPlot(reducedTerms,
            title = "L01 response - Biological Process")
dev.off()

# Cellular component 
simMatrix <- calculateSimMatrix(GO.df$GO.ID,
                                orgdb="org.Dr.eg.db", ## using Danio rerio
                                ont="CC",
                                method="Wang")
scores <- setNames(-log10(GO.df$padj), GO.df$GO.ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=sim.threshold,
                                orgdb="org.Dr.eg.db")

write.csv(reducedTerms,
          file.path(here::here(), "02_Results", new_dir, "ReducedTerms_L01_CC.csv"))

L01.CC.plot <- scatterPlot(simMatrix, reducedTerms) +
  ggtitle("Cellular component") + 
  theme(plot.title = element_text(size=12)) 

pdf(file = file.path(here::here(), "02_Results", new_dir, "L01_CC_treemap.pdf"),
    width = 4, height = 6)
treemapPlot(reducedTerms,
            title = "L01 response - Cellular component")
dev.off()



GO.df <- go.l02

# Molecular function
simMatrix <- calculateSimMatrix(GO.df$GO.ID,
                                orgdb="org.Dr.eg.db", ## using Danio rerio
                                ont="MF",
                                method="Wang")
scores <- setNames(-log10(GO.df$padj), GO.df$GO.ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=sim.threshold,
                                orgdb="org.Dr.eg.db")

write.csv(reducedTerms,
          file.path(here::here(), "02_Results", new_dir, "ReducedTerms_L02_MF.csv"))

L02.MF.plot <- scatterPlot(simMatrix, reducedTerms) +
  ggtitle("Molecular Function") + 
  theme(plot.title = element_text(size=12)) 

pdf(file = file.path(here::here(), "02_Results", new_dir, "L02_MF_treemap.pdf"),
    width = 4, height = 6)
treemapPlot(reducedTerms,
            title = "Group L02 response - Molecular function")
dev.off()


# Biological processes
simMatrix <- calculateSimMatrix(GO.df$GO.ID,
                                orgdb="org.Dr.eg.db", ## using Danio rerio
                                ont="BP",
                                method="Wang")
scores <- setNames(-log10(GO.df$padj), GO.df$GO.ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=sim.threshold,
                                orgdb="org.Dr.eg.db")
write.csv(reducedTerms,
          file.path(here::here(), "02_Results", new_dir, "ReducedTerms_L02_BP.csv"))


L02.BP.plot <- scatterPlot(simMatrix, reducedTerms) +
  ggtitle("Biological Process") + 
  theme(plot.title = element_text(size=12)) 

pdf(file = file.path(here::here(), "02_Results", new_dir, "L02_BP_treemap.pdf"),
    width = 4, height = 6)
treemapPlot(reducedTerms,
            title = "L02 response - Biological Process")
dev.off()

# Cellular component 
simMatrix <- calculateSimMatrix(GO.df$GO.ID,
                                orgdb="org.Dr.eg.db", ## using Danio rerio
                                ont="CC",
                                method="Wang")
scores <- setNames(-log10(GO.df$padj), GO.df$GO.ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=sim.threshold,
                                orgdb="org.Dr.eg.db")

write.csv(reducedTerms,
          file.path(here::here(), "02_Results", new_dir, "ReducedTerms_L02_CC.csv"))

L02.CC.plot <- scatterPlot(simMatrix, reducedTerms) +
  ggtitle("Cellular component") + 
  theme(plot.title = element_text(size=12)) 

pdf(file = file.path(here::here(), "02_Results", new_dir, "L02_CC_treemap.pdf"),
    width = 4, height = 6)
treemapPlot(reducedTerms,
            title = "L02 response - Cellular component")
dev.off()


## Export GO fig ----

GO.long.fig <- ggarrange(LongTerm.BP.plot, LongTerm.MF.plot, LongTerm.CC.plot,
                         ncol = 3, nrow = 1,
                         labels = LETTERS) +
  theme(plot.background = element_rect(fill = "white"))

GO.short.5.fig <- ggarrange(ShortTerm.5.BP.plot, ShortTerm.5.MF.plot, NULL,
                          ncol = 3, nrow = 1,
                          labels = LETTERS) +
  theme(plot.background = element_rect(fill = "white"))

GO.short.75.fig <- ggarrange(ShortTerm.75.BP.plot, ShortTerm.75.MF.plot, ShortTerm.75.CC.plot,
                            ncol = 3, nrow = 1,
                            labels = LETTERS) +
  theme(plot.background = element_rect(fill = "white"))



GO.comm.fig <- ggarrange(Common.BP.plot, Common.MF.plot, NULL,
                         ncol = 3, nrow = 1, labels = c("A", "B", "")) +
  theme(plot.background = element_rect(fill = "white"))

GO.long.fig
GO.short.5.fig
GO.short.75.fig
GO.comm.fig

ggsave(file.path(here::here(), "02_Results", new_dir, "GO.long.fig.png"),
       GO.long.fig,
       width = 6, height = 3, scale = 2.2)

ggsave(file.path(here::here(), "02_Results", new_dir, "GO.short.5.fig.png"),
       GO.short.5.fig,
       width = 6, height = 3, scale = 2.2)

ggsave(file.path(here::here(), "02_Results", new_dir, "GO.short.75.fig.png"),
       GO.short.75.fig,
       width = 6, height = 3, scale = 2.2)

ggsave(file.path(here::here(), "02_Results", new_dir, "GO.comm.fig.png"),
       GO.comm.fig,
       width = 6, height = 3, scale = 2.2)

ggsave(file.path(here::here(), "02_Results", new_dir, "Venn_all.png"),
       venn.all,
       width = 6, height = 3, scale = 1.75)

## Check for number of transcript involved in optic cup formation --------------

annot.file <- read.csv("00_Data/03b_Annotations/annot_df.csv", header = T, row.names = 1)

DETs.all <-  read.csv(file.path(here(), "02_Results", new_dir, "List_DETs.csv"))


df.common <- read.csv(file.path(here::here(), "02_Results", new_dir, "ReducedTerms_Common_BP.csv"),
                    row.names = 1) %>%
  mutate(Comparison = "SL_common")
df.LongTerm <- read.csv(file.path(here::here(), "02_Results", new_dir, "ReducedTerms_LongTerm_BP.csv"),
                    row.names = 1) %>%
  mutate(Comparison = "LongTerm")
df.short5 <- read.csv(file.path(here::here(), "02_Results", new_dir, "ReducedTerms_ShortTerm_5_BP.csv"),
                    row.names = 1) %>%
  mutate(Comparison = "Short5")
df.short7.5 <- read.csv(file.path(here::here(), "02_Results", new_dir, "ReducedTerms_ShortTerm_75_BP.csv"),
                      row.names = 1) %>%
  mutate(Comparison = "Short7.5")

df.temp <- rbind(df.common, df.LongTerm, df.short5, df.short7.5)

temp.2 <- subset(df.temp, parent == "GO:0003408")
temp.2$go %>% unique %>% length

annot.optic <- subset(annot.file, GO %in% temp.2$go)     

annot.dets <- subset(DETs.all, Transcrip_name %in% annot.optic$transcript_ID)
table(annot.dets$Compare)

annot.dets$Transcrip_name %>% unique %>% length
                
list.optic <- list(long = subset(annot.dets, Compare == "LongTerm") %>% pull(Transcrip_name),
                   short5 = subset(annot.dets, Compare == "ShortTerm_5") %>% pull(Transcrip_name),
                   shorty = subset(annot.dets, Compare == "ShortTerm_7.5") %>% pull(Transcrip_name)) 

plot(euler(list.optic, shape = "circle"),
     quantities = list(T, cex = 1),
     legend = list(T, cex = 0.75, side = "right"),
     # fill =  hcl.colors(7, "RdYlBu"),
     fill = c("#CA0020", "#0571B0", "blue4"),
     alpha = 0.75, edges = F) 



# comm.bp <- subset(res.GO.enrich.comm, ontology == "BP" & GO.ID %in% temp.2$go)
# dim(comm.bp)
# sum(comm.bp$Significant)


# ## Check for number of transcript involved in heat acclimation -----------------
# 
# 
# test <- subset(res.GO.enrich.long, ontology == "BP" & GO.ID ==  "GO:0008150")
# dim(test)
# 
# test <- subset(res.GO.enrich.short.75, ontology == "BP" & GO.ID ==  "GO:0009266")
# dim(test)
# 
# 
# df.temp.5 <- read.csv(file.path(here::here(), "02_Results", new_dir, "ReducedTerms_ShortTerm_5_BP.csv"), 
#                     row.names = 1)
# 
# temp.2 <- subset(df.temp.5, parent == "GO:0009266")
# 
# df.temp.75 <- read.csv(file.path(here::here(), "02_Results", new_dir, "ReducedTerms_ShortTerm_75_BP.csv"), 
#                       row.names = 1)
# 
# temp.2 <- subset(df.temp.75, parent == "GO:0009266")
# 
# 

# Specifically check for Heat-Shock-Protein ------------------------------------

# Expression level of HSP genes ------------------------------------------------

rm(list = ls())
library(tidyverse)


list.DETs <- read.csv("02_Results/07_DETs_results_clean/List_DETs.csv",
                      row.names = 1)


list.go.hsp <- c("GO:0070846",
                 "GO:1990565",
                 "GO:0051879",
                 "GO:0051008",
                 "GO:0030544",
                 "GO:0110078",
                 "GO:0001018"
                 # "GO:0030192", ## obsolete
                 # "GO:0030191",
                 # "GO:0008077",
                 # "GO:0010538",
                 # "GO:0010539",
                 # "GO:0010545",
                 # "GO:0010546" 
) 


## In annotation files 
annot.file <- read.csv("00_Data/03b_Annotations/annot_df.csv", header = T, row.names = 1)

### Transcript associated to HSP
hsp <- subset(annot.file, GO %in% list.go.hsp)

# Gene cluster
gene.cluster <- read.csv("02_Results/07_DETs_results_clean/Longterm_Gene_cluster.csv",
                         row.names = 1) %>%
  select(genes, cluster_new) %>%
  distinct %>%
  left_join(., hsp, by = c("genes" = "transcript_ID"))


## counts and metadata
meta.all <- read.csv(file.path(here::here(), "00_Data", "metaData_subset_20nov2024.csv"),
                     row.names = 1) 
counts.all <- read.csv(file.path(here::here(), "00_Data", "countsRNA_subset_20nov2024.csv"),
                       row.names = 1)
# Design
library(DESeq2)

design <- as.formula(~ 1)
ddsObj <- DESeqDataSetFromMatrix(
  countData = counts.all,
  colData = meta.all,
  design = design
)

counts.vst.all <- vst(ddsObj, blind = T)

# test <- t(assay(counts.vst.all))

hsp.vst <- counts.vst.all[c(hsp$transcript_ID)] %>% assay %>% data.frame %>%
  mutate(transcript_ID = row.names(.)) %>%
  pivot_longer(cols = -"transcript_ID") %>%
  mutate(Samples = gsub("X", "", name)) %>%
  left_join(., meta.all %>%
              mutate(Samples = paste(Seq_batch, ID_sample, sep = "."))) %>%
  mutate(DETs = ifelse(transcript_ID %in% list.DETs$Transcrip_name, "Yes", "No")) %>%
  left_join(., gene.cluster %>% select(genes, cluster_new) %>% distinct,
            by = c("transcript_ID" = "genes"))


# Number of transcript associated to HPS, but not DE
hsp.vst %>% filter(DETs == "No") %>% pull(transcript_ID) %>% unique %>% length

# Number of transcript associated to HPS, and DE
hsp.vst %>% filter(DETs == "Yes") %>% pull(transcript_ID) %>% unique %>% length

# Modules DETs associated to HSP
hsp.vst %>% filter(DETs == "Yes") %>% select(transcript_ID, cluster_new) %>%
  distinct %>% pull(cluster_new) %>% table

plot.HSP <- ggplot(hsp.vst %>% filter (Delta_short == 0) %>%
         mutate(group_names = recode(cluster_new, 
                                     "Group L01 (1490)" = "Group L01\n(12 DETs)",
                                     "Group L02 (1228)" = "Group L02\n(11 DETs)", 
                                     "Group L05 (180)" = "Group L05\n(2 DETs)",
                                     .missing = "Not differentially expressed\n(100 transcripts)")),
       aes(x = Temperature_acclim, y = value)) +
  geom_violin(aes(group = acclim), alpha = 0.5, fill = "grey") +
  geom_smooth(#aes(group = acclim),
    col = "red", 
    linetype="dashed", method = "gam", formula = y ~ poly(x, 2)) +
  scale_x_continuous(name = "Long-term temperature (°C)", breaks = c(2.5, 5, 7.5, 10)) +
  xlab("Long-term temperature (°C)") +
  ylab("Normalized transcript abundance") +
  theme_bw() +
  facet_wrap(~ group_names, ncol = 4, scales = "free_y") +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, vjust = 0)) +
  theme(plot.margin = margin(.5, .5, .5, .5, "cm")) 

plot.HSP

ggsave(file.path(here::here(), "02_Results", "07_DETs_results_clean", "Fig_HSP_genes.png"),
       plot.HSP,
       width = 6.5, height = 2.5, scale = 1.5)
