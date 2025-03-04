# Info -------------------------------------------------------------------------

# Individual sex effect on gene expression
# Differential gene expression analysis
# Sebastes faciatus in common garden


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

new_dir = "08_Sex_effect"
dir.create(file.path(here::here(), "02_Results", new_dir),
           showWarnings = FALSE)


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


## Subset samples sequenced in both batches ----

ID.duplicated <- (metaData.init %>% filter(duplicated(.[["ID_sample"]])))$ID_sample

metaData.batches <- subset(metaData.init, ID_sample %in% ID.duplicated)

counts.batches <- subset(counts.init, select = metaData.batches$ID)
counts.batches <- round(counts.batches, 0)

metaData.batches <- metaData.batches[match(colnames(counts.batches), metaData.batches$ID),]
metaData.batches$Seq_batch <- as.factor(metaData.batches$Seq_batch)


## subsample samples - remove duplicated individuals from the two sequencing batches ----
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

# Design
design <- as.formula(~ 1)
ddsObj <- DESeqDataSetFromMatrix(
  countData = counts.all,
  colData = meta.all,
  design = design
)

counts.vst.all <- vst(ddsObj, blind = T)

rm("metaData")


metaData <- meta.all
metaData$Sex <- substring(metaData$Sex_genet,1,1)

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

# RDA sex
rda.sex.temp <- capscale(count.dist ~ metaData$Sex * metaData$Temperature_acclim,
                         na.action = "na.omit")
anova.cca(rda.sex.temp, by = "margin")
anova.cca(rda.sex.temp, by = "terms")
RsquareAdj(capscale(count.dist ~ Sex:Temperature_acclim + Condition(Temperature_acclim + Sex), metaData))

rda.sex <- capscale(count.dist ~ metaData$Sex,
                        na.action = "na.omit")
anova.cca(rda.sex, by = "margin")
RsquareAdj(rda.sex)



# PCA plot

## Principal Component Analyis - PCA ----
pca.sex <- rda(t(count.vst))

pca.df.sex <- data.frame(scores(pca.sex, display = "sites"),
                         Sex_genet = metaData$Sex_genet,
                         Sex = metaData$Sex,
                         Temperature = metaData$Temperature_acclim %>% as.factor)

percentVar <- round(summary(pca.sex)$cont$importance[2,1:2] * 100, 2)

plot.sex.plot <- ggplot(pca.df.sex, aes(PC1, PC2, color = Temperature,
                                        shape = Sex)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey") +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(name = "Temperature (°C)",
                     values = c("#0571B0", "#92C5DE", "#F4A582", "#CA0020")) +
  scale_shape_manual(values = c(15, 17),
                     labels = c("Female", "Male")) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank()
  ) +
  theme(plot.margin = margin(.5, .5, .5, .5, "cm")) 

plot.sex.plot



## Genetic variation between sex -----------------------------------------------

## Retrieve genotype data from dd-RADseq analysis ------------------------------

SexGenotype <- read.csv(file.path(here::here(), "00_Data", "00_Data_Info",
                                  "SexGenotypes.csv"))

ref.ind <- read.table(file.path(here::here(), "00_Data", "00_Data_Info",
                              "FAS_4groups_samples.txt"))
names(ref.ind) <- c("ID_GQ", "POP")

pop.data.update <- read_csv(file = file.path(here::here(), "00_Data", "pop.data.csv")) 

ref.ind <- ref.ind %>% left_join(pop.data.update %>% dplyr::select(ID_GQ, Sexe_visuel, Longueur_mm))
ref.ind$POP <- "Natural population"
ref.ind <- subset(ref.ind, Longueur_mm > 220)

RNA.ind <- data.frame(ID_GQ = meta.all$Numero_unique_specimen,
                      Sexe_visuel = NA, POP = "Common garden",
                      Longueur_mm = NA) %>% 
  rbind(ref.ind)

df.gen.sex <- left_join(RNA.ind, SexGenotype)


## PCA on snps -----------------------------------------------------------------

## Missing value

minorAlleleCounts <- df.gen.sex[,5:10]

dataForPCA <- minorAlleleCounts %>%
  apply(MARGIN = 2, function(x){
    missing <- is.na(x)
    if (sum(missing) > 0) {
      x[missing] <- sample(x[!missing], sum(missing),replace = TRUE)
    }
    x
  })

row.names(dataForPCA) <- df.gen.sex$ID_GQ

res.pca.gen <- dataForPCA %>%
  rda()


pca.gen.sex <- data.frame(scores(res.pca.gen, display = "sites"),
                          Sex = df.gen.sex$Sexe_visuel,
                          Origin = df.gen.sex$POP)

percentVar <- round(summary(res.pca.gen)$cont$importance[2,1:2] * 100, 2)

plot.gen.plot <- ggplot(pca.gen.sex, aes(PC1, PC2, color = Origin, shape = Sex)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey") +
  geom_point(size = 4, alpha = 0.88) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(values = c("dodgerblue2", "chartreuse3")) +
  scale_shape_manual(values = c(17, 15), na.value = 18) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank()
  ) +
  theme(plot.margin = margin(.5, .5, .5, .5, "cm")) 

plot.gen.plot

plot.sex <- ggarrange(plot.gen.plot,
                      plot.sex.plot, 
                      nrow = 2,
                      labels = LETTERS)


ggsave(file.path(here::here(), "02_Results", new_dir, "PCA_sexEffect.png"),
       plot.sex,
       width = 5, height = 5, scale = 1.2)

