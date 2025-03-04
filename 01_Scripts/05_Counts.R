# Info --------------------------------------------------------------------

# featureCounts tool to get the gene/transcript counts
# 
# CL - May 2022

# Library -----------------------------------------------------------------

rm(list = ls())
gc()
library(here)
library(tidyverse)
library(vegan)
library(Biostrings)

# Counts ----

list.bam <- paste(file.path(here(), "00_Data", "02_Hisat2"),
                  list.files(file.path(here(), "00_Data", "02_Hisat2"), 
                             pattern = "hisat2.bam$"),
                  sep = "/")

annot.gtf <- paste(file.path(here(), "00_Data", "03b_StringTie_merge"),
                   list.files(file.path(here(), "00_Data", "03b_StringTie_merge"), 
                              pattern = ".gtf"),
                   sep = "/")

# count table
### /!\ choisir pour multi-mapping et overlap (-M -O --fraction) si nescessaire

# file.out <- paste(file.path(here(), "00_Data", "04_FeatureCounts"),
#                   "Count_transcript.txt", sep = "/")
# file.out2 <- paste(file.path(here(), "00_Data", "04_FeatureCounts"),
#                    "Count_gene.txt", sep = "/")

file.out <- paste(file.path(here(), "00_Data", "04_FeatureCounts"),
                  "Count_transcript_No_MultipleMapp_NO_overlapp.txt", sep = "/")

file.out2 <- paste(file.path(here(), "00_Data", "04_FeatureCounts"),
                   "Count_gene_No_MultipleMapp_NO_overlapp.txt", sep = "/")

cmd.transcript <- paste("-a", annot.gtf,
                        "-F", 'GTF',
                        "-g", 'transcript_id',
                        "-T", 20,
                        "-p", "--countReadPairs",
                        "-C", #  Do not count read pairs that have their two ends mapping to  different  chromosomes or mapping to same chromosome but on different strands.
                    #    "-M -O --fraction", 
                        "-o", file.out,
                        list.bam)

cmd.gene <- paste("-a", annot.gtf,
                  "-F", 'GTF',
                  "-g", 'gene_id',
                  "-T", 20,
                  "-p", "--countReadPairs",
                  "-C", #  Do not count read pairs that have their two ends mapping to  different  chromosomes or mapping to same chromosome but on different strands.
              #    "-M -O --fraction", 
                  "-o", file.out2,
                  list.bam)

system2("featureCounts", cmd.transcript)
system2("featureCounts", cmd.gene)


# check for count summary when not considering -M and -O

df.summary <- read.table("00_Data/04_FeatureCounts/Count_transcript_No_MultipleMapp_NO_overlapp.txt.summary",
                   header = T)
colnames(df.summary) <- colnames(df.summary) %>% 
  str_replace(., "X.media.genobiwan.Extra_Storage.Projets.REDTANK_Transcripto.00_Data.02_Hisat2.",
              "") %>%
  str_replace(., "_hisat2.bam", "")

df.summary.melted <- reshape2::melt(df.summary)

colourCount <- length(unique(df.summary.melted$Status))
getPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))

plot.fig <- ggplot(df.summary.melted, aes(x = value, fill = Status)) +
  geom_histogram(bins = 50) +
  facet_wrap(~ Status, scales = "free_y", ncol = 5) +
  theme_classic() +
  scale_fill_manual(values = getPalette(colourCount))+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(hjust = 0),
        legend.position = "none")
plot.fig

## 2 libraries showed more over-represented sequences than the other libraries
subset.df <- subset(df.summary.melted, Status == "Unassigned_MultiMapping")
head(subset.df[order(subset.df$value, decreasing = T),])

## 1828.ARN_22_00034 and 1828.ARN_22_00004

# Identification of the over-represented sequences
count.MO <- read.table("00_Data/04_FeatureCounts/Count_transcript.txt", header = T, row.names = 1)[,-(1:5)]
colnames(count.MO) <- colnames(count.MO) %>% 
  str_replace(., "X.media.genobiwan.Extra_Storage.Projets.REDTANK_Transcripto.00_Data.02_Hisat2.",
              "") %>%
  str_replace(., "_hisat2.bam", "")

#Hellinger transformation
count.MO.hell <-decostand(t(count.MO), method="hellinger") 

pca.count <- rda(count.MO.hell)
plot(pca.count, scaling = 1)

# % variation explained by the two first axis
round(100*(summary(pca.count)$cont$importance[2, 1:2]), 2)

# subst transcript overreprensetend in the two outlier libraries
score.pca.count <- data.frame(scores(pca.count, display = "species"))
score.pca.count <- score.pca.count[order(score.pca.count$PC1, decreasing = T),]
head(score.pca.count, n = 10)

transcript.outlier <- rownames(score.pca.count[1:2,])

# Blastn for over-represented sequences

fasta.transcript <- readDNAStringSet("00_Data/03b_StringTie_merge/Sfaciatus_transcript_seq.fasta")

sequence <- fasta.transcript[transcript.outlier[10]]
paste(sequence)

# seq #1: SEFA.16369.1 -> 28s Sebastes umbrosus
# seq #2: SEFA.16370.3 -> 28s Sebastes umbrosus
# seq #7: SEFA.17172.9 -> 28s Sebastes umbrosus
# seq #10: SEFA.17171.18 -> 28s Sebastes umbrosus

seq.temp <- fasta.transcript[transcript.outlier]
seq.temp

writeXStringSet(seq.temp,
                file.path(here(),
                          "00_Data",
                          "04_FeatureCounts",
                          "OverRepres_seq.fasta.gz"),
                append=FALSE,
                compress=T, compression_level=NA, format="fasta")

