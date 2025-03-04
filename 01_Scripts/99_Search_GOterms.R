# BiocManager::install("preprocessCore")
# BiocManager::install("impute")
# 

rm(list = ls())

library(tidyverse)

annot.file <- read.csv("00_Data/03b_Annotations/annot_df.csv", header = T, row.names = 1)

# GO terms associated to rhodopsin and opsin were retreived from GeneCards database (website)
rhodopsine <- subset(annot.file,GO %in% c("GO:0004930",
                                          "GO:0005502",
                                          "GO:0005515",
                                          "GO:0008020",
                                          "GO:0009881",
                                          "GO:0004930",
                                          "GO:0008020",
                                          "GO:0009881",
                                          "GO:0042802"))


LT.enrich <- read.csv("02_Results/07_DETs_results_clean/DETs_long_term.csv") %>%
  filter(padj < 0.05,
         abs(log2FoldChange) >= 1) 

ST.enrich <- read.csv("02_Results/07_DETs_results_clean/DETs_short_term.csv") %>%
  filter(padj < 0.05,
         abs(log2FoldChange) >= 1)


dets <- rbind(LT.enrich, ST.enrich)

# dets <- read.csv("02_Results/07_DETs_results_clean/List_DETs.csv")
# dets.eye <- subset(dets, Transcrip_name %in% rhodopsine$transcript_ID)
# 
# dets.eye$Transcrip_name %>% unique %>% length

dets.eye <- subset(dets, X %in% rhodopsine$transcript_ID)

dets.eye %>% group_by(Compare) %>%
  reframe(N_transcript = length(unique(X)))


xx <- read.csv("02_Results/07_DETs_results_clean/Go_Enrich_Comm_LongShort.csv") %>% 
  filter(padj <= 0.05) %>%
  filter(GO.ID %in% c("GO:0004930",
                      "GO:0005502",
                      "GO:0005515",
                      "GO:0008020",
                      "GO:0009881",
                      "GO:0004930",
                      "GO:0008020",
                      "GO:0009881",
                      "GO:0042802"))

yy <- read.csv("02_Results/07_DETs_results_clean/Go_Enrich_LongTermExclu.csv") %>% 
  filter(padj <= 0.05) %>%
  filter(GO.ID %in% c("GO:0004930",
                      "GO:0005502",
                      "GO:0005515",
                      "GO:0008020",
                      "GO:0009881",
                      "GO:0004930",
                      "GO:0008020",
                      "GO:0009881",
                      "GO:0042802"))

zz <- read.csv("02_Results/07_DETs_results_clean/Go_Enrich_ShortTermExclu_5C.csv")  %>% 
  filter(padj <= 0.05) %>%
  filter(GO.ID %in% c("GO:0004930",
                      "GO:0005502",
                      "GO:0005515",
                      "GO:0008020",
                      "GO:0009881",
                      "GO:0004930",
                      "GO:0008020",
                      "GO:0009881",
                      "GO:0042802"))

ww <- read.csv("02_Results/07_DETs_results_clean/Go_Enrich_ShortTermExclu_7.5C.csv") %>% 
  filter(padj <= 0.05) %>%
  filter(GO.ID %in% c("GO:0004930",
                      "GO:0005502",
                      "GO:0005515",
                      "GO:0008020",
                      "GO:0009881",
                      "GO:0004930",
                      "GO:0008020",
                      "GO:0009881",
                      "GO:0042802"))


## Search for GO terms associated to Heat-Shock_Proteins -----------------------

rm(list = ls())
library(tidyverse)
library(ggwordcloud)


list.go.hsp <- c("GO:0070846",
                 "GO:1990565",
                 "GO:0051879",
                 "GO:0051008",
                 "GO:0030544",
                 "GO:0110078",
                 "GO:0001018",
                 "GO:0030192", 
                 "GO:0030191",
                 "GO:0008077",
                 "GO:0010538",
                 "GO:0010539",
                 "GO:0010545",
                 "GO:0010546" 
) 


## In annotation files 
annot.file <- read.csv("00_Data/03b_Annotations/annot_df.csv", header = T, row.names = 1)

hsp <- subset(annot.file, GO %in% list.go.hsp)

## --> 125 transcripts associated to HSP


xx <- read.csv("02_Results/07_DETs_results_clean/Go_Enrich_Comm_LongShort.csv") %>% 
  filter(GO.ID %in% list.go.hsp) %>%
  filter(padj <= 0.05)

yy <- read.csv("02_Results/07_DETs_results_clean/Go_Enrich_LongTermExclu.csv")  %>% 
  filter(GO.ID %in% list.go.hsp) %>%
  filter(padj <= 0.05)

zz <- read.csv("02_Results/07_DETs_results_clean/Go_Enrich_ShortTermExclu_5C.csv")  %>% 
  filter(GO.ID %in% list.go.hsp) %>%
  filter(padj <= 0.05)

ww <- read.csv("02_Results/07_DETs_results_clean/Go_Enrich_ShortTermExclu_7.5C.csv")  %>% 
  filter(GO.ID %in% list.go.hsp) %>%
  filter(padj <= 0.05)

l01 <- read.csv("02_Results/07_DETs_results_clean/Go_Enrich_long_L01.csv")  %>% 
  filter(GO.ID %in% list.go.hsp) %>%
  filter(padj <= 0.05)

l02 <- read.csv("02_Results/07_DETs_results_clean/Go_Enrich_long_L02.csv")  %>% 
  filter(GO.ID %in% list.go.hsp) %>%
  filter(padj <= 0.05)





## Presence in DETs ?

LT.dets <- read.csv("02_Results/07_DETs_results_clean/DETs_long_term.csv") %>%
  filter(padj < 0.05,
         abs(log2FoldChange) >= 1) 

ST.dets <- read.csv("02_Results/07_DETs_results_clean/DETs_short_term.csv") %>%
  filter(padj < 0.05,
         abs(log2FoldChange) >= 1)


dets <- rbind(LT.dets, ST.dets)

dets.hsp <- subset(dets, X %in% hsp$transcript_ID)

## --> HSP seems to be absent from DETs... but detection in GO enrich 
##     (because of parent GO terms ?)



# Search for specific word 

df.test <- read.csv("02_Results/07_DETs_results_clean/Go_Enrich_long_L02.csv") %>% 
  filter(padj <= 0.05)

df.test <- read.csv("02_Results/07_DETs_results_clean/Go_Enrich_long_L01.csv") %>% 
  filter(padj <= 0.05)

df.test <- read.csv("02_Results/07_DETs_results_clean/Go_Enrich_Comm_LongShort.csv") %>% 
  filter(padj <= 0.05)

df.test <- read.csv("02_Results/07_DETs_results_clean/Go_Enrich_LongTermExclu.csv") %>% 
  filter(padj <= 0.05)

df.test <- read.csv("02_Results/07_DETs_results_clean/Go_Enrich_ShortTermExclu_5C.csv") %>% 
  filter(padj <= 0.05)

df.test <- read.csv("02_Results/07_DETs_results_clean/Go_Enrich_ShortTermExclu_7.5C.csv") %>% 
  filter(padj <= 0.05)

df.test$Term[df.test$Term %>% as.vector %>% grep("stress", ., ignore.case = TRUE)]


# Expression level of HSP genes ------------------------------------------------

rm(list = ls())
library(tidyverse)

list.DETs <- read.csv("02_Results/07_DETs_results_clean/List_DETs.csv",
                      row.names = 1)

# dets.long <- list.DETs %>% filter(Compare == "LongTerm")
# dets.short.5 <- list.DETs %>% filter(Compare == "ShortTerm_5")
# dets.short.7.5 <- list.DETs %>% filter(Compare == "ShortTerm_7.5")

list.go.hsp <- c("GO:0070846",
                 "GO:1990565",
                 "GO:0051879",
                 "GO:0051008",
                 "GO:0030544",
                 "GO:0110078",
                 "GO:0001018",
                 "GO:0030192", 
                 "GO:0030191",
                 "GO:0008077",
                 "GO:0010538",
                 "GO:0010539",
                 "GO:0010545",
                 "GO:0010546" 
) 


## In annotation files 
annot.file <- read.csv("00_Data/03b_Annotations/annot_df.csv", header = T, row.names = 1)

hsp <- subset(annot.file, GO %in% list.go.hsp) 
## --> 125 transcripts associated to HSP


## subset DETs belonging to HSP family

dets.hsp <- list.DETs %>% filter(Transcrip_name %in% hsp$transcript_ID)

table(dets.hsp$Compare)

## 26 DETs, belonging to wich cluster/modules ?

# list of transcript LO2
gene.cluster <- read.csv("02_Results/07_DETs_results_clean/Longterm_Gene_cluster.csv",
                     row.names = 1) %>%
  filter(genes %in% dets.hsp$Transcrip_name) %>%
  select(genes, cluster_new) %>%
  distinct %>%
  left_join(., hsp, by = c("genes" = "transcript_ID"))

table(gene.cluster %>% select(cluster_new, genes) %>% distinct() %>% pull(cluster_new))

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
  mutate(DETs = ifelse(transcript_ID %in% dets.hsp$Transcrip_name, "Yes", "No")) %>%
  left_join(., gene.cluster %>% select(genes, cluster_new) %>% distinct,
            by = c("transcript_ID" = "genes"))


ggplot(hsp.vst %>% filter (Delta_short == 0) %>%
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
  facet_wrap(~ group_names, ncol = 4) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, vjust = 0)) +
  theme(plot.margin = margin(.5, .5, .5, .5, "cm")) 

# 


test <- hsp.vst %>% filter (Delta_short == 0) %>%
  mutate(group_names = recode(cluster_new, 
                              "Group L01 (1490)" = "Group L01 (12 DETs)",
                              "Group L02 (1228)" = "Group L02 (11 DETs)", 
                              "Group L05 (180)" = "Group L05 (2 DETs)",
                              .missing = "81 not-DETs"))
