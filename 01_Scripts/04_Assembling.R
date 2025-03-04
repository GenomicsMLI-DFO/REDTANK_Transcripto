# Info --------------------------------------------------------------------

# Assembling RNA-Seq alignments into potential transcripts
# (de novo ref transcritpome)
# Structural annotation of the genome
# 
# CL - May 2022

# Library -----------------------------------------------------------------

rm(list = ls())
gc()
library(here)
library(tidyverse)

library(Biostrings)
library(polyester)
library(ballgown)


# Potential transcript for each sample ----

folder.in <- file.path(here(), "00_Data", "02_Hisat2")

files.in <- paste(folder.in, 
                 list.files(folder.in, pattern = "hisat2.bam$"), sep="/")

NumCores = 20

parallel::mclapply(files.in,
                   FUN = function(i){
                     
                     file.out <- i %>%
                       str_replace("02_Hisat2", "03a_StringTie") %>%
                       str_replace("_hisat2.bam", ".gtf") 
                     
                     cmd <- paste("-o",
                                  file.out,
                                  "--conservative", #Assembles transcripts in a conservative mode. Same as -t -c 1.5 -f 0.05
                                  "-p", 2, # number of processing threads (CPUs)
                                  i)
                     
                     system2("stringtie", cmd)
                   },
                   mc.cores = NumCores
)


#  Assemble transcripts into a non-redundant set of transcripts ----

folder.out <- file.path(here(), "00_Data", "03b_Annotations")
files.in <- paste(file.path(here(), "00_Data", "03a_StringTie"),
                  list.files(file.path(here(), "00_Data", "03a_StringTie"), pattern = ".gtf"),
                  sep = "/")

cmd.merge <- paste("-p", 20, 
                   "-l", "SEFA", # Prefix for the name of the output transcript: Sebastes faciatus
                   "--merge",
                   "-m", 50, # min transcript length
                   "-c", 0, # min coverage
                   "-F", 1, # min FPKM
                   "-T", 1, # min TPM
                   "-f", 0.01, # isoform fraction
                   "-g", 250, # Minimum locus gap separation value
                   "-o", 
                   paste(folder.out, "Sfaciatus_GeneAnnotation.gtf", sep = "/"), 
                   files.in) 

system2("stringtie", cmd.merge)

# Extract sequence from GTF > fasta
#BiocManager::install("ballgown")

gtf <-  gffRead(file.path(here(),
                          "00_Data", 
                          "03b_Annotations", 
                          'Sfaciatus_GeneAnnotation.gtf'))

seqs <- readDNAStringSet(file.path(here(),
                                   "00_Data", 
                                   "00_Genome_Ref", 
                                   "ADN_SF_S_21_04830.asm.bp.p_ctg.fa"))

tempo <- names(seqs)
tempo2 <- strsplit(tempo, " ")

names <- mat.or.vec(0,0)

for (i in 1:length(tempo2)){
  liste.tempo <- tempo2[[i]]
  name.tempo <- liste.tempo[1]
  names <- c(names, name.tempo)
}

names(seqs) <- names

test1 <- seq_gtf(gtf, seqs, exononly = TRUE, idfield = "transcript_id", 
                 attrsep = "; ")

writeXStringSet(test1, 
                file.path(here(), 
                          "00_Data", 
                          "03b_Annotations",
                          "Sfaciatus_transcript_seq.fasta.gz"), 
                append=FALSE,
                compress=T, compression_level=NA, format="fasta")


# Unzip fasta files

system2("gunzip", 
        file.path(here(), 
                  "00_Data", 
                  "03b_Annotations",
                  "Sfaciatus_transcript_seq.fasta.gz"))

# Gene annotation statistics ----

featureTable <- read.delim(paste(file.path(here(), "00_Data", "03b_Annotations"),
                                 "/Sfaciatus_GeneAnnotation.gtf", sep = ""),
                           comment.char = "#", header = F, row.names = NULL)
names(featureTable) <- c("seqname", "source", "feature", "start", "end",
                          "score", "strand", "frame", "attribute")

featureTable <- separate(featureTable, "attribute", into = c("gene_id",
                                                             "transcript_id",
                                                             "exon_number", NA), 
                         sep = ";")

featureTable$gene_id <- gsub("gene_id ", "", featureTable$gene_id)
featureTable$transcript_id <- gsub("transcript_id ", "", featureTable$transcript_id)
featureTable$exon_number <- gsub("exon_number ", "", featureTable$exon_number)

# number of genes
length(unique(featureTable$gene_id))

# number of transcript
length(unique(featureTable$transcript_id))

# number of transcript per genes
gene.transcript <- featureTable %>% 
  group_by(gene_id) %>% 
  summarise(nb_transcript = length(unique(transcript_id)))
summary(gene.transcript$nb_transcript)

# number of exon per genes
gene.exon <- featureTable[!is.na(featureTable$exon_number),] %>% 
  group_by(gene_id) %>% 
  summarise(nb_exon = as.numeric(max(exon_number)))
            
summary(gene.exon$nb_exon)

table(featureTable$strand, featureTable$feature)

featureTable$length <- featureTable$end - featureTable$start+1

## Select subset of features having "transcript" as "feature" 
transcript <- subset(featureTable, feature == "transcript")
n50.transcript <- Biostrings::N50(transcript$length)
summary(transcript$length)

Fig.transcript <- ggplot(transcript, aes(x = length)) +
  geom_histogram(bins = 100, fill = "lightgrey", col = "black") +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) {10^x}),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_bw() +
  geom_vline(xintercept = n50.transcript, col = "red", linetype = "dashed") +
  annotate(geom = "text",
           label = "N50 = 35,346 bp",
           x = n50.transcript + 5000,
           y = 1000, hjust = 0) +
  ylab("Transcript count") +
  xlab("Transcript length (bp, log10 scale)")
  
Fig.transcript

ggsave(file.path(here::here(), "02_Results", "denovo_transcripts_length.png"),
       Fig.transcript,
       width = 4, height = 3, scale = 1.5)

# genes lengths including only the exons or with the introns
gene.exon <- featureTable %>% 
  group_by(gene_id) %>% 
  subset(., !exon_number == " ") %>%
  summarise(size.gene.exon = sum(length),
            size.gene.all = max(end) - min(start) +1)

Biostrings::N50(gene.exon$size.gene.exon)
Biostrings::N50(gene.exon$size.gene.all)

