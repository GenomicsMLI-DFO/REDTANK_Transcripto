# Info --------------------------------------------------------------------

# Functional annotation of de novo transcritpome
# 
# CL - May 2022

# Library -----------------------------------------------------------------

rm(list = ls())
gc()
library(here)
library(tidyverse)

# Blastx nucl to proteom --------------------------------------------------

# add Blast to PATH
old_path <- Sys.getenv("PATH")
blast_path <- "/home/genobiwan/Documents/Programs/ncbi-blast-2.12.0+-src/c++/ReleaseMT/bin/"

Sys.setenv(PATH = paste(blast_path, old_path, sep = ":"))

# Make uniprot_sprot database
cmd <- paste("-in", file.path(here(), "00_Data", "00_Uniprot_db", "uniprot_sprot.fasta"),
             "-out", file.path(here(), "00_Data", "00_Uniprot_db", "Uniprot_db", "uniprot_sprot"),
             "-dbtype", "prot", 
             "-parse_seqids")

system2("makeblastdb", cmd)

# Blasx against uniprot_sprot

fasta.in <- file.path(here(), 
                      "00_Data", 
                      "03b_Annotations",
                      "Sfaciatus_transcript_seq.fasta")

file.out <- file.path(here(), 
                      "00_Data", 
                      "03b_Annotations",
                      "Sfaciatus_transcript_seq.swissprot")

uniprot_sprot <- file.path(here(), "00_Data", "00_Uniprot_db", "Uniprot_db", "uniprot_sprot")
  
cmd <- paste("-query", fasta.in,
             "-db", uniprot_sprot,
             "-num_threads", 20,
             "-mt_mode", 1, # multi-thread mode, split by queries
             "-evalue", 1e-3,
             "-max_target_seqs", 1,
             "-outfmt", 6,
             "-out", file.out)

system2("blastx", cmd)

# Get annotation information from uniprot --------------------------------------

# Global variables
swissprot_hits <- read.table(file.path(here(), "00_Data", "03b_Annotations", "Sfaciatus_transcript_seq.swissprot"),
                             header = F)
annot_folder <- file.path(here(), "00_Data", "00_Uniprot_db", "uniprot_feature_temp")

# Download info from uniprot for each hit 
df.hits <- data.frame(feature = swissprot_hits$V1,
                      hit = swissprot_hits$V2)

list.wget <- paste("http://www.uniprot.org/uniprot/", 
                   df.hits$hit,
                   ".txt", sep = "")
list.out <- paste(annot_folder, "/",
                  df.hits$feature,
                  ".info", sep = "")

list.cmd <- lapply(1:length(list.wget), function(i){
  c(list.wget[i], list.out[i])
})

parallel::mclapply(list.cmd,
                   FUN = function(i){
            
                     cmd <- paste("-q", 
                                  "-O",
                                  "-",
                                  i[1],
                                  ">",
                                  i[2])
                     
                     system2("wget", cmd)
                   } ,
                   mc.cores = 20)

# Retrieve GO terms from downloaded info (copy the following line in terminal)

##  for i in /media/genobiwan/Extra_Storage/Projets/REDTANK_Transcripto/00_Data/00_Uniprot_db/uniprot_feature_temp/*.info; do echo -e "$(basename $i | perl -pe 's/\.info//')\t$(cat $i | grep -E '^DR\s+GO' | awk '{print $3}' | perl -pe 's/\n//')"; done > /media/genobiwan/Extra_Storage/Projets/REDTANK_Transcripto/00_Data/03b_Annotations/all_go_annotations.csv

# # Annotation file
# 
# path.in <- file.path(here(), "00_Data", "00_Uniprot_db", "uniprot_feature_temp")
# path.out <- file.path(here(), "00_Data", "03b_Annotations")
# 
# list.info <- list.files(path.in)
# 
# parallel::mclapply(list.info,
#                    FUN = function(i){
# 
#                      transcript <- i %>% str_replace(., ".info", "")
#                      gene <- i %>% strsplit(., "[.]") %>% unlist(.) %>% .[2] %>% paste("SEFA", ., sep = ".")
#                      
#                      
#                    } ,
#                    mc.cores = 20)

