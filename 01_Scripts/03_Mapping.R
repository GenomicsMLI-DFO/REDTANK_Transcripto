# Info --------------------------------------------------------------------

# Mapping reads on reference genome, with spliced alignment (Hisat2)
# 
# CL - May 2022

# Library -----------------------------------------------------------------

rm(list = ls())
gc()
library(here)
library(tidyverse)

# add Hisat2 to PATH
old_path <- Sys.getenv("PATH")
hisat_path <- "/home/genobiwan/Documents/Programs/hisat2-2.2.1"

Sys.setenv(PATH = paste(hisat_path, old_path, sep = ":"))

# Add python env 
Sys.setenv(PATH = paste(c("/home/genobiwan/Documents/PythonVenv/GenoBaseEnv/bin",
                          Sys.getenv("PATH")),
                        collapse = .Platform$path.sep))

# Build index ----
ref_genome <- paste(file.path(here::here(), "00_Data", "00_Genome_Ref"),
                    "/ADN_SF_S_21_04830.asm.bp.p_ctg.fa", 
                    sep = "")

index_genome <- paste(file.path(here::here(), "00_Data", "00_Genome_Ref"),
                    "/SFaciatus_ref", 
                    sep = "")

cmd <- paste(ref_genome, 
             index_genome)
system2("hisat2-build", cmd)

# Map reads ----

folder.in = file.path(here::here(), "00_Data", "01_TrimGalore")
folder.out = file.path(here::here(), "00_Data", "02_Hisat2")

## Match R1 and R2 
R1.list <- paste(folder.in, list.files(folder.in, pattern = "1.fq.gz"), sep="/")
R2.list <- paste(folder.in, list.files(folder.in, pattern = "2.fq.gz"), sep="/")

paired.list <- lapply(1:length(R1.list), function(i){
  c(R1.list[i], R2.list[i])
})

## Set number of cores to be used
NumCores = 20 

## Path to ref genome index
ref_index <- file.path(here::here(), "00_Data", "00_Genome_Ref", "SFaciatus_ref")

## Hisat2
parallel::mclapply(paired.list,
                   FUN = function(i){
                     
                     name1 <- i[1] %>%
                       str_replace(folder.in, "") %>%
                       str_replace("_val_1.fq.gz", ".bam") %>%
                       substr(., 5, 9)
                     
                     name2 <- i[1] %>%
                       str_replace(folder.in, "") %>%
                       str_replace("R1_val_1.fq.gz", "hisat2") %>%
                       substr(., (nchar(.)-18), (nchar(.)))
                     
                     file.out.summary <- paste(folder.out, "/", name1, name2, "_mapping_summary.txt", sep = "")

                     cmd <- paste("-p", 2, # nb of CPUs for alignment
                                  "-x", ref_index,
                                  "-1", i[1],
                                  "-2", i[2], 
                                  "--summary-file", file.out.summary,
                                  "| samtools sort --no-PG -l 0 -T",
                                  paste(folder.out, "/", name1, name2, sep = ""),
                                  "-O bam", ## sort bam files
                                  "| samtools view --no-PG -O bam -@ 2 -o", # compress 
                                  paste(folder.out, "/", name1, name2, ".bam", sep = ""))
                     
                     system2("hisat2", cmd)
                   },
                   mc.cores = NumCores
                   )


## Index the bam files
gc()
folder = file.path(here::here(), "00_Data", "02_Hisat2")
file.in <- paste(folder, list.files(folder, pattern = ".bam"), sep="/")

parallel::mclapply(file.in,
                   FUN = function(i){
                     
                     file.out <- i %>%
                       str_replace(".bam", ".bam.bai") 
                     
                     cmd <- paste("index",
                                  i,
                                  file.out)
                     
                     system2("samtools", cmd)
                   },
                   mc.cores = NumCores
)


# Move mapping_summary.txt files into results folder

cmd <- paste(file.path(here::here(), "00_Data", "02_Hisat2"),
             "/*mapping_summary.txt",
             " ",
             file.path(here::here(), "02_Results", "02_Mapping_summary"),
             sep = "")
system2("mv", cmd)
