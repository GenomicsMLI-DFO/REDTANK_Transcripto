# Info --------------------------------------------------------------------

# Trim reads for low quality and adaptators (Trim Galore!)
# 
# CL - May 2022

# Library -----------------------------------------------------------------

rm(list = ls())

library(here)
library(tidyverse)

# 1. Reads cleaning process (Trim Galore!)--------------------------------------

folder.in = file.path(here::here(), "00_Data", "00_RawData")
folder.out = file.path(here::here(), "00_Data", "01_TrimGalore")

R1.list <- list.files(folder.in, pattern = "R1.fastq.gz")
R2.list <- list.files(folder.in, pattern = "R2.fastq.gz")

paired.list <- lapply(1:length(R1.list), function(i){
  c(R1.list[i], R2.list[i])
})

# Add python env 
Sys.setenv(PATH = paste(c("/home/genobiwan/Documents/PythonVenv/GenoBaseEnv/bin",
                          Sys.getenv("PATH")),
                        collapse = .Platform$path.sep))

# Process Trim Galore!
## Auto-detecting adapter type

parallel::mclapply(paired.list,
                   FUN = function(i){
                     cmd <- paste("--phred33", 
                                  "--quality", 30, 
                                  "--length", 20,  # min length
                                  "--clip_R1", 12, # trim 12 bp at 5'
                                  "--clip_R2", 12, 
                                  "--output_dir", folder.out,
                                  "--paired", 
                                  paste(folder.in, "/", i[1], sep = ""),
                                  paste(folder.in, "/", i[2], sep = ""),
                                  sep = " ")
                     
                     system2("trim_galore", cmd)
                     
                   } ,
                   mc.cores = 20)


## Read ----------------------------------------------------------------------

# Define new folders
gc()

rm(folder.in)
rm(folder.out)

folder.in = file.path(here::here(), "00_Data", "01_TrimGalore")
folder.out = file.path(here::here(), "02_Results", "01_FastQC", "02_TrimGalore")

# How many files?
Nfiles <- length(list.files(folder.in, pattern = ".fq.gz"))
cat("Performing a FastQC analysis on", Nfiles,"files \n")

# Process fastqc
parallel::mclapply(list.files(folder.in, full.name = T, pattern = ".fq.gz"),
                   FUN = function(i){cmd <- paste("--outdir", folder.out, i, "-q")
                   system2("fastqc", cmd)
                   } ,
                   mc.cores = 20)

## Process multiQC

for (s in c("R1", "R2")) {
  cmd <- paste(list.files(folder.out, full.names = T) %>%
                 str_subset(paste0("_", s)) %>%
                 str_subset(".zip"),
               "--outdir", file.path(folder.out, "MultiQC_report"),
               "--filename", paste0("multiqc_report_", s, ".html"),
               "-f" # to rewrite on previous data
  )
  
  system2("multiqc", cmd)
}
