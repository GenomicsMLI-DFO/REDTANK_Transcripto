# Info --------------------------------------------------------------------

# Raw sequence quality assessment
# 
# CL - May 2022

# Library -----------------------------------------------------------------

library(here)
library(tidyverse)

# 1. RAW quality assesment (fastqc) ---------------------------------------

folder.in = file.path(here::here(), "00_Data", "00_RawData")
folder.out = file.path(here::here(), "02_Results", "01_FastQC", "01_RawData")

## FastQC ----
Nfiles <- length(list.files(folder.in, pattern = ".fastq.gz"))
cat("Performing a FastQC analysis on", Nfiles,"files \n")

parallel::mclapply(list.files(folder.in, full.name = T, pattern = "fastq"),
                   FUN = function(i){cmd <- paste("--outdir", folder.out, i, "-q")
                   system2("fastqc", cmd)
                   } ,
                   mc.cores = 20)

## MultiQC ----

# Add python env 
Sys.setenv(PATH = paste(c("/home/genobiwan/Documents/PythonVenv/GenoBaseEnv/bin",
                          Sys.getenv("PATH")),
                        collapse = .Platform$path.sep))

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

