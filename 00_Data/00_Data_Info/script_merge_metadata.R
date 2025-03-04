# Info -------------------------------------------------------------------------
#
# Script to merge the different metadata files 
# CL - 18/05/2022
#

library(dplyr)

# Files ------------------------------------------------------------------------

lab <- read.csv("00_Data/00_Data_Info/metaData_lab.csv")
lab2 <- read.csv("00_Data/00_Data_Info/metaData_lab2.csv")
extraction <- read.csv("00_Data/00_Data_Info/metaData_extraction.csv")
gq <- read.csv("00_Data/00_Data_Info/metaData_GQ.csv")

df.extract <- left_join(lab2, extraction)
df.lab <- left_join(df.extract, lab, by = c("Numero_reception_specimen" = "ID_lab"))

df.seq <- left_join(gq, df.lab, by = c("ID_Sample" = "Numero_unique_extrait"))

metaData <- data.frame(ID_sample = df.seq$ID_Sample,
                       ID_lab = df.seq$Numero_reception_specimen,
                       Lib = df.seq$Lib,
                       Seq_batch = df.seq$Seq_batch,
                       Numero_unique_groupe = df.seq$Numero_unique_groupe,
                       Numero_unique_specimen = df.seq$Numero_unique_specimen, 
                       Weight = df.seq$weight,
                       Length = df.seq$length_f,
                       TAG = df.seq$tag,
                       Tank_exp = df.seq$tank.exp,
                       Tank_sample = df.seq$tank.sample,
                       Temperature_acclim = df.seq$Temp.acclim,
                       Temperature_sampling = df.seq$Temp.sample,
                       Delta_acclim = df.seq$Delta.acclim,
                       Delta_temperature = df.seq$Temp.acclim - df.seq$Temp.sample,
                       DatHr_euthan = df.seq$DatHr_euthan, 
                       deltaT.ARN = df.seq$deltaT.ARN,
                       Time = df.seq$Time)

write.csv(metaData, "00_Data/00_Data_Info/metaData_all.csv", row.names = F)
