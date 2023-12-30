library(MicrobiomeGS2)
library(stringr)
library(dplyr)
library(ggplot2)
library(data.table)
library(cplexAPI)


###HRGM
Metadata <- fread("Data/REPR_Genomes_metadata.tsv")
relGenomes <- Metadata[`Completeness (%)`>= 85 & `Contamination (%)` <=2 & !grepl("^d__Archaea", `GTDB Taxonomy`), `HRGM name`]

models <- fetch_model_collection("/mnt/nuuk/2022/HRGM/models/", IDs = relGenomes)
