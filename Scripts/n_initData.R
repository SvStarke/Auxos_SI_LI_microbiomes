library(data.table)
library(ggplot2)
library(vegan)
library(MicrobiomeGS2)
library(ggpubr)
library(ggtext)
# library(ComplexHeatmap)
# library(circlize)

#------------------------------------------------------------------------------#
# Central parameters
#------------------------------------------------------------------------------#
param <- list(
  hrgm.minCompleteness = 85,
  hrgm.maxContamination = 2,
  mero.minBitscore = 100,
  mero.minSCov = 0.95,
  mero.families = c("A" = "Aspartic (A) Peptidases",
                    "C" = "Cysteine (C) Peptidases",
                    "M" = "Metallo (M) Peptidases",
                    "N" = "Asparagine (N) Peptide Lyases",
                    "S" = "Serine (S) Peptidases",
                    "T" = "Threonine (T) Peptidases",
                    "P" = "Mixed (P) Peptidases",
                    "U" = "Unknown (U) Catalytic Type"),
  mero.familyColors = c("A" = "#1b9e77", 
                        "C" = "#d95f02",
                        "M" = "#7570b3",
                        "N" = "#cccccc",
                        "S" = "#e7298a",
                        "T" = "#66a61e",
                        "P" = "#e6ab02",
                        "U" = "#232323"),
  sibo.colors = c("SIBO" = "#fc8d62",
                  "non-SIBO" = "#8da0cb"),
  git.colors = c("SI-Duodenum" = "#fdd49e",
                 "SI-Jejunum" = "#fc8d59",
                 "SI-FD" = "#d7301f",
                 "LI (Stool)" = "#7f0000"),
  modelPath = "/mnt/nuuk/2021/HRGM/models_20230430/",
  essentialAA = c("Val","Met","Leu","Ile","Trp","Phe","Lys","His","Thr")
)

#------------------------------------------------------------------------------#
# Meta data HRGM genomes
#------------------------------------------------------------------------------# 
db_hrgm <- fread("Data/REPR_Genomes_metadata.tsv")
db_hrgm <- db_hrgm[`Completeness (%)` >= param$hrgm.minCompleteness & `Contamination (%)` <= param$hrgm.maxContamination & !grepl("^d__Archaea", `GTDB Taxonomy`)]
db_hrgm <- db_hrgm[!grepl("^NISW", `HRGM name`)] # removes the non-HRGM models

#------------------------------------------------------------------------------#
# Load gapseq models and predict auxotrophies
#------------------------------------------------------------------------------#
models <- fetch_model_collection(param$modelPath, IDs = db_hrgm$`HRGM name`)
auxos <- predict_auxotrophies(models, min.growth = 1e-7,
                              min.growth.fraction = 1e-7)
auxos <- do.call("rbind", auxos)
auxos[which(is.na(auxos), arr.ind = TRUE)] <- 1
auxos <- abs(auxos - 1) # now: 1 = auxotrophy

#------------------------------------------------------------------------------#
# Load and filter peptidase predictions including SignalP predictions
#------------------------------------------------------------------------------# 
mero <- fread("Data/merops_HRGM.m8.gz")
colnames(mero) <- c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
                       "evalue", "bitscore")
mero[, model := substr(qaccver, 1, 16)]
mero <- mero[grepl("^HRGM",qaccver)]

# filter to translocated peptidases only (signalP-Predictions)
sgp <- fread("Data/signalP_predictions.tsv.gz")
mero <- mero[qaccver %in% sgp$gene]

# filter by subject (Merops seqs) coverage
mero <- merge(mero, fread("Data/merops_scan_filt_length.csv"), by.x = "saccver", by.y = "mernum")
mero[, scov := (send - sstart)/len]
mero <- mero[scov >= param$mero.minSCov]
mero <- mero[bitscore >= param$mero.minBitscore]

mero_fam <- fread("Data/domain2.csv.gz")
mero <- merge(mero, mero_fam[,.(mernum, code, protein)], by.x = "saccver", by.y="mernum")
rm(mero_fam)
mero[,family := substr(code, 1, 1)]
mero$tmp <- 1

# remove duplicate hits of single hrgm genes to multiple merops IDs
mero <- mero[order(qaccver, -bitscore)]
mero <- mero[!duplicated(qaccver)]

mero_mat <- dcast(mero,  code ~ model, fun.aggregate = sum, value.var = "tmp")
mero[, tmp := NULL]
tmp <- mero_mat$code 
mero_mat <- as.matrix(mero_mat[,-1])
rownames(mero_mat) <- tmp
rm(tmp)

hrgmNoPepd <- db_hrgm$`HRGM name`[!(db_hrgm$`HRGM name` %in% colnames(mero_mat))]
tmpmat <- matrix(0, ncol = length(hrgmNoPepd), nrow = nrow(mero_mat),
                 dimnames = list(rownames(mero_mat),hrgmNoPepd))
mero_mat <- cbind(mero_mat, tmpmat)
mero_mat <- mero_mat[,db_hrgm$`HRGM name`]
mero_mat_bin <- mero_mat > 0
rm(tmpmat)

#------------------------------------------------------------------------------#
# Meta data cohorts
#------------------------------------------------------------------------------#

# Reimagine
db_meta_reim <- fread("Data/Metadata.csv")
db_meta_reim <- db_meta_reim[,.(sample, location.main = Isolation_Source,
                                location, PatNo, Group1, Group2, Group3)]
db_meta_reim[, location := ifelse(location == "Stool","LI (Stool)", paste0("SI-",location))]
db_meta_reim[location == "SI-Ileum", location := "SI-FD"]

# SIBO cobort
db_meta_sibo <- fread("Data/SIBO.csv")
db_meta_sibo <- db_meta_sibo[,.(sample, SIBO = V41)]
db_meta_sibo[SIBO == "SIBO-replicate", SIBO := "SIBO"]

#------------------------------------------------------------------------------#
# Relative abundance data
#------------------------------------------------------------------------------#

# Reimagine
mic <- new("Microbiome",
           uniq.table.file = "Data/asv_tab.tsv",
           model.mapping.file = "Output/asvs_to_HRGM.m8",
           sample.description.file = "Data/Metadata.csv",
           uniq.table.format = "R_table"
)
mic@model.mapping <- mic@model.mapping[!grepl("^NISW",target.label)] # removes the non-HRGM hits
mic <- filter_mapping(mic, method.resolve.multiple = "first")
mic <- create_model_table(mic)
mic <- filter_samples(mic, min.seqs = 1000, max.unclassified = 0.3)
abun_reim <- mic@model.table[-1,]
abun_reim <- abun_reim[rownames(abun_reim) %in% db_hrgm$`HRGM name`,]
abun_reim_rel <- t(t(abun_reim)/colSums(abun_reim))
rm(mic)

# SIBO cohort
mic <- new("Microbiome",
           uniq.table.file = "Data/asv_tab_SIBO.tsv",
           model.mapping.file = "Output/asvs_to_HRGM_SIBO.m8",
           sample.description.file = "Data/SIBO.csv",
           uniq.table.format = "R_table"
)
mic@model.mapping <- mic@model.mapping[!grepl("^NISW",target.label)] # removes the non-HRGM hits
mic <- filter_mapping(mic, method.resolve.multiple = "first")
mic <- create_model_table(mic)
mic <- filter_samples(mic, min.seqs = 1000, max.unclassified = 0.3)
abun_sibo <- mic@model.table[-1,]
abun_sibo <- abun_sibo[rownames(abun_sibo) %in% db_hrgm$`HRGM name`,]
abun_sibo_rel <- t(t(abun_sibo)/colSums(abun_sibo))
rm(mic)

#------------------------------------------------------------------------------#
# Helper functions
#------------------------------------------------------------------------------#
export_plot <- function(plot, path, width, height) {
  print(path)
  # Export PDF
  ggsave(paste0(path,".pdf"), plot = plot, width = width, height = height, device = cairo_pdf)
  
  # Export SVG
  ggsave(paste0(path,".svg"), plot = plot, width = width, height = height, device = svg)
  
  # Export PNG
  ggsave(paste0(path,".png"), plot = plot, width = width, height = height, device = png)
  
  return(TRUE)
}

pval_labs <- function(pval, ns.lab = "") {
  res <- ifelse(pval < 0.0001, "****",
                ifelse(pval < 0.001, "***",
                       ifelse(pval < 0.01, "**",
                              ifelse(pval < 0.05, "*",ns.lab))))
}
