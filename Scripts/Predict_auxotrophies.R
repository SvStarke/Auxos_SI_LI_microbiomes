#Predict Auxotrophies
model.auxo <- predict_auxotrophies(models, min.growth = 1e-12, min.growth.fraction = 1e-12)


Auxotrophie <- data.frame(model.auxo)
# head(Auxotrophie) 
# summary(Auxotrophie)
# str(Auxotrophie)
# is.data.frame(Auxotrophie)

#column und rows tauschen
Auxotroph <- t(Auxotrophie)
Auxotroph[which(is.na(Auxotroph), arr.ind = T)] <- 1
#is.matrix(Auxotroph)
#data frame erzeugen
Auxotrophy <- data.frame(Auxotroph)
# is.data.frame(Auxotrophy)
# str(Auxotrophy)
# Auxotrophy
Genome <- rownames(Auxotrophy)
Auxotrophy$Genomes <- Genome
# ----
Auxotrophy <- as.data.table(Auxotrophy)

#filter genomes for completeness and contamination
Metadata <- fread("Data/REPR_Genomes_metadata.tsv")
relGenomes <- Metadata[`Completeness (%)`>= 85 & `Contamination (%)` <=2 & !grepl("^d__Archaea", `GTDB Taxonomy`), `HRGM name`]
Auxotrophy <- filter(Auxotrophy, Genomes %in% relGenomes)


