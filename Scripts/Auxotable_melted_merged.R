##################    create a dataframe with genome infos #####################
#table with auxotrophies first melted and then merged to all auxotrophies in one column


Auxotrophy_2 <- melt(Auxotrophy, id.vars = "Genomes",
                     value.name = "Prototrophy", variable.name = "Compound")

################only for HRGM###############################################
Metadata <- fread("/mnt/nuuk/2021/HRGM/REPR_Genomes_metadata.tsv")
if(length(models) == 3687)
  print(Auxotrophy_2 <- merge(Auxotrophy_2, Metadata, by.x = "Genomes", by.y = "HRGM name"))
if(length(models) == 3687)
  print(Auxotrophy_2[, phylum := str_match(`GTDB Taxonomy`, "p__.*;c__")[,1]])
if(length(models) == 3687)
  print(Auxotrophy_2[, phylum := gsub("p__|;c__","", phylum)])
if(length(models) == 3687)
  print(Auxotrophy_2[, phylum := gsub("_C$","", phylum)])
