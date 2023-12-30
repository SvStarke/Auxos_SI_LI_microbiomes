library(MicrobiomeGS2)
library(ggplot2)
library(parallel)
library(ggpubr)
library(rstatix)
library(tidyverse)

n.cores <- detectCores()
cl <- makeCluster(min(n.cores-1, 3))
cl_big <- makeCluster(n.cores - 1)

min_rel_abun <- 0.001 # minimum relative abundance of species to be included in commFBA

mic <- new("Microbiome",
           uniq.table.file = "Data/asv_tab.tsv",
           model.mapping.file = "Output/asvs_to_HRGM.m8",
           sample.description.file = "Data/Metadata.csv",
           uniq.table.format = "R_table"
)

mic <- filter_mapping(mic, method.resolve.multiple = "first")
mic <- create_model_table(mic)
mic <- filter_samples(mic, min.seqs = 1000, max.unclassified = 0.3)

mic@model.table

# featch relevant models
models <- fetch_model_collection("/mnt/nuuk/2022/HRGM/models/",
                                 IDs = rownames(mic@model.table[-1,]))

models <- parLapply(cl_big, models,
                    function(mod) {
                      require(sybil)
                      
                      der <- deadEndMetabolites(mod)
                      mod <- rmReact(mod, react = der$der)
                      
                      return(mod)
                    })
stopCluster(cl_big)

clusterExport(cl, c("models","mic","min_rel_abun")) 

cFBA_out <- parLapply(cl, 1:ncol(mic@model.table),
                      fun = function(i) {
                        require(MicrobiomeGS2)
                        
                        n <- ncol(mic@model.table)
                        
                        rel_models <- rownames(mic@model.table[which(mic@model.table[,i] >0),])
                        rel_models <- rel_models[rel_models != "_unclassified"]
                        rel_abun <- mic@model.table[rel_models, i]
                        rel_abun <- rel_abun/sum(rel_abun)
                        rel_abun <- rel_abun[rel_abun >= min_rel_abun]
                        
                        return(rel_abun)
                      })

stopCluster(cl)
names(cFBA_out) <- colnames(mic@model.table)

##create a dataframe
test <- plyr::ldply(cFBA_out, rbind)
test <- as.data.table(test)
rownames(test) <- test$.id

test2 <- melt(test)

test2 <- test2[value >0,]
colnames(test2) <- c("Sample","Genomes","Freq")
#data <- aggregate(test2$value, by= list(test2$.id), FUN= function(x) {NROW(x)})
#data


metadata <- fread("Data/Metadata.csv")

##create file 
reimagine <- merge(metadata, test2, by.x="sample", by.y="Sample")

##analazying the data
# ##number of genomes per sample
# data <- aggregate(reimagine$Freq, by= list(reimagine$sample, reimagine$location), FUN= function(x) {NROW(x)})
# data                  

##mean numer of genomes per location
#mean_genomes <- aggregate(data$x, by=list(data$Group.2), FUN = mean)
