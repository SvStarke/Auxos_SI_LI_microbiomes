library(MicrobiomeGS2)
library(ggplot2)
library(parallel)
library(ggpubr)
library(rstatix)

n.cores <- detectCores()
cl <- makeCluster(min(n.cores-1, 3))

min_rel_abun <- 0.001 # minimum relative abundance of species to be included in commFBA
# SIBO <- fread("Data/SIBO_2.csv")
# SIBO <- as.data.table(SIBO)
# colnames(SIBO)[1] <- "sample"
# write.csv(SIBO, file = "Data/SIBO.csv")

mic <- new("Microbiome",
           uniq.table.file = "Data/asv_tab_SIBO.tsv",
           model.mapping.file = "Output/asvs_to_HRGM_SIBO.m8",
           sample.description.file = "Data/SIBO.csv",
           uniq.table.format = "R_table"
)

mic <- filter_mapping(mic, method.resolve.multiple = "first")
mic <- create_model_table(mic)
mic <- filter_samples(mic, min.seqs = 1000, max.unclassified = 0.3)

mic@model.table

SIBO <- fread("Data/SIBO.csv")
SIBO

clusterExport(cl, c("mic")) 

cFBA_out <- parLapply(cl, 1:ncol(mic@model.table),
                      fun = function(i) {
                        require(MicrobiomeGS2)
                        
                        n <- ncol(mic@model.table)
                        
                        rel_models <- rownames(mic@model.table[which(mic@model.table[,i] >0),])
                        rel_models <- rel_models[rel_models != "_unclassified"]
                        rel_abun <- mic@model.table[rel_models, i]
                        rel_abun <- rel_abun/sum(rel_abun)
                        
                        return(rel_abun)
                      })

stopCluster(cl)
names(cFBA_out) <- colnames(mic@model.table)
warnings()
# ##create a dataframe
test <- plyr::ldply(cFBA_out, rbind)
test <- as.data.table(test)
rownames(test) <- test$.id
# 
test2 <- melt(test)
test2 <- test2[value >0,]
colnames(test2) <- c("Sample","Genomes","Freq")
#data <- aggregate(test2$value, by= list(test2$.id), FUN= function(x) {NROW(x)})
#data


metadata <- fread("Data/SIBO.csv")
colnames(metadata)[34] <- "Name2"
colnames(metadata)[12] <- "Name3"

##create file 
SIBO_study <- merge(metadata, test2, by.x="sample", by.y="Sample")

