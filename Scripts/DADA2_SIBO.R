library(dada2); packageVersion("dada2")

path <- "/home/svenja/sratoolkit.3.0.0-ubuntu64/fastqSIBO/fastq_trim"
list.files(path)

### sort the data

fnFs <- sort(list.files(path, pattern="_1.trimmed.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.trimmed.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

### read quality

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

### filter and trim

filtFs <- file.path(path, "filtered", paste0(sample.names, "_1.trimmed.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_2.trimmed.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,230),
maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
compress=TRUE, multithread=TRUE)
head(out)

### learn the error rates

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

### dereplication

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

### sample inference

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]

### merge paired reads

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])

### construct sequence table

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

### remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

###
##Assign numbers to ASV
ASVs <- paste0("asv",1:ncol(seqtab.nochim))
asv.tab <- seqtab.nochim
colnames(asv.tab) <- ASVs

##track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
#names(asv.tab)[names(asv.tab) == "rn"] <- "Samples"
#View(asv_tab)
#asv.tab1 <- as.matrix(asv.tab[,-1])
#View(asv.tab1)
#rownames(asv.tab1) <- asv.tab$rn
#asv.tab1 <- as.table(t(asv.tab1))



# write ASV representative sequences (temporary!)
asv_tmp_fna <- tempfile()
asv_fasta <- DNAStringSet(colnames(seqtab.nochim))
names(asv_fasta) <- ASVs
writeXStringSet(asv_fasta, filepath = asv_tmp_fna)

# assign taxonomy
asv.tax <- assignTaxonomy(seqtab.nochim, "/home/svenja/sratoolkit.3.0.0-ubuntu64/Fastq-files/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
#View(asv.tax)
#taxa <- data.table(taxa)
#summary(taxa$Genus)
#table <- count(taxa, "Genus")
#sort(table$freq, decreasing = TRUE)

#all(colnames(seqtab.nochim) == rownames(asv.tax))
rownames(asv.tax) <- ASVs


###remove non-bacterial taxonomic classifications
#asv_tax <- data.table(asv_tax, keep.rownames = TRUE)
#asv_tax <- asv_tax[Kingdom != "Archaea"]
#asv_tax <- asv_tax[Kingdom != "Eukaryota"]
#names(asv_tax)[names(asv_tax) == "rn"] <- "ASV"
##ind_unknown <- which(is.na(asv_tax) | grepl("^Unknown$|^unknown$", asv_tax),arr.ind = T)
#asv_tax[ind_unknown] <- "Unknown"

###NOCH ANPASSEN!!!
# Removing non-bacterial ASVs
cat("Removing non-bacterial ASVs ...\n")
ind_NA <- which(is.na(asv.tax[,"Kingdom"]))
ind_noBac <- which(asv.tax[,"Kingdom"] != "Bacteria" | asv.tax[,"Order"] == "Chloroplast" | asv.tax[,"Family"] == "Mitochondria")
ind_noBac <- c(ind_NA, ind_noBac)
if(length(ind_noBac) > 0) {
  asv.tax <- asv.tax[-ind_noBac,]
  asv.tab <- asv.tab[,-ind_noBac]
  asv_fasta <- asv_fasta[-ind_noBac]
}
cat("\tremoved",length(ind_noBac),"ASVs\n")
cat("\tremaining ASVs:",length(asv_fasta),"\n")

# Filter samples
cat("Filter samples by total counts ...\n")
asv.tab <- t(asv.tab)
View(asv.tab)

min_seq_count <- 500
min_seq_count_SIBO <- 200
ind_spl_rm <- which(colSums(asv.tab) < min_seq_count_SIBO)
if(length(ind_spl_rm) != 0) {
  cat("Removing", length(ind_spl_rm), "samples, which have less than",min_seq_count_SIBO, "reads:")
  cat("\n\t",paste(colnames(asv.tab)[ind_spl_rm], collapse = "\n\t"),"\n")
  asv.tab <- asv.tab[,-ind_spl_rm]
}
writeLines(names(ind_spl_rm), con = "Data/removed_samples_by_counts.txt")

# aggregating counts by taxonomy predictions
ind_unknown <- which(is.na(asv.tax) | grepl("^Unknown$|^unknown$", asv.tax),arr.ind = T)
asv.tax[ind_unknown] <- "Unknown"

tax.tab <- list()
tax.levels <- colnames(asv.tax)[-1]
for(taxlvl in 1:length(tax.levels)) {
  if(taxlvl == 1) {
    lvls_tmp <- asv.tax[,2]
  } else {
    lvls_tmp <- apply(asv.tax[,2:(taxlvl+1)],1,function(x) paste(x, collapse = ";"))
  }
  
  grps <- unique(lvls_tmp)
  mat_tmp <- matrix(0, ncol = ncol(asv.tab), nrow = length(grps))
  rownames(mat_tmp) <- grps
  colnames(mat_tmp) <- colnames(asv.tab)
  
  for(grp_i in grps) {
    if(sum(lvls_tmp == grp_i) > 1) {
      mat_tmp[grp_i,] <- colSums(asv.tab[lvls_tmp == grp_i,])
    } else {
      mat_tmp[grp_i,] <- asv.tab[lvls_tmp == grp_i,]
    }
  }
  tax.tab[[tax.levels[taxlvl]]] <- mat_tmp
}

###export files
# Export processed data
cat("Exporting processed data tables ...\n")
write.table(asv.tab, file = "Data/asv_tab_SIBO.tsv", sep = "\t", quote = F)
#write.table(asv.tab_rarefied, file = "data/dada/asv_tab_rarefied.tsv", sep = "\t", quote = F)
write.table(asv.tax, file = "Data/asv_tax_SIBO.tsv", sep = "\t", quote = F)
writeXStringSet(asv_fasta, filepath = "Data/asv_seqs_SIBO.fna")
#fwrite(alpha_div, file = "data/dada/DT_alpha_diversity.tsv", sep = "\t", quote = F)
for(i in names(tax.tab)) {
  write.table(tax.tab[[i]], file = paste0("Data/tax_tab_SIBO_",i,".tsv"), sep = "\t", quote = F)
}










### track reads through the pipeline

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(merges, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

### assign taxonomy

taxa <- assignTaxonomy(seqtab.nochim, "~/tax/silva_nr_v128_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "~/tax/silva_species_assignment_v128.fa.gz")
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

### evaluate accuracy

unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE)
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")

