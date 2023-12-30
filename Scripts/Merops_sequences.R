library(data.table)

sequence1 <- fread("/home/svenja/workspace/Merops/meropsweb121/sequence_part1.csv")
sequence2 <- fread("/home/svenja/workspace/Merops/meropsweb121/sequence_part2.csv")
sequence3 <- fread("/home/svenja/workspace/Merops/meropsweb121/sequence_part3.csv")
sequence4 <- fread("/home/svenja/workspace/Merops/meropsweb121/sequence_part4.csv")

sequence <- rbind(sequence1, sequence2, sequence3,sequence4)

fwrite(sequence, file = "/home/svenja/workspace/Merops/meropsweb121/sequence_all.csv")

organism <- fread("/home/svenja/workspace/Merops/meropsweb121/organism.csv")

organism_bac <- organism[kingdom == "Bacteria", ]
vec_bac <- organism_bac$merops_taxonomy_id
sequence <- as.data.frame(sequence)
sequence_bac <- sequence[sequence$merops_taxonomy_id %in% vec_bac, ]

fwrite(sequence_bac, file = "/home/svenja/workspace/Merops/meropsweb121/sequence_bac.csv")


library(Biostrings)
g <- readDNAStringSet("/Users/svenjabusche/Downloads/HRGM_Genome_4163.fna")
g

amino <- translate(g, genetic.code=GENETIC_CODE, no.init.codon=FALSE,
                   if.fuzzy.codon="error")

out24a <- "~/Desktop/out24a.fasta"

writeXStringSet(amino,out24a, format = "fasta")

###convert csv to fasta files
gr <- GRanges(sequence_bac$merops_taxonomy_id, IRanges(sequence_bac$sequence))





######       Biostrings      #######
library(Biostrings)
seq = sequence_bac$seq
names(seq) = sequence_bac$id
dna = DNAStringSet(seq)
writeXStringSet(dna, "foo.fasta")

dna

head(sequence_bac)



#######                           Blastp                             ##########

listF<-list.files("Trich_prot_fasta/")

fa<-paste0("Trich_prot_fasta/",listF[i])

makeblastdb(fa, dbtype = "prot", args="")

bl <- blast("Trich_prot_fasta/Tri5640_1_GeneModels_FilteredModels1_aa.fasta", type="blastp")

seq <- readAAStringSet("NDRkinase/testSeq.txt")

cl <- predict(bl, seq)
