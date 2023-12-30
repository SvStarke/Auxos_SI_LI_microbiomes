#### Proteolytic activity SIBO vs non-SIBO   ####

source("Scripts/Mapping_ASV2HRGM_SIBO.R")

source("Scripts/Predict_auxotrophies.R")

source("Scripts/Auxotable_melted_merged.R")


#####  peptidases in SIBO vs non-SIBO

HRGM_m8 <- fread("Data/all.m8")
HRGM_m8


colnames(HRGM_m8) <- c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
                       "evalue", "bitscore")
HRGM_m8

## shorten length of HRGM_Genome name 
HRGM_m8$qaccver <- substr(HRGM_m8$qaccver, 1, 16)


genome_names_SIBO <- names(models)

##filter for peptidase sequences of genomes included in blastp result

HRGM_m8_SIBO <- HRGM_m8[qaccver %in% genome_names_SIBO]

SIBO_study


##peptidases found

domain <-fread("Data/domain.csv")
domain

domain_results <- merge(domain, HRGM_m8_SIBO, by.x="mernum", by.y="saccver")


##delete all inhibitors

domain_results <- domain_results[type != "inhibitor",]
domain_results$family <- substr(domain_results$code, 1,1)



### add Metainformation
SIBO_pep <- merge(SIBO_study, domain_results, by.x="Genomes", by.y="qaccver")
unique(SIBO_pep$code)
SIBO_pep$V41[SIBO_pep$V41 == "SIBO-replicate"] <- "SIBO"
SIBO_pep410 <- ggplot(SIBO_pep, aes(sample, code)) +
  geom_point(aes (colour = factor(family)), size =0.5) +
  facet_grid(~ V41, scales = "free", space = "free") +
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "white"),
        strip.text.x = element_text(colour = "black", size =11),
        strip.text.y= element_text(colour = "black", size =11),
        axis.title.x = element_text(size = 13, colour = "black"),
        axis.title.y = element_text(size = 13, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 13, colour = "black")) +
  labs(x="Samples", y="Peptidases", colour="Family") +
  scale_colour_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00"),
                      labels =  c("A" = "Aspartic peptidase", "C" = "Cysteine peptidase", "M" = "Metallopeptidase", "S" = "Serine peptidase",
                                  "T" = "Threonine peptidase", "U" = "Unknown catalytic type"),
                      breaks = c("A", "M", "C", "S", "T", "U")) +guides(colour = guide_legend(override.aes = list(size=3)))
SIBO_pep410
ggsave("Output/SIBO_peptidases.pdf", plot = SIBO_pep410,
       width = 12, height = 7)
# ggsave("Output/SIBO_peptidases.png", plot = SIBO_pep410,
#        width = 12, height = 7)

pep <- unique(SIBO_pep$code)
pep_list <- list()
k <- 1

for(pi in pep) {
  tmp_pep <- SIBO_pep[code == pi]
  vectortmp <- tmp_pep$V41
  vectortmp <- unique(vectortmp)
  if(length(vectortmp) == 2) {
    tmp_pep$V41 <- "Both"
  }
  pep_list[[k]] <- tmp_pep
  k <- k+1
}

pep_shared_SI_LI <- rbindlist(pep_list)

#
# ###visualization
# shared_pep <- ggplot(pep_shared_SI_LI, aes(V41, code)) +
#          geom_col(aes(colour = family)) +
#   xlab("") +
#   ylab("Peptidases") +
#   scale_colour_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7"),
#                       labels =  c("A" = "Aspartic peptidase", "C" = "Cysteine peptidase", "M" = "Metallopeptidase", "S" = "Serine peptidase",
#                                   "T" = "Threonine peptidase")) +
#   theme(axis.text.x = element_text(colour = "black"),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         panel.grid =  element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black")) +
#   guides(col = guide_legend(title = "Family",override.aes = list(shape = 15, size = 8, fill = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7")))) +
#   theme(legend.position = "none")
# 
# ggsave("Output/SIBO_pept_non_SIBO_both.pdf", plot = shared_pep,
#        width = 7, height = 5)


# pept <- unique(pep_shared_SI_LI$code)
# loc <- unique(pep_shared_SI_LI$V41)
# pep_list1 <- list()
# k <- 1
# 
# for(loci in loc) {
#   vector_pep <- as.numeric(length(pept))
#   tmp_loc <- pep_shared_SI_LI[V41 == loci]
#   tmp_loc_pep <- unique(tmp_loc$code)
#   length_tmp_pep_loc <- as.numeric(length(tmp_loc_pep))
#   perc <- (length_tmp_pep_loc/vector_pep)*100
#   loc_pep <- data.table(Location = loci,
#                         Perc = perc)
#   pep_list1[[k]] <- loc_pep
#   k <- k+1
#   
# }
# 
# perc_pep_loc <- rbindlist(pep_list1)
# perc_pep_loc
# perc_pep_loc$Nr <- 1
# 
# per_pep_loc <- ggplot(perc_pep_loc, aes(Nr,Perc)) +
#   geom_col(aes(fill = Location)) +
#   scale_fill_manual(values = c("#DDCC77", "#000000", "#888888"),
#                     labels =  c("Both" = "Both", "SIBO" = "SIBO", "non-SIBO" = "non-SIBO")) +
#   xlab("") +
#   ylab("Number of peptidases [%]") +
#   theme_bw() +
#   theme(axis.text.x = element_blank(),
#         panel.grid = element_blank(),
#         #panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         axis.ticks.x = element_blank())
# per_pep_loc
# 
# ggsave("Output/Percentage_Peptidases_SIBO.pdf", plot = per_pep_loc,
#        width = 4, height = 3)

###PCoA Plot
peptidases_bray_curtis1 <- SIBO_pep[,c("code", "Sample")]
t(peptidases_bray_curtis1)

pep <- unique(SIBO_pep$code)
library(dplyr)
library(vegan)
library(ape)
pep_list <- list()
k <- 1

for(pi in pep) {
  pep_tmp <- SIBO_pep[code == pi]
  pep_SI_LI <- pep_tmp %>% count(Sample, V41)
  pep_SI_LI$code <- pi
  pep_list[[k]] <- pep_SI_LI
  k <- k+1
}
pep_all_BC <- rbindlist(pep_list)
new_pep_all <- dcast(pep_all_BC, V41 + Sample ~ code, value.var = "n")
new_pep_all[is.na(new_pep_all)] <- 0
tmp2 <- new_pep_all[,-c(1,2)]
dist <- vegdist(tmp2, method = "bray")
testx <- pcoa(dist)
dist2 <- as.matrix(dist)

biplot(testx)
bray_curtis_pcoa_df <- data.frame(pcoa1 = testx$vectors[,1], 
                                  pcoa2 = testx$vectors[,2])
bray_curtis_plot <- ggplot(data = bray_curtis_pcoa_df, aes(x=pcoa1, y=pcoa2, colour = new_pep_all$V41)) +
  geom_point() +
  labs(x = "PCoA1",
       y = "PCoA2") +
  theme(title = element_text(size = 10)) +
  guides(colour = guide_legend(title = "Disease state")) +
  scale_color_manual(values = c("#661100", "#999933")) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  stat_ellipse() +
  theme(axis.title.y = element_text(colour = "black", size = 13, margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 13, margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=10, colour = "black", vjust = 1,  angle = 0, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size =10, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) 

bray_curtis_plot

ggsave("Output/PcoA_plot_SIBO.pdf", plot = bray_curtis_plot,
       width = 4, height = 3)

#statistical evaluation PERMANOVA
library(vegan)
adonis2(dist2 ~ V41, data = new_pep_all)


# ### 10 most abundant peptidases in SIBO vs non-SIBO
# SIBO_state <- unique(SIBO_pep$V41)
# 
# list <- list()
# k <- 1
# 
# for(si in SIBO_state) {
#   print(si)
#   x1 <- SIBO_pep[V41 == si]
#   sam <- unique(x1$sample)
#   for(sami in sam) {
#     x2 <- x1[sample == sami]
#     x3 <- data.table(abunpepti = names(sort(table(x2$code),decreasing=TRUE)[1:10]),
#                      Disease = si,
#                      sample = sami)
#     list[[k]] <- x3
#     k <- k+1
#   }
# }
# 
# pep_10_SIBO <- rbindlist(list)
# pep_10_SIBO
# 
# disease <- unique(pep_10_SIBO$Disease)
# list <- list()
# k <- 1
# 
# for(si in disease) {
#   x5 <- pep_10_SIBO[Disease == si]
#   x6 <- data.table(abunpepti = names(sort(table(x5$abunpepti),decreasing=TRUE)[1:10]),
#                    Disease = si)
#   list[[k]] <- x6
#   k<- k+1
# }
# 
# pep_10_SIBO <- rbindlist(list)
# pep_10_SIBO
# 
# 
# top10 <- ggplot(pep_10_SIBO, aes(Disease, abunpepti)) +
#   geom_tile(aes(fill = Disease), colour = "white", lwd = 1.5, linetype = 1) +
#   theme(axis.text.x = element_text(angle = 90, colour = "black", vjust = 0.5),
#         axis.text.y = element_text(colour = "black")) +
#   coord_equal() +
#   theme(panel.background = element_blank(),
#         axis.line.x = element_line(colour = "black"),
#         axis.line.y = element_line(colour = "black"),
#         axis.title.y =element_blank(),
#         axis.title.x = element_blank()) +
#   scale_fill_manual(breaks = levels(pep_10_SIBO$Disease),
#                     values = c("#000000", "#000000")) +
#   facet_grid("10 most found peptidases" ~., space="free") +
#   theme(strip.text.x = element_text(size = 11, colour = "black"))
# top10
# 
# ggsave("Output/SIBO_top10_peptidases.pdf", plot = top10,
#        width = 12, height = 7)
# 

# ####  abundance of peptidases ####
# ###boxplot
# 
# 
# domain_results_filt <- domain_results[,c("mernum","qaccver")]
# library(tidyverse)
# sum_pep <- domain_results_filt %>%
#   # count number of rows for each combination of server_id and protocol
#   group_by(mernum, qaccver) %>%
#   tally() %>%
#   # pivot the protocol names over the columns
#   pivot_wider(names_from=mernum, values_from=n)
# 
# sum_pep <- as.data.frame(sum_pep)
# sum_pep[is.na(sum_pep)] <- 0
# rownames(sum_pep) <- sum_pep$qaccver
# sum_pep <- sum_pep[,-1]
# sum_pep$count_pep <- rowSums(sum_pep)
# sum_pep$Genomes <- rownames(sum_pep)
# 
# sum_pep_filt <- sum_pep[,c(120,121)]
# 
# SIBO_sum_pep_filt <- merge(SIBO_study, sum_pep_filt, by.x = "Genomes", by.y="Genomes")
# #View(SIBO_sum_pep_filt)
# 
# SIBO_sum_pep_filt_2 <- SIBO_sum_pep_filt[ ,sum(count_pep*Freq), by = sample]
# colnames(SIBO_sum_pep_filt_2) <- c("sample", "Freq_pep")
# colnames(SIBO)[34] <- "Name2"
# colnames(SIBO)[12] <- "Name3"
# SIBO_study_sum_pep_sample <- merge(SIBO_sum_pep_filt_2, SIBO, by.x="sample", by.y="sample")
# 
# ggplot(SIBO_study_sum_pep_sample, aes(V41, Freq_pep))+
#   geom_bar(stat = "identity")
# 
# wilcox.test(Freq_pep ~ V41, data = SIBO_study_sum_pep_sample)

# SIBO_study_sum_pep_sample$V41[SIBO_study_sum_pep_sample$V41 == "SIBO-replicate"] <- "SIBO"
# 
# SIBO_pep_freq <- ggplot(SIBO_study_sum_pep_sample, aes(V41,Freq_pep, fill = V41))+
#   geom_boxplot(outlier.shape = NA) +
#   stat_compare_means(method = "wilcoxon", label.x = 0.5) +
#   xlab("Disease state") +
#   ylab("Abundance-weighted average of peptidase") +
#   scale_fill_manual(breaks = SIBO_study_sum_pep_sample$V41,
#                     values = c("#D55E00", "#009E73")) +
#   theme_bw() +
#   theme(axis.text.x = element_text(colour = "black"),
#         axis.text.y = element_text(colour = "black"),
#         axis.title.x = element_text(colour = "black", size = 11,  margin = margin(t = 10, r = 0, b = 0, l = 0)),
#         axis.title.y = element_text(colour = "black", size = 11, margin = margin(t = 0, r = 10, b = 0, l = 0))) +
#   theme(legend.position = "none")
# 
# SIBO_pep_freq
# 
# ggsave("Output/SIBO_freq_peptidases_Wilcox.pdf", plot = SIBO_pep_freq,
#        width = 6, height = 5)
# ggsave("Output/SIBO_freq_peptidases_t.test.pdf", plot = SIBO_pep_freq,
#        width = 6, height = 5)
# 
# SIBO_study_sum_pep_sample_nonSIBO <- SIBO_study_sum_pep_sample[V41 == "non-SIBO"]
# 
# SIBO_study_sum_pep_sample_SIBO <- SIBO_study_sum_pep_sample[V41 == "SIBO-replicate"]

# library("rcompanion")
# 
# hist_all <- plotNormalHistogram(SIBO_study_sum_pep_sample$Freq_pep,main = "all")
# hist_SIBO <- plotNormalHistogram(SIBO_study_sum_pep_sample_SIBO$Freq_pep, main = "SIBO")
# hist_non_SIBO <- plotNormalHistogram(SIBO_study_sum_pep_sample_nonSIBO$Freq_pep, main = "non-SIBO")
#save visualizations manually 

# ###### fisher test for peptidases in SIBO vs non-SIBO
# SIBO_pep1 <- SIBO_pep
# SIBO_pep1$V41[SIBO_pep1$V41 == "SIBO"] <- "ill-SIBO"
# SIBO_study1 <- SIBO_study 
# SIBO_study1$V41[SIBO_study1$V41 == "SIBO-replicate"] <- "ill-SIBO"
# 
# pep <- unique(SIBO_pep1$code)
# #View(SIBO_study)
# 
# k <- 1
# list <- list()
# for(pi in pep) {
#   tmp_pep <- SIBO_pep1[code == pi]
#   test6 <- table(tmp_pep$code, tmp_pep$V41)
#   if(ncol(test6) == 2) {
#     testa <- test6[,1]
#     testp <- test6[,2]
#     fSIBO <- nrow(SIBO_study1[SIBO_study1$V41 == "ill-SIBO"])
#     tSIBO <- nrow(SIBO_study1[SIBO_study1$V41 == "non-SIBO"])
#     nfSIBO <- fSIBO - testa
#     ntSIBO <- tSIBO - testp
#     nopep <- c(nfSIBO, ntSIBO)
#     newtable <- rbind(test6, nopep)
#     if(sum(newtable) == 3410) { #quality control
#       fishertest <- fisher.test(newtable)
#       fisher <- data.table(Peptidase = pi,
#                            fisher.p = fishertest$p.value,
#                            fisher.or = fishertest$estimate)
#       list[[k]] <- fisher
#       k <- k+1
#     }
#   }
# }
# 
# fisher_pep_SIBO <- rbindlist(list) 
# fisher_pep_SIBO
# fisher_pep_SIBO[, padjust := p.adjust(fisher.p, method = "fdr")]
# fisher_pep_SIBO[padjust < 0.05, sign.label := "*"]
# fisher_pep_SIBO$family <- substr(fisher_pep_SIBO$Peptidase, 1,1)
# fisher_pep_SIBO$SIBO <- "SIBO"
# fisher_pep_SIBO[, fisher.or.log2 := log2(fisher.or)]
# 
# pep_SIBo_nonSIBO <- ggplot(fisher_pep_SIBO, aes(SIBO,Peptidase, fill = fisher.or.log2)) +
#   geom_tile() +
#   geom_point(aes(shape = sign.label), size = 0.5, show.legend = FALSE) +
#   scale_shape_manual(values = 8, na.translate = FALSE) +
#   scale_fill_gradient2(high = "#ca0020", mid = "white", low = "#0571b0") +
#   labs(shape = "",fill = expression(log[2]~'(odds ratio)')) +
#   theme(legend.position = "bottom",
#         legend.justification = 	1,
#         axis.text.y = element_text(color = "black", angle =0, hjust = 1, size = 9),
#         axis.text.x = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.ticks.x= element_blank()) +
#   theme(axis.title.x = element_blank()) +
#   theme(panel.background = element_blank()) +
#   ggh4x::facet_nested("Family of Peptidases" + family ~., scales = "free", space = "free",nest_line = element_line(linetype = 1)) +
#   labs(y="Peptidases") 
# 
# pep_SIBo_nonSIBO 
# 
# ggsave("Output/Fisher_pep_SIBO.pdf", plot = pep_SIBo_nonSIBO,
#        width = 8.5, height = 3)



test1 <- HRGM_m8_SIBO[,c(1,2)]



library(dplyr)
library(tidyverse)
test3 <- test1 %>%
  group_by(qaccver, saccver) %>%
  summarise(Count = n())

test3


test3 <- test1 %>%
  # count number of rows for each combination of server_id and protocol
  group_by(saccver, qaccver) %>%
  tally() %>%
  # pivot the protocol names over the columns
  pivot_wider(names_from=saccver, values_from=n)

test3 <- as.data.frame(test3)
test3[is.na(test3)] <- 0
rownames(test3) <- test3$qaccver
test3 <- test3[,-1]
test3$count_pep <- rowSums(test3)
test3$Genomes <- rownames(test3)

Auxotrophy$count <- rowSums(Auxotrophy == 0)

Auxo_pep <- merge(Auxotrophy, test3, by.x="Genomes", by.y= "Genomes")
Auxo_pep

auxo_pep_all <- merge(Auxo_pep, SIBO_study, by.x="Genomes", by.y="Genomes")


new1 <- auxo_pep_all[ ,sum(count*Freq), by = sample]
colnames(new1) <- c("sample", "Abundance_Auxos")
new2 <- auxo_pep_all[ ,sum(count_pep*Freq), by = sample]
colnames(new2) <- c("sample", "Abundance_Pept")
new <- merge(new1, new2, by.x="sample", by.y="sample")
new

##SIBO
library(ggpubr)
colnames(SIBO)[34] <- "Name2"
colnames(SIBO)[12] <- "Name3"
new4 <- merge(SIBO, new, by.x="sample", by.y="sample")

###wilcox-Test 

#abundance-weighted average of auxotrophies
new4$diseasestatus <- ifelse(new4$V41 == "SIBO-replicate", 1, 0)
new4$V41[new4$V41 == "SIBO-replicate"] <- "SIBO"
new4$diseasestatus <- as.numeric(new4$diseasestatus)

wilcox.test(Abundance_Auxos ~diseasestatus, data = new4)

abun_auxo_SIBO <- ggplot(new4, aes(V41, Abundance_Auxos)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.3, width = 0.3) +
  stat_compare_means(method = "wilcox", label.x = 1.5,aes(label = sprintf("p = %5.4f", as.numeric(..p.format..)))) +
  theme_bw() +
  xlab("Disease state") +
  ylab("Abundance-weighted average of auxotrophies")
  
abun_auxo_SIBO

ggsave("Output/Fig_abun_auxos_SIBO.pdf", plot = abun_auxo_SIBO,
       width = 7, height = 4.5)

#abundance-weighted average of peptidases
wilcox.test(Abundance_Pept ~ diseasestatus, data = new4)

abun_pept_SIBO <- ggplot(new4, aes(V41, Abundance_Pept)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.3, width = 0.3) +
  stat_compare_means(method = "wilcox", label.x = 1.5,aes(label = sprintf("p = %5.4f", as.numeric(..p.format..)))) +
  theme_bw() +
  xlab("Disease state") +
  ylab("Abundance-weighted average of peptidases")
abun_pept_SIBO


ggsave("Output/Fig_abun_pept_SIBO.pdf", plot = abun_pept_SIBO,
       width = 7, height = 4.5)


# ###quotient abundance-weighted average of auxotrophies/peptidases
# new4$quotient <- new4$Abundance_Auxos/new4$Abundance_Pept
# 
# wilcox.test(quotient ~ diseasestatus, data = new4)
# 
# abun_pept_auxos_SIBO <- ggplot(new4, aes(V41, quotient)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(alpha = 0.3, width = 0.3) +
#   stat_compare_means(method = "wilcox", label.x = 0.5) +
#   theme_bw() +
#   xlab("Disease state") +
#   ylab("Abundance-weighted average of auxotrophies/
#        Abundance-wegihted average of peptidases")
# abun_pept_auxos_SIBO
# 
# 
# ggsave("Output/Fig_abun_auxos_pept_SIBO.pdf", plot = abun_pept_auxos_SIBO,
#        width = 7, height = 4.5)
# 
# ##barplot



# ###Co-correlation
# library(cocor)
# SIBO <- new4[new4$V41 == "SIBO"]
# nSIBO <- new4[new4$V41 == "non-SIBO"]
# new5 <- list(SIBO, nSIBO)
# cocor(~Abundance_Auxos+ Abundance_Pept | Abundance_Auxos + Abundance_Pept, new5)
# cor(SIBO$Abundance_Auxos, SIBO$Abundance_Pept)
# cor(nSIBO$Abundance_Auxos, nSIBO$Abundance_Pept)
