#### protelytic activity HRGM ###
source("Scripts/init_models_HRGM.R")

source("Scripts/Predict_auxotrophies.R")

source("Scripts/Auxotable_melted_merged.R")

library(tidyverse)

HRGM_m8 <- fread("Data/HRGM_all.m8")
HRGM_m8
colnames(HRGM_m8) <- c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
                    "evalue", "bitscore")
HRGM_m8

## shorten length of HRGM_Genome name 
HRGM_m8$qaccver <- substr(HRGM_m8$qaccver, 1, 16)

##  

#HRGM_merops <- merge(meta, HRGM_m8, by.x= "Genomes", by.y= "qaccver")


##peptidases found

domain <-fread("Data/domain.csv")
domain

domain_results <- merge(domain, HRGM_m8, by.x="mernum", by.y="saccver")

##display number of peptidases in HRGM
number_pept_HRGM <- unique(domain_results$mernum)

##delete all inhibitors
library(tidyverse)
domain_results <- domain_results[type != "inhibitor",]
#View(domain_results)
domain_results$family <- substr(domain_results$code, 1,1)
##get count of peptidases overall in HRGM catalogue
number_genomes_pep <- unique(domain_results$qaccver)

##visualization
##visualization

z <- domain_results %>%
  arrange(family, qaccver) %>%
  mutate(mernum = fct_inorder(mernum)) %>%
  ggplot(domain_results, mapping= aes(x=qaccver, y=mernum, colour= family, size= 0.3))+
  geom_point( size =0.3) +
  theme(axis.text.x=element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  ylab("Peptidases") +
  xlab("HRGM genomes") +
  scale_colour_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00"),
                      labels =  c("A" = "Aspartic peptidase", "C" = "Cysteine peptidase", "M" = "Metallopeptidase", "S" = "Serine peptidase",
                                  "T" = "Threonine peptidase", "U" = "Unknown catalytic type"),
                      breaks = c("A", "M", "C", "S", "T", "U"))
z <- z + labs(colour = "Family")
z
# 
# ggsave("Output/HRGM_peptidases.pdf", plot = z,
#        width = 12, height = 7)

###preediction of secreted peptidases
prediction_SP <- read.delim("Data//prediction_results.txt")
# prediction_SP_1 <- read.delim("/mnt/nuuk/2022/HRGM/Protein_sequences_Reimagine/SignalP/prediction_results.txt")


prediction_SP$Mernum <- substr(prediction_SP$X..ID, 1, 10)

secpep <- merge(domain_results, prediction_SP, by.x= "mernum", by.y = "Mernum")
secpep <- as.data.table(secpep)
secpep_new <- secpep[secpep$Prediction != "NO_SP"]

secr_pep_HRGM <- ggplot(secpep_new, aes(qaccver, code))+
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, colour = "black"),
        axis.text.y = element_text(colour = "black")) +
  labs( x= "Genomes", y ="Peptidases")
secr_pep_HRGM 
# ggsave("Output/HRGM_peptidases_secreted.pdf", plot = secr_pep_HRGM,
#        width = 12, height = 7)

###analyze peptidases per genome in relation to auxotrophies per genome
domain_results_filt <- domain_results[,c("mernum","qaccver")]
library(tidyverse)
sum_pep <- domain_results_filt %>%
  # count number of rows for each combination of server_id and protocol
  group_by(mernum, qaccver) %>%
  tally() %>%
  # pivot the protocol names over the columns
  pivot_wider(names_from=mernum, values_from=n)

sum_pep <- as.data.frame(sum_pep)
sum_pep[is.na(sum_pep)] <- 0
rownames(sum_pep) <- sum_pep$qaccver
sum_pep <- sum_pep[,-1]
sum_pep$count_pep <- rowSums(sum_pep)
sum_pep$Genomes <- rownames(sum_pep)

#View(domain_results)
#get auxotrophy info
Auxotrophy$count <- rowSums(Auxotrophy == 0)

Aux_pep <- merge(Auxotrophy, sum_pep, by.x="Genomes", by.y= "Genomes")
Aux_pep
Aux_pep<- Aux_pep[, c("Genomes", "count", "count_pep")]
#View(Aux_pep)

Aux_pep$count_pep <- as.numeric(Aux_pep$count_pep)
Aux_pep$count <- as.numeric(Aux_pep$count)
cor.test(Aux_pep$count_pep, Aux_pep$count, method = "kendall")


library(ggpubr)
# numb_auxos_pep <- ggplot(Aux_pep, aes(count, count_pep))+
#   geom_point() +
#   geom_smooth()+
#   xlab("Number of auxotrophies") +
#   ylab("Number of peptidases") +
#   theme_bw() +
#   theme(panel.background = element_blank()) +
#   stat_cor(method = "kendall", size = 3)
#   
# numb_auxos_pep 
# 
# ggsave("Output/Corr_HRGM_auxo_pep.pdf", plot = numb_auxos_pep,
#        width = 12, height = 7)


library(ggh4x)
# numb_pep_HRGM <- ggplot(Aux_pep, aes(x=count_pep, y=Genomes))+
#   geom_point(aes(colour = factor(count)), size = 0.4) +
#   xlab("HRGM genomes") +
#   ylab("Number of peptidases") +
#   theme(axis.text.y = element_blank(),
#         axis.ticks.y = element_blank()) +
#   scale_color_manual(values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
#                                 "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888",
#                                 "#E69F00", "#56B4E9", "#009E73", 
#                                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
#   #facet_grid(. ~count, space="free") +
#   theme(legend.position = "none") +
#   #ggh4x::facet_nested(.~ "Number of auxotrophies" + count,  nest_line = element_line(linetype = 1)) +
#   theme(strip.background=element_rect(color="grey30", fill="grey90"))
# numb_pep_HRGM


# ggsave("Output/HRGM_peptidases_numb.pdf", plot = numb_pep_HRGM,
#        width = 12, height = 7)


# ### get percentage of peptidases per genome in HRGM catalogue
# pept <- unique(Aux_pep$count_pep)
# perc_pepti <- list()
# k <- 1
# 
# for(pepti in pept) {
#   peptid <- Aux_pep[Aux_pep$count_pep == pepti]
#   perc_pep <- (nrow(peptid) / 3687) * 100
#   perc_pept <- data.table(Perc_pept = perc_pep,
#                           Pept_count = pepti)
#   perc_pepti[[k]] <- perc_pept
#   k <- k + 1
#   
# }
# perc_peptid <- rbindlist(perc_pepti)
# perc_peptid
# Perc_pept <- ((3687 - nrow(Aux_pep))/3687) * 100
# Pept_count <- 0
# zero_pept3 <- data.frame(Perc_pept, Pept_count)
# perc_peptid1 <- rbind(perc_peptid, zero_pept3)
# 
# ###visualization
# pept_perc <- ggplot(perc_peptid1, aes(Pept_count, Perc_pept)) +
#   geom_bar(stat = "identity") +
#   labs(y = "Peptidases per genome [%]", x = "Number of peptidases") +
#   theme_bw()
# 
# ggsave("Output/HRGM_peptidases_numb_perc.pdf", plot = pept_perc,
#        width = 12, height = 7)

# 
# tab <- table(Pep = Aux_pep$count_pep, Auxos = Aux_pep$count)
# tab <- data.frame(tab)
# #View(tab)


##correlation number of auxotrophies with peptidases
# Aux_peptidases_log_reg <- merge(Aux_pep, domain_results, by.x= "Genomes", by.y="qaccver")
# pep <- unique(Aux_peptidases_log_reg$code)
# l <- list()
# k <- 1
# 
# for(pepi in pep){
#   tmp <- Aux_peptidases_log_reg[code == pepi]
#   log.model <- glm(count ~ Prototrophy, data = tmp, family = 'binomial')
#   tmp_3 <- data.table(code = pepi,
#                       Peptidase = rownames(pvalue),
#                       pvalue = coef(summary(log.model))[,4],
#                       Estimate  = coef(summary(log.model))[,1])
#   l[[k]] <- tmp_3
#   k <- k+1
# }
# 


#########auxotrophies and peptidases


auxo_pep <- merge(domain_results, Auxotrophy_2, by.x="qaccver", by.y="Genomes")

##filter for abundant peptidases
test <- table(auxo_pep$code)
test2 <- data.table(test)
filt_pep <- test2[test2$N > 500]
filt_pep1 <- filt_pep$V1

##delete auxotrophies
auxo_pep3 <- auxo_pep[auxo_pep$code %in% filt_pep1]
auxo_pep3 <- auxo_pep3[Compound != "Gly"]
auxo_pep3 <- auxo_pep3[Compound != "Glu"]
auxo_pep3 <- auxo_pep3[Compound != "Asp"]
auxo_pep3 <- auxo_pep3[Compound != "Ala"]

auxo_pep3$Prototrophy <- as.factor(auxo_pep3$Prototrophy)
# aux <- unique(auxo_pep3$Compound)
# l <- list()
# k <- 1
# 
# for(a in aux) {
#     tmp_2 <- auxo_pep3[Compound == a]
#     log.model <- glm(Prototrophy ~ code, data = tmp_2, family = 'binomial')
#     pvalue <- data.frame(coef(summary(log.model))[,4]) #get code for peptidase
#     tmp_3 <- data.table(AA = a,
#                         Peptidase = rownames(pvalue),
#                         pvalue = coef(summary(log.model))[,4],
#                         Estimate  = coef(summary(log.model))[,1])
#     l[[k]] <- tmp_3
#     k <- k+1
# }
# 
# auxo_pep_lin <- rbindlist(l)
# auxo_pep_lin
# 
# auxo_pep_lin <- auxo_pep_lin[auxo_pep_lin$Peptidase != "(Intercept)"]
# auxo_pep_lin$Peptidase <- substring(auxo_pep_lin$Peptidase, 5)
# auxo_pep_lin[, padjust := p.adjust(pvalue, method = "fdr")]
# auxo_pep_lin[padjust < 0.05, sign.label := "*"]
# 
# ###visualization
# auxo_lin <- ggplot(auxo_pep_lin, aes(AA, Peptidase, fill = Estimate)) +
#   geom_tile() +
#   geom_point(aes(shape = sign.label), size = 0.5, show.legend = FALSE) +
#   scale_shape_manual(values = 8, na.translate = FALSE) +
#   scale_fill_gradient2(high = "#ca0020", mid = "white", low = "#0571b0") +
#   labs(x = "Auxotrophy", y = "Peptidases", shape = "") +
#   theme(legend.position = "bottom",
#         legend.justification = 	1,
#         axis.text.x = element_text(color = "black", angle = 45, hjust = 1, size = 9),
#         axis.text.y = element_text(color = "black", size = 9)) +
#   theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b= 0, l = 0))) +
#   theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b= 0, l = 0))) +
#   theme(panel.background = element_blank())
# 
# 
# auxo_lin  
# 
# ggsave("Output/HRGM_peptidases_log_reg.pdf", plot = auxo_lin,
#        width = 8, height = 7)

###fisher test


test4 <- auxo_pep[code =="S24.003"]
test5 <- test4[Compound == "Trp"]
test5

##filter for abundant peptidases
test <- table(auxo_pep$code)
test2 <- data.table(test)
filt_pep <- test2[test2$N > 500]
filt_pep1 <- filt_pep$V1
auxo_pep$familypep <-  substr(auxo_pep$code, 1,4)

##delete auxotrophies
auxo_pep4 <- auxo_pep[auxo_pep$code %in% filt_pep1]
auxo_pep4 <- auxo_pep4[Compound != "Gly"]
auxo_pep4 <- auxo_pep4[Compound != "Glu"]
auxo_pep4 <- auxo_pep4[Compound != "Asp"]
auxo_pep4 <- auxo_pep4[Compound != "Ala"]

pep <- unique(auxo_pep4$code)
aa <- unique(auxo_pep4$Compound)
k <- 1
list <- list()
for(p in pep) {
  for(a in aa) {
    tmp_aux_p <- auxo_pep4[code == p]
    tmp_aux_p2 <- tmp_aux_p[Compound == a]
    test6 <- table(tmp_aux_p2$code, tmp_aux_p2$Prototrophy)
    if(ncol(test6) == 2) {
    testa <- test6[,1]
    testp <- test6[,2]
    Auxotrophy_tmp <- Auxotrophy_2[Compound == a]
    auxo <- nrow(Auxotrophy_tmp[Auxotrophy_tmp$Prototrophy == 0])
    proto <- nrow(Auxotrophy_tmp[Auxotrophy_tmp$Prototrophy == 1])
    nauxo <- auxo - testa 
    nproto <- proto - testp
    nopep <- c(nauxo, nproto)
    newtable <- rbind(test6, nopep)
    if(sum(newtable) == 3687) { #quality control
    fishertest <- fisher.test(newtable)
    fisher <- data.table(AA = a,
                         Peptidase = p,
                         fisher.p = fishertest$p.value,
                         fisher.or = fishertest$estimate)
    list[[k]] <- fisher
    k <- k+1
    }
    }
  }
}

fisher_pep <- rbindlist(list)
fisher_pep
fisher_pep[, fisher.padj := p.adjust(fisher.p, method = "fdr")]
fisher_pep[, fisher.or.log2 := -log2(fisher.or)]
fisher_pep[fisher.padj < 0.05, sign.label1 := "Padj < 0.05"]

fisher_pep$family <- substr(fisher_pep$Peptidase, 1,1)

p <- ggplot(fisher_pep, aes(AA, Peptidase, fill = -log2(fisher.or))) +
  geom_tile() +
  geom_point(aes(shape = sign.label1), size = 1) +
  scale_fill_gradient2(high = "#ca0020", mid = "white", low = "#0571b0") +
  scale_shape_manual(values = 8, na.translate = FALSE) +
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(x = "Auxotrophy", y = "Peptidases", shape = "",
       fill = expression(log[2]~'(odds ratio)')) +
  theme(legend.position = "bottom",
        legend.justification = 	1,
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(color = "black", hjust = 1, size = 8),
        axis.ticks.y = element_blank()) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b= 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b= 0, l = 0))) +
  theme(panel.background = element_blank()) +
  ggh4x::facet_nested("Family of Peptidases" + family  ~., scales = "free", space="free", nest_line = element_line(linetype = 1))
  #ggh4x::facet_nested( "Family of Peptidases" + family~.,  nest_line = element_line(linetype = 1))
p

# ggsave("Output/Peptidases_fisher_auxo.pdf", plot = p,
#        width = 10, height =15)

# 
# 
# 
# 
# 
# auxo_pep4 <- auxo_pep
# auxo_pep4 <- auxo_pep4[Compound != "Gly"]
# auxo_pep4 <- auxo_pep4[Compound != "Glu"]
# auxo_pep4 <- auxo_pep4[Compound != "Asp"]
# auxo_pep4 <- auxo_pep4[Compound != "Ala"]
# 
# pep <- unique(auxo_pep4$familypep)
# aa <- unique(auxo_pep4$Compound)
# k <- 1
# list <- list()
# for(p in pep) {
#   for(a in aa) {
#     tmp_aux_p <- auxo_pep4[familypep == p]
#     tmp_aux_p2 <- tmp_aux_p[Compound == a]
#     test6 <- table(tmp_aux_p2$familypep, tmp_aux_p2$Prototrophy)
#     if(ncol(test6) == 2) {
#       testa <- test6[,1]
#       testp <- test6[,2]
#       Auxotrophy_tmp <- Auxotrophy_2[Compound == a]
#       auxo <- nrow(Auxotrophy_tmp[Auxotrophy_tmp$Prototrophy == 0])
#       proto <- nrow(Auxotrophy_tmp[Auxotrophy_tmp$Prototrophy == 1])
#       nauxo <- auxo - testa 
#       nproto <- proto - testp
#       nopep <- c(nauxo, nproto)
#       newtable <- rbind(test6, nopep)
#       if(sum(newtable) == 3687) { #quality control
#         fishertest <- fisher.test(newtable)
#         fisher <- data.table(AA = a,
#                              Peptidase = p,
#                              fisher.p = fishertest$p.value,
#                              fisher.or = fishertest$estimate)
#         list[[k]] <- fisher
#         k <- k+1
#       }
#     }
#   }
# }
# 
# fisher_pep <- rbindlist(list)
# fisher_pep
# fisher_pep[, fisher.padj := p.adjust(fisher.p, method = "fdr")]
# fisher_pep[, fisher.or.log2 := -log2(fisher.or)]
# fisher_pep[fisher.padj < 0.05, sign.label1 := "Padj < 0.05"]
# 
# fisher_pep$family <- substr(fisher_pep$Peptidase, 1,1)
# 
# p_fam <- ggplot(fisher_pep, aes(AA, Peptidase, fill = -log2(fisher.or))) +
#   geom_tile() +
#   geom_point(aes(shape = sign.label1), size = 1) +
#   scale_fill_gradient2(high = "#ca0020", mid = "white", low = "#0571b0") +
#   scale_shape_manual(values = 8, na.translate = FALSE) +
#   scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
#   labs(x = "Auxotrophy", y = "Peptidases", shape = "",
#        fill = expression(log[2]~'(odds ratio)')) +
#   theme(legend.position = "bottom",
#         legend.justification = 	1,
#         axis.text.x = element_text(color = "black", angle = 45, hjust = 1, size = 8),
#         axis.text.y = element_text(color = "black", hjust = 1, size = 8),
#         axis.ticks.y = element_blank()) +
#   theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b= 0, l = 0))) +
#   theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b= 0, l = 0))) +
#   theme(panel.background = element_blank()) +
#   ggh4x::facet_nested("Family of Peptidases" + family  ~., scales = "free", space="free", nest_line = element_line(linetype = 1))
# #ggh4x::facet_nested( "Family of Peptidases" + family~.,  nest_line = element_line(linetype = 1))
# p_fam
# 
# ggsave("Output/Peptidasesfamilies_fisher_auxo.pdf", plot = p_fam,
#        width = 10, height =15)

# Auxo_pep5 <- merge(Auxotrophy, domain_results, by.x= "Genomes", by.y="qaccver")
# Auxo_pep5$peptidase_Status <- 1
# 
# 
# 
# Auxotrophy$peptidase_Status <- 0
# 
# 
# pep <- unique(filt_pep$V1)
# k <- 1
# l <- list()
# 
# for(pi in pep) {
#   tmp_aux_p <- Auxo_pep5[code == pi]
#   Genomes_0 <- unique(tmp_aux_p$Genomes)
#   tmp_aux_p2 <- tmp_aux_p[,c(1, 23, 24)]
#   
#   Auxo_pep6 <- Auxotrophy[!Auxotrophy$Genomes %in% Genomes_0]
#   Auxo_pep7 <- Auxo_pep6[,c(22,23,24)]
#   
#   new2 <- rbind(Auxo_pep7,tmp_aux_p2)
#   new2 <- new2[!duplicated(new2),]
#   
#   log.model <- glm(peptidase_Status ~ count, data = new2, family = 'binomial')
#   pvalue <- data.frame(coef(summary(log.model))[,4]) #get code for peptidase
#   tmp_3 <- data.table(Peptidase = pi,
#                       names = rownames(pvalue),
#                       pvalue = coef(summary(log.model))[,4],
#                       Zscore  = coef(summary(log.model))[,3])
#   l[[k]] <- tmp_3
#   k <- k+1
# }
# tmp_4 <- rbindlist(l)
# tmp_4
# 
# tmp_4 <- tmp_4[tmp_4$names != "(Intercept)"]
# tmp_4[, padjust := p.adjust(pvalue, method = "fdr")]
# tmp_4[padjust < 0.05, sign.label := "*"]
# tmp_4$family <- substr(tmp_4$Peptidase, 1,1)
# 
# 
# ###visualization
# auxo_log <- ggplot(tmp_4, aes(Peptidase, names, fill = Zscore)) +
#   geom_tile() +
#   geom_point(aes(shape = sign.label), size = 0.5, show.legend = FALSE) +
#   scale_shape_manual(values = 8, na.translate = FALSE) +
#   scale_fill_gradient2(high = "#ca0020", mid = "white", low = "#0571b0") +
#   labs(shape = "") +
#   theme(legend.position = "bottom",
#         legend.justification = 	1,
#         axis.text.x = element_text(color = "black", angle = 45, hjust = 1, size = 9),
#         axis.text.y = element_blank()) +
#   theme(axis.title.y = element_blank()) +
#   theme(axis.title.x = element_blank()) +
#   theme(panel.background = element_blank()) +
#   ggh4x::facet_nested(~ "Family of Peptidases" + family, scales = "free", space = "free",nest_line = element_line(linetype = 1))
# 
# 
# auxo_log
# 
# ggsave("Output/Log.Regres.pdf", plot = auxo_log,
#        width = 8.5, height = 3)
