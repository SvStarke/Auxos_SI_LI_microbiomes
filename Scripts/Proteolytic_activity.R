library(data.table)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(vegan)
library(rlist)
library(ape)

##load other scripts
source("Scripts/Mapping_ASV2HRGM.R")

source("Scripts/Predict_auxotrophies.R")


##proteolytic activity
test <- fread("Data/all.m8")
test

colnames(test) <- c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
   "evalue", "bitscore")
test

test$qaccver <- substr(test$qaccver, 1, 16)

reimagine_merops <- reimagine

reimagine_merops <- merge(reimagine, test, by.x= "Genomes", by.y= "qaccver")

#View(reimagine_merops)

ggplot(reimagine_merops, aes(sample, saccver)) +
  geom_point(aes(colour = factor(location)))

##peptidases found

domain <-fread("Data/domain.csv")
domain

domain_results <- merge(domain, test, by.x="mernum", by.y="saccver")

#View(domain_results)
##delete all inhibitors

domain_results <- domain_results[type != "inhibitor",]

domain_results_p <- merge(reimagine, domain_results, by.x="Genomes", by.y = "qaccver")
domain_results_p$family <- substr(domain_results_p$code, 1,1)

##visualization
##visualization
# 
# t <- ggplot(domain_results_p, aes(sample, mernum)) +
#   geom_point(aes (colour = factor(family)), size =0.4) +
#   theme_bw() +
#   theme(axis.text.x=element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks = element_blank()) +
#   ylab("Peptidases") +
#   xlab("Samples") +
#   scale_colour_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00"),
#                       labels =  c("A" = "Aspartic peptidase", "C" = "Cysteine peptidase", "M" = "Metallopeptidase", "S" = "Serine peptidase",
#                                   "T" = "Threonine peptidase", "U" = "Unknown catalytic type"),
#                       breaks = c("A", "M", "C", "S", "T", "U"))
# 
# t <- t + labs(colour = "family")  + facet_grid(. ~ location, scales = "free", space="free")
# t
# 
# ggsave("Output/Family.pdf", plot = t,
#        width = 12, height = 7)


##new figure with sorted colours in plot
domain_results_p$location[domain_results_p$location == "Ileum"] <- "Farthest distance"
domain_results_p$location[domain_results_p$location == "Stool"] <- "Colon"
domain_results_p$location_sort <- factor(domain_results_p$location, levels=c("Duodenum", "Jejunum", "Farthest distance", "Colon"))
t_new <- domain_results_p %>%
  arrange(family, sample) %>%
  mutate(mernum = fct_inorder(mernum)) %>%
  ggplot(domain_results_p, mapping= aes(x=sample, y=mernum, colour= family))+
  geom_point(size = 0.5)  +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "white"),
        strip.text.x = element_text(colour = "black", size = 13),
        axis.title.x = element_text(size = 13, colour = "black"),
        axis.title.y = element_text(size = 13, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 13, colour = "black")) +
  ylab("Peptidases") +
  xlab("Samples") +
  scale_colour_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00"),
                      labels =  c("A" = "Aspartic peptidase", "C" = "Cysteine peptidase", "M" = "Metallopeptidase", "S" = "Serine peptidase",
                                  "T" = "Threonine peptidase", "U" = "Unknown catalytic type"),
                      breaks = c("A", "M", "C", "S", "T", "U"))
t_new <- t_new + labs(colour = "Family")  + facet_grid(. ~ location_sort, scales = "free", space="free") +
  guides(colour = guide_legend(override.aes = list(size=3)))
t_new


# ggsave("Output/Family_sorted.pdf", plot = t_new,
#        width = 12, height = 7)

##calculate Bray Curtis distance
peptidases_bray_curtis <- domain_results_p
peptidases_bray_curtis$location <- as.character(peptidases_bray_curtis$location)
peptidases_bray_curtis$location[peptidases_bray_curtis$location == "Duodenum"] <- "Small intestine"
peptidases_bray_curtis$location[peptidases_bray_curtis$location == "Jejunum"] <- "Small intestine"
peptidases_bray_curtis$location[peptidases_bray_curtis$location == "Farthest distance"] <- "Small intestine"
#peptidases_bray_curtis$location[peptidases_bray_curtis$location == "Ileum"] <- "Small intestine"
peptidases_bray_curtis$location[peptidases_bray_curtis$location == "Colon"] <- "Large intestine"

#ggplot(peptidases_bray_curtis, aes(location, code, colour = family)) +
  geom_bar(stat = "identity")

pep <- unique(peptidases_bray_curtis$code)
pep_list <- list()
k <- 1

for(pi in pep) {
  tmp_pep <- peptidases_bray_curtis[code == pi]
  vectortmp <- tmp_pep$location
  vectortmp <- unique(vectortmp)
  if(length(vectortmp) == 2) {
    tmp_pep$location <- "Both"
  }
  pep_list[[k]] <- tmp_pep
  k <- k+1
}

pep_shared_SI_LI <- rbindlist(pep_list)

pept <- unique(pep_shared_SI_LI$code)
loc <- unique(pep_shared_SI_LI$location)
pep_list1 <- list()
k <- 1

for(loci in loc) {
  vector_pep <- as.numeric(length(pept))
  tmp_loc <- pep_shared_SI_LI[location == loci]
    tmp_loc_pep <- unique(tmp_loc$code)
    length_tmp_pep_loc <- as.numeric(length(tmp_loc_pep))
    perc <- (length_tmp_pep_loc/vector_pep)*100
    loc_pep <- data.table(Location = loci,
                          Perc = perc)
    pep_list1[[k]] <- loc_pep
    k <- k+1
  
}

perc_pep_loc <- rbindlist(pep_list1)
perc_pep_loc
perc_pep_loc$Nr <- 1

# per_pep_loc <- ggplot(perc_pep_loc, aes(Nr,Perc)) +
#   geom_col(aes(fill = Location)) +
#   scale_fill_manual(values = c("#DDCC77", "#000000", "#888888"),
#                       labels =  c("Small Intestine" = "Small intestine", "Large intestine" = "Large intestine", "Both" = "Both")) +
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
# ggsave("Output/Percentage_Peptidases.pdf", plot = per_pep_loc,
#        width = 4, height = 3)

# stacked_bar_plot <- ggplot(pep_shared_SI_LI, aes(location, code)) +
#   geom_col(aes(colour = family)) +
#   xlab("Location") +
#   ylab("Peptidase") +
#   scale_colour_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00"),
#                       labels =  c("A" = "Aspartic peptidase", "C" = "Cysteine peptidase", "M" = "Metallopeptidase", "S" = "Serine peptidase",
#                                   "T" = "Threonine peptidase", "U" = "Unknown catalytic type")) +
#   theme(axis.text.x = element_text(colour = "black"),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         panel.grid =  element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black")) +
#   guides(col = guide_legend(title = "Family",override.aes = list(shape = 15, size = 8, fill = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2"))))+
#   theme(legend.position = "none")
# ggsave("Output/Stacked_bar_Pep_location.pdf", plot = stacked_bar_plot,
#        width = 8, height = 5)


# peptidases_bray_curtis1 <- peptidases_bray_curtis[,c("code", "PatNo")]
# t(peptidases_bray_curtis1)

pep <- unique(peptidases_bray_curtis$code)

pep_list <- list()
k <- 1

for(pi in pep) {
  pep_tmp <- peptidases_bray_curtis[code == pi]
  pep_SI_LI <- pep_tmp %>% count(PatNo, location)
  pep_SI_LI$code <- pi
  pep_list[[k]] <- pep_SI_LI
  k <- k+1
}
pep_all_BC <- rbindlist(pep_list)
new_pep_all <- dcast(pep_all_BC, location + PatNo ~ code, value.var = "n")
new_pep_all[is.na(new_pep_all)] <- 0
tmp2 <- new_pep_all[,-c(1,2)]
dist <- vegdist(tmp2, method = "bray")
testx <- pcoa(dist)
dist_per <- as.matrix(dist)
biplot(testx)
bray_curtis_pcoa_df <- data.frame(pcoa1 = testx$vectors[,1], 
                                  pcoa2 = testx$vectors[,2])
bray_curtis_plot <- ggplot(data = bray_curtis_pcoa_df, aes(x=pcoa1, y=pcoa2, colour = new_pep_all$location)) +
  geom_point() +
  stat_ellipse() +
  labs(x = "PCoA1",
       y = "PCoA2") +
  theme(title = element_text(size = 10),
        legend.position = "bottom") +
  guides(colour = guide_legend(title = "Location")) +
  scale_color_manual(values = c("#661100", "#999933", "#117733", "332288")) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size =11)) +
  theme(axis.title.y = element_text(colour = "black", size = 13, margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 13, margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=10, colour = "black", vjust = 1,  angle = 0, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size =10, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) 

bray_curtis_plot


# ggsave("Output/PCoA_all4.pdf", plot = bray_curtis_plot,
       # width = 8, height = 5)

#statistical evaluation PERMANOVA
adonis2(dist_per ~ location, data = new_pep_all)


# ###all four locations
# 
# peptidases_bray_curtis <- domain_results_p
# 
# pep <- unique(peptidases_bray_curtis$code)
# 
# pep_list <- list()
# k <- 1
# 
# for(pi in pep) {
#   pep_tmp <- peptidases_bray_curtis[code == pi]
#   pep_SI_LI <- pep_tmp %>% count(PatNo, location)
#   pep_SI_LI$code <- pi
#   pep_list[[k]] <- pep_SI_LI
#   k <- k+1
# }
# pep_all_BC <- rbindlist(pep_list)
# new_pep_all <- dcast(pep_all_BC, location + PatNo ~ code, value.var = "n")
# new_pep_all[is.na(new_pep_all)] <- 0
# tmp2 <- new_pep_all[,-c(1,2)]
# dist <- vegdist(tmp2, method = "bray")
# testx <- pcoa(dist)
# 
# biplot(testx)
# bray_curtis_pcoa_df <- data.frame(pcoa1 = testx$vectors[,1], 
#                                   pcoa2 = testx$vectors[,2])
# bray_curtis_plot_all4 <- ggplot(data = bray_curtis_pcoa_df, aes(x=pcoa1, y=pcoa2, colour = new_pep_all$location)) +
#   geom_point() +
#   labs(x = "PCoA1",
#        y = "PCoA2") +
#   theme(title = element_text(size = 10),
#         legend.position = "right") +
#   guides(colour = guide_legend(title = "Location")) +
#   scale_color_manual(values = c("#661100", "#999933", "#117733", "332288")) +
#   theme(panel.background = element_blank(),
#         axis.line = element_line(colour = "black")) +
#   theme(legend.text = element_text(size=10)) +
#   theme(legend.title = element_text(size =11)) +
#   theme(axis.title.y = element_text(colour = "black", size = 13, margin = margin(0,10,0,0)))+
#   theme(axis.title.x = element_text(colour = "black", size = 13, margin = margin(10,0,0,0))) +
#   theme(axis.text.x = element_text(size=10, colour = "black", hjust = 1,  angle = 0, margin = margin(10,0,0,0))) +
#   theme(axis.text.y = element_text(size =10, colour = "black")) +
#   theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) 
# 
# bray_curtis_plot_all4
# 
# ggsave("Output/PCoA_all4.pdf", plot = bray_curtis_plot_all4,
#        width = 8, height = 5)

### sort peptidases to segments 
# t1 <- ggplot(domain_results_p, aes(sample, mernum)) +
#   geom_point(aes (colour = factor(family)), size = 0.5) +
#   theme_bw() +
#   theme(axis.text.x=element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks = element_blank(),
#         panel.grid.major = element_blank()) +
#   ylab("Peptidases") +
#   xlab("Samples") +
#   scale_colour_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00"),
#                       labels =  c("A" = "Aspartic peptidase", "C" = "Cysteine peptidase", "M" = "Metallopeptidase", "S" = "Serine peptidase",
#                                   "T" = "Threonine peptidase", "U" = "Unknown catalytic type"),
#                       breaks = c("A", "C", "M", "S", "T", "U"))
# 
# t1 <- t1 + labs(colour = "family")  + facet_grid(family ~ location_sort,  scales = "free")
# t1
# 
# 
# ggsave("Output/Family_per_segment.pdf", plot = t1,
#        width = 12, height = )

# ##location
# domain_results_p$location <- factor(domain_results_p$location, levels=c("Duodenum", "Jejunum", "Farthest distance", "Colon"))
# level_order <-c("Duodenum", "Jejunum", "Farthest distance", "Colon")
# u <- ggplot(domain_results_p, aes(sample, mernum)) +
#   geom_point(aes(colour = factor(location))) +
#   theme(axis.text.x=element_blank(),
#         axis.text.y=element_blank(),
#         legend.text = element_blank(),
#         axis.ticks = element_blank())+
#   ylab("Peptidases") +
#   xlab("Samples") +
#   scale_color_manual(values = c("#F0E442", "#0072B2", "#D55E00", "#CC79A7")) 
# u <- u+ theme(legend.position = "none") + facet_grid(. ~ location, scales = "free", space="free")
# 
# ggsave("Output/Location.pdf", plot = u,
#        width = 12, height = 6)

# 
# p <- ggplot(domain_results_p, aes(family, mernum)) +
#   geom_point(aes(colour = factor(location))) +
#   theme(axis.text.y=element_blank()) +
#   ylab("Peptidases") +
#   xlab("Samples")
# 
# p+ labs(colour="location")

##correlation between the number of auxotrophies (abundance-weighted) and the number of peptidases per sample

test1 <- test[,c(1,2)]



library(dplyr)
# test3 <- test1 %>%
#         group_by(qaccver, saccver) %>%
#           summarise(Count = n())




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

auxo_pep_all <- merge(Auxo_pep, reimagine, by.x="Genomes", by.y="Genomes")


new1 <- auxo_pep_all[ ,sum(count*Freq), by = sample]
new2 <- auxo_pep_all[ ,sum(count_pep*Freq), by = sample]
new <- merge(new1, new2, by.x="sample", by.y="sample")
new

cor.test(new$V1.x, new$V1.y)

source("Scripts/Reimagine_study.R")

new_w_regions <- merge(new, meta, by.x="sample", by.y="sample")
new_w_regions
new_w_regions$location[new_w_regions$location == "Ileum"] <- "Farthest distance"
new_w_regions$location[new_w_regions$location == "Stool"] <- "Colon"

###Jejunum
new_w_regions_j <- new_w_regions[location == "Jejunum"]

j <- ggplot(new_w_regions_j, aes(V1.y, V1.x)) +
  geom_smooth(method=lm, col = "black") +
  geom_point(shape = 21) +
  theme_bw() +
  theme(axis.title.x = element_text(size=8),
        axis.title.y = element_text(size =8)) +
  xlab("Abundance-weighted average of peptidases") +
  ylab("Abundance-weighted average of auxotrophies") +
  facet_grid(. ~ location, scales = "free", space="free") +
  stat_cor(method = "spearman", size = 3, cor.coef.name = "rho",p.accuracy = 0.0001)
j

####Duodenum
new_w_regions_d <- new_w_regions[location == "Duodenum"]

d <- ggplot(new_w_regions_d, aes(V1.y, V1.x)) +
  geom_smooth(method=lm, col = "black") +
  geom_point(shape = 21) +
  theme_bw() +
  theme(axis.title.x = element_text(size=8),
        axis.title.y = element_text(size =8)) +
  xlab("Abundance-weighted average of peptidases") +
  ylab("Abundance-weighted average of auxotrophies") +
  facet_grid(. ~ location, scales = "free", space="free") +
  stat_cor(method = "spearman", size = 3, cor.coef.name = "rho",p.accuracy = 0.0001)
d

###Ileum
new_w_regions_i <- new_w_regions[location == "Farthest distance"]

i <- ggplot(new_w_regions_i, aes(V1.y, V1.x)) +
  geom_smooth(method=lm, col = "black") +
  geom_point(shape = 21) +
  theme_bw() +
  theme(axis.title.x = element_text(size=8),
        axis.title.y = element_text(size =8)) +
  xlab("Abundance-weighted average of peptidases") +
  ylab("Abundance-weighted average of auxotrophies") +
  facet_grid(. ~ location, scales = "free", space="free")+
  stat_cor(method = "spearman", size = 3, cor.coef.name = "rho", label.y = 6.8,p.accuracy = 0.0001)
  
i

###Stool
new_w_regions_s <- new_w_regions[location == "Colon"]

s <- ggplot(new_w_regions_s, aes(V1.y, V1.x)) +
  geom_smooth(method=lm, col = "black") +
  geom_point(shape = 21) +
  theme_bw() +
  theme(axis.title.x = element_text(size=8),
        axis.title.y = element_text(size =8))+
  xlab("Abundance-weighted average of peptidases") +
  ylab("Abundance-weighted average of auxotrophies") +
  facet_grid(. ~ location, scales = "free", space="free") +
  stat_cor(method = "spearman", size = 3, cor.coef.name = "rho",p.accuracy = 0.0001)

s

##combine all to one figure
fi_all <- ggarrange(d,j,i,s,
                   labels = c("b"),
                   ncol=2, nrow= 2, common.legend = FALSE)
fi_all

# ggsave("Output/corr_pep_auxo.pdf", plot = fi_all,
#        width = 8, height = 7)

# auxo_pep_all_d <- auxo_pep_all[location == "Duodenum"]
# auxo_pep_all_j <- auxo_pep_all[location == "Jejunum"]
# auxo_pep_all_i <- auxo_pep_all[location == "Ileum"]
# auxo_pep_all_s <- auxo_pep_all[location == "Stool"]
# 
# auxo_pepmean_d <- aggregate(auxo_pep_all_d$Freq, by = list(Peptidases = auxo_pep_all_d$count_pep, location = auxo_pep_all_d$location), FUN = mean)
# 
# auxo_pepmean_j <- aggregate(auxo_pep_all_j$Freq, by = list(Peptidases = auxo_pep_all_j$count_pep, location = auxo_pep_all_j$location), FUN = mean)
# 
# auxo_pepmean_i <- aggregate(auxo_pep_all_i$Freq, by = list(Peptidases = auxo_pep_all_i$count_pep, location = auxo_pep_all_i$location), FUN = mean)
# 
# auxo_pepmean_s <- aggregate(auxo_pep_all_s$Freq, by = list(Peptidases = auxo_pep_all_s$count_pep, location = auxo_pep_all_s$location), FUN = mean)
# 
# all_pepmean <- rbind(auxo_pepmean_d, auxo_pepmean_j, auxo_pepmean_i, auxo_pepmean_s)
# 
# #View(auxo_pepmean)
# ggplot(all_pepmean, aes(Peptidases, x, fill = location)) +
#   geom_bar(stat = "identity") +
#   facet_grid(. ~ location, scales = "free", space="free")+
#   theme_bw() +
#   theme(legend.position = "none") +
#   xlab("Number of Peptidases") +
#   ylab("Average Frequency of MAGs") +
#   xlim(1,30)
#   

# ### get substrate specifity
# substrate <-fread("/home/svenja/workspace/Merops/meropsweb121/Substrate_search.csv")
# substrate
# 
# subs_pep <- merge(domain_results_p, substrate, by.x="code", by.y="code")
# subs_pep
# 
# ##filter for substrates of bacteria
# unique(subs_pep$organism)
# filt_organims <- c("Homo sapiens", "Bos taurus", "Sus scrofa" ,"Vespula lewisii", "Mus musculus", 
#                    "Saccharomyces cerevisiae","Schistosoma japonicum","Dictyostelium discoideum",
#                    "Spirulina platensis", "Kluyveromyces lactis", "Drosophila melanogaster","Prochloron didemni")
# subs_pep_filt <- subs_pep[! subs_pep$organism %in% filt_organims]
# #View(subs_pep_filt)

###count avergae number of genoems with peptidases along gi-tract
# jej <- reimagine[location == "Jejunum",]
# 
# new2 <- auxo_pep_all[ ,sum(count_pep*Freq), by = sample]
# sample <- unique(reimagine$sample)
# location <- unique(reimagine$location)
# k <- 1
# pep_list <- list()
# 
# 
# 
# for ( i in sample) {
#   j <- reimagine[sample == i, ]
#   row <- nrow(j)
#   l <- unique(j$location)
#   tmp_pep <- domain_results_p[sample == i, ]
#   row_pep <- nrow(tmp_pep)
#   table <- data.table(nrows = row,
#                       sample = i,
#                       nrows_pep = row_pep,
#                       location = l)
#   pep_list[[k]] <- table
#   k <- k+1
# }
# pep_numb <- rbindlist(pep_list)
# pep_numb$perc <- (pep_numb$nrows_pep / pep_numb$nrows) * 100
# 
# pep_total <- aggregate(pep_numb$perc, by = list(pep_numb$location), FUN = median)
# pep_total
# 
# ###average number of peptidases per sample and location
# new_w_regions
# pep_total <- aggregate(new_w_regions$V1.y, by = list(new_w_regions$location), FUN = median)
# pep_total
# 
# pep_abun4 <- ggplot(pep_total, aes(Group.1, x, fill = Group.1)) +
#   geom_boxplot() +
#   theme_bw() +
#   xlab("Location") +
#   ylab("Abundance-weighted average of peptidases per location") +
#   theme(legend.position = "none") +
#   scale_fill_manual(values=c("#D55E00", "#E69F00", "#56B4E9", "#009E73")) +
#   theme(axis.text.x = element_text(colour = "black"),
#         axis.text.y = element_text(colour = "black"))
# 
# ggsave("Output/pep_location_abun_weigh.pdf", plot = pep_abun4,
#               width = 8, height = 7)


#########auxotrophies and specific peptiases per location

# source("Scripts/Auxotable_melted_merged.R")
# 
# 
# test4 <- merge(test, reimagine, by.x="qaccver", by.y= "Genomes")
# 
# domain_results
# gen <- unique(test4$qaccver)
# relAA <- unique(Auxotrophy_2$Compound)
# 
# p <- list()
# k <- 1
# 
# 
# for (geni in gen) {
#   print(geni) 
#   for (AAi in relAA) {
#     x <- test4[qaccver == geni]
#     y <- Auxotrophy_2[Compound == AAi]
#     z <- merge(x,y, by.x = "qaccver", by.y = "Genomes")
#     t <- z[Prototrophy == 0]
#     p[[k]] <- t
#     k <- k +1
#   }
# }
# 
# u <- rbindlist(p) 
# 
# domain <-fread("Data/domain.csv")
# domain
# 
# domain_results <- merge(domain, u, by.x="mernum", by.y="saccver")
# 
# 
# ##filter for 10 most abundant peptidases per location
# domain_results$location[domain_results$location == "Ileum"] <- "Farthest distance"
# domain_results$location[domain_results$location == "Stool"] <- "Colon"
# 
# loc <- unique(domain_results$location)
# 
# list <- list()
# k <- 1
# 
# for(loci in loc) {
#   print(loci)
#   x1 <- domain_results[location == loci]
#   sam <- unique(x1$sample)
#   for(sami in sam) {
#     x2 <- x1[sample == sami]
#     x3 <- data.table(abunpepti = names(sort(table(x2$code),decreasing=TRUE)[1:10]),
#                      location = loci,
#                      sample = sami)
#     list[[k]] <- x3
#     k<- k+1
#   }
# }
# 
# pep_abun <- rbindlist(list)
# pep_abun
# 
# pep_abun_2 <- na.omit(pep_abun)
# 
# 
# list <- list()
# k <- 1
# 
# for(loci in loc) {
#   print(loci)
#   x1 <- pep_abun_2[location == loci]
#   x2 <- data.table(abunpepti = names(sort(table(x1$abunpepti),decreasing=TRUE)[1:10]),
#                      location = loci)
#     list[[k]] <- x2
#     k<- k+1
#   }
# 
# 
# abun <- rbindlist(list)
# abun$family <- substr(abun$abunpepti, 1,1)
# 
# top10 <- ggplot(abun, aes(abunpepti, location)) +
#   geom_tile(aes(fill = location), colour = "white", lwd = 1.5, linetype = 1) +
#   theme(axis.text.x = element_text(angle = 90, colour = "black", vjust = 0.5),
#         axis.text.y = element_text(colour = "black")) +
#   theme(panel.background = element_blank(),
#         axis.line.x = element_line(colour = "black"),
#         axis.line.y = element_line(colour = "black"),
#         axis.title.y =element_blank(),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size = 12),
#         axis.text.y = element_text(size = 12)) +
#   scale_fill_manual(breaks = levels(abun$location),
#                     values = c("#000000", "#000000", "#000000",  "#000000")) +
#   labs(x = "10 most abundant peptidases") +
#   scale_y_discrete(limits = c("Colon", "Farthest distance", "Jejunum", "Duodenum")) +
#   ggh4x::facet_nested(~ "Family of Peptidases" + family, scales = "free", space = "free",nest_line = element_line(linetype = 1)) +
#   theme(strip.text.x = element_text(size = 13, colour = "black"))
# top10
# ggsave("Output/top10_pep_Reimagine.pdf", plot = top10,
#        width = 5, height =2.5)

### aufsummierte Frequenz der Peptidasen abundance-weighted avergae of peptidases per location

new3 <- merge(new2, meta, by.x="sample", by.y="sample")
new3$location[new3$location == "Ileum"] <- "Farthest distance"
new3$location[new3$location == "Stool"] <- "Colon"
##calculate pairwise wilcoxon test
kruskal.test(V1 ~ location,  data = new3)
#pairwise.wilcox.test(new3$V1, new3$location, p.adjust.method = "BH",
                     #paired = FALSE)
###visualization 
#my_comparisons <- list(c("Duodenum", "Ileum"), c("Duodenum", "Jejunum"), c("Duodenum", "Stool"),c("Ileum", "Jejunum"), c("Ileum", "Stool"), c("Jejunum", "Stool"))
new3$location_sort <- factor(new3$location, levels=c("Duodenum", "Jejunum", "Farthest distance", "Colon"))
all_pep_freq <- ggplot(new3, aes(location_sort,V1, fill = location_sort))+
  geom_boxplot(outlier.shape = NA) +
  xlab("Location") +
  ylab("Abundance-weighted average\nof #peptidases per genome") +
  scale_fill_manual(values = c("#cccccc", "#cccccc", "#cccccc", "#cccccc")) +
  geom_jitter(alpha = 0.3, width = 0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.title.x = element_text(colour = "black", size = 13,  margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(colour = "black", size = 13, margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(legend.position = "none")
all_pep_freq

# ggsave("Output/Pepfreq_Reimagine.pdf", plot = all_pep_freq,
#        width = 7, height =5)


# new2 <- auxo_pep_all[ ,sum(count_pep*Freq), by = c(sample, family)]









       