###abundancies REIMAGINE study

##load scripts
source("Scripts/Mapping_ASV2HRGM.R")

source("Scripts/Predict_auxotrophies.R")

source("Scripts/Auxotable_melted_merged.R")

#######                 analyzing GROUP 1:= Dudodenum & Stool             ######

#filter for data of group 1
reimagine$location[reimagine$location == "Ileum"] <- "Farthest distance"
reimagine$location[reimagine$location == "Stool"] <- "Colon"
group1 <- reimagine[reimagine$Group1 == TRUE & Group2 == FALSE & Group3 == FALSE, ]
group1


#filter for values that are in both samples
group1_d <- group1[location == "Duodenum", ]
group1_s <- group1[location == "Colon", ]

v <- intersect(group1_d$PatNo,group1_s$PatNo)
group1 <- group1[PatNo %in% v, ]

##add Auxotrophy information

sub <- unique(group1$sample)
relAA <- unique(Auxotrophy_2$Compound)

p <- list()
k <- 1


for (subi in sub) {
  print(subi) 
  for (AAi in relAA) {
    x <- group1[sample == subi]
    y <- Auxotrophy_2[Compound == AAi]
    z <- merge(x,y, by.x = "Genomes", by.y = "Genomes")
    t <- z[Prototrophy == 0]
    p[[k]] <- t
    k <- k +1
  }
}

u <- rbindlist(p) 

sumfreq <- aggregate(u$Freq, by=list(sample=u$sample, AA=u$Compound, location =u$location), FUN=sum)
sumfreq <- as.data.table(sumfreq)

sumfreq[AA %in% c("Val","Met","Leu","Ile","Trp","Phe","Lys","His","Thr"), is.essential := "essential (human)"]
sumfreq[is.na(is.essential), is.essential := "not essential (human)"]

##wilcox test
AA <- unique(sumfreq$AA)
k <- 1
dtt <- list()

for(i in AA) {
  wilcox <- sumfreq[AA == i]
  wilcox_test <-wilcox.test(x ~ location, data = wilcox, exact = FALSE)
  dt <- data.table(p.value = wilcox_test$p.value,
                   AA = i)
  dtt[[k]] <- dt
  k <- k+1
}


wilcox <- rbindlist(dtt)
wilcox$padjust = p.adjust(wilcox$p.value, method = "fdr")
wilcox[padjust < 0.05, sign.label1 := "*"]

psignf <-wilcox[sign.label1 == "P < 0.05", ]
psignf[AA %in% c("Val","Met","Leu","Ile","Trp","Phe","Lys","His","Thr"), is.essential := "essential (human)"]
psignf[is.na(is.essential), is.essential := "not essential (human)"]

##visualization
ü <- ggplot(sumfreq[AA != "Gly"], aes(AA, x*100, fill = location)) +
  geom_boxplot(outlier.shape = NA) +
  #geom_jitter(alpha = 0.05, width = 0.2, color = "black") +
  ylab("Relative abundance of auxotrophies [%]")+
  xlab("Amino acids") +
  theme(axis.line = element_line(size=0.2, colour = "black")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.title.y = element_text(colour = "black", size = 13, margin = margin(0,10,0,0)))+
  theme(axis.title.x = element_text(colour = "black", size = 13, margin = margin(10,0,0,0))) +
  theme(axis.text.x = element_text(size=10, colour = "black", vjust = 1,  angle = 0, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size =10, colour = "black")) +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
  theme(strip.text.x = element_text(size=11, colour = "black")) +
  guides(fill = guide_legend(title = "Location")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size =11))+
  theme(legend.position = "bottom") +
  facet_grid(.~ is.essential, scales = "free_x", space = "free_x") +
  theme(strip.text.x = element_text(size=11)) +
  stat_compare_means(label = "p.signif", method = "wilcox.test", hide.ns = TRUE, symnum.args = list(cutpoints = c(0,0.05,Inf), symbols = c( "*", "ns")))
ü
##cysteine not FDR adjusted pvalue <0.05!!

# ggsave("Output/Group1.pdf", plot = ü,
#        width = 9, height = 6)

##abundance-weighted number of auxotrophies per genome
# ##load scripts
# 
# source("Scripts/Auxotable_not_melted.R")
# 
# Auxotrophy[is.na(count), count:= 0]
# 
# numb_auxos_reimagine <- merge(group1, Auxotrophy, by.x="Genomes", by.y = "Genomes")
# count_auxos_reimagine <- numb_auxos_reimagine[ ,sum(count*Freq), by = list(sample,location,PatNo)]
# count_auxos_reimagine
# 
# ##visualization
# q <- ggplot(count_auxos_reimagine,aes(location, V1, fill = location)) +
#   geom_boxplot() +
#   ylab("Abundance-weighted average\n of auxotrophies per MAG")+
#   xlab("") +
#   theme(axis.line = element_line(size=0.2, colour = "black")) +
#   theme(panel.background = element_rect(fill="white", colour= "white")) +
#   theme(axis.title.y = element_text(colour = "black", size = 12, margin = margin(0,10,0,0)))+
#   theme(axis.title.x = element_text(colour = "black", size = 12, margin = margin(10,0,0,0))) +
#   theme(axis.text.x = element_text(size=12, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
#   theme(axis.text.y = element_text(size =10, colour = "black")) +
#   theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
#   #guides(fill = guide_legend(title = "Essentiality")) +
#   scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
#   theme(legend.position = "none") 
# 
# 
# q <- q +stat_compare_means( aes(label = ..p.signif..), 
#                                   label.x = 1.5, label.y = 12)
# q
# 
# ggsave("Output/Group1_Numb_auxos.pdf", plot = q,
#        width = 7, height = 6)
# ggsave("Output/Group1_Numb_auxos.png", plot = q,
#        width = 7, height = 6)
# 
# ###number of auxotrophies per person
# auxos_onesample <- ggplot(count_auxos_reimagine, aes(PatNo, y=V1, group = location)) +
#   geom_line(aes(color=location), size =0.2) +
#   geom_point(aes(color=location), size =1) +
#   theme(legend.position = "right",
#         legend.justification = 	0.5,
#         axis.text.x = element_text(color = "black", angle = 90, vjust = 0.2, size = 10),
#         axis.text.y = element_text(color = "black", size = 10)) +
#   theme(axis.text.x = element_blank()) +
#   theme(axis.ticks.x = element_blank() ) +
#   theme(axis.title.y = element_text(size=12)) +
#   theme(axis.title.x = element_text(size=12)) +
#   theme(panel.background = element_blank()) +
#   labs(y="Abundance-weighted average of auxotrophies per MAG",x = "Samples") +
#   theme(axis.line.x = element_line(color="black", size = 0.4),
#         axis.line.y = element_line(color = "black", size =0.4)) +
#   scale_color_manual(values = c("Duodenum" = "#E69F00", "Stool" = "#56B4E9"))
# 
# ggsave("Output/Group1_auxos_onesample.pdf", plot = auxos_onesample,
#        width = 7, height = 5)
# ggsave("Output/Group1_auxos_onesample.png", plot = auxos_onesample,
#        width = 7, height = 5)
# 
#######         analyzing GROUP 2:= Dudodenum & Ileum  & Jejunum         ######
group2 <- reimagine[reimagine$Group2 == TRUE & Group1 == FALSE & Group3 == FALSE, ]
group2

group2$location

##filter vor samples from all three locations
group2_d <- group2[location == "Duodenum", ]
group2_j <- group2[location == "Jejunum", ]
group2_i <- group2[location == "Farthest distance", ]

v2 <- intersect(intersect(group2_d$PatNo, group2_j$PatNo),group2_i$PatNo)
group2 <- group2[PatNo %in% v2, ]

##add Auxotrophy information

sub <- unique(group2$sample)
relAA <- unique(Auxotrophy_2$Compound)

p <- list()
k <- 1


for (subi in sub) {
  print(subi)
  for (AAi in relAA) {
    x <- group2[sample == subi]
    y <- Auxotrophy_2[Compound == AAi]
    z <- merge(x,y, by.x = "Genomes", by.y = "Genomes")
    t <- z[Prototrophy == 0]
    p[[k]] <- t
    k <- k +1
  }
}

u2 <- rbindlist(p)

sumfreq2 <- aggregate(u2$Freq, by=list(sample=u2$sample, AA=u2$Compound, location =u2$location), FUN=sum)
sumfreq2 <- as.data.table(sumfreq2)

sumfreq2[AA %in% c("Val","Met","Leu","Ile","Trp","Phe","Lys","His","Thr"), is.essential := "essential (human)"]
sumfreq2[is.na(is.essential), is.essential := "not essential (human)"]


# ##kruskal-wallis-tets
# AA <- unique(sumfreq2$AA)
# k <- 1
# dtt <- list()
# 
# for(i in AA) {
#   Krusk_test <- sumfreq2[AA == i]
#   krusk_test <- kruskal.test(x ~ location, data = Krusk_test)
#   dt <- data.table(p.value = krusk_test$p.value,
#              AA = i)
#   dtt[[k]] <- dt
#   k <- k+1
# }
# 
# krusk2 <- rbindlist(dtt)
# krusk2$padjust = p.adjust(krusk2$p.value, method = "fdr")
# krusk2[padjust < 0.05, sign.label1 := "P < 0.05"]
# 
# 
# #boxplot
# ü2 <- ggplot(sumfreq2[AA != "Gly"], aes(AA, x*100, fill = location)) +
#   geom_boxplot(outlier.shape = NA) +
#   #geom_jitter(alpha = 0.05, width = 0.2, color = "black") +
#   ylab("Relative abundance of auxotrophies [%]")+
#   xlab("Amino acids") +
#   theme(axis.line = element_line(size=0.2, colour = "black")) +
#   theme(panel.background = element_rect(fill="white", colour= "white")) +
#   theme(axis.title.y = element_text(colour = "black", size = 10, margin = margin(0,10,0,0)))+
#   theme(axis.title.x = element_text(colour = "black", size = 10, margin = margin(10,0,0,0))) +
#   theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
#   theme(axis.text.y = element_text(size =8, colour = "black")) +
#   theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
#   guides(fill = guide_legend(title = "Location")) +
#   scale_fill_manual(values = c("#E69F00", "white", "red")) +
#   #theme(legend.position = "none") +
#   theme(legend.text = element_text(size=8)) +
#   theme(legend.title = element_text(size =10, face = "bold")) +
#   facet_grid(.~ is.essential, scales = "free_x", space = "free_x") +
#   theme(strip.text.x = element_text(size=10))
#   #stat_compare_means(method = "kruskal", label.y = 100, label = "p.signif", hide.ns=TRUE)
# ü2
# 
# ggsave("Output/Group2.pdf", plot = ü2,
#        width = 7, height = 5.5)
# ggsave("Output/Group2.png", plot = ü2,
#        width = 7, height = 5.5)
# 
# ##abundance-weighted number of auxotrophies per genome
# ##load scripts
# 
# source("Scripts/Auxotable_not_melted.R")
# 
# Auxotrophy[is.na(count), count:= 0]
# 
# numb2_auxos_reimagine <- merge(group2, Auxotrophy, by.x="Genomes", by.y = "Genomes")
# count2_auxos_reimagine <- numb2_auxos_reimagine[ ,sum(count*Freq), by = list(sample,location,PatNo)]
# count2_auxos_reimagine
# 
# 
# ##visualization
# q2 <- ggplot(count2_auxos_reimagine,aes(location, V1, fill = location)) +
#   geom_boxplot() +
#   ylab("Abundance-weighted average of auxotrophies per MAG")+
#   xlab("") +
#   theme(axis.line = element_line(size=0.2, colour = "black")) +
#   theme(panel.background = element_rect(fill="white", colour= "white")) +
#   theme(axis.title.y = element_text(colour = "black", size = 12, margin = margin(0,10,0,0)))+
#   theme(axis.title.x = element_text(colour = "black", size = 12, margin = margin(10,0,0,0))) +
#   theme(axis.text.x = element_text(size=12, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
#   theme(axis.text.y = element_text(size =10, colour = "black")) +
#   theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
#   #guides(fill = guide_legend(title = "Essentiality")) +
#   scale_fill_manual(values = c("#009E73", "#F0E442", "#0072B2")) +
#   theme(legend.position = "none") 
# q2
# 
# ggsave("Output/Group2_Numb_auxos.pdf", plot = q2,
#        width = 7, height = 6)
# ggsave("Output/Group2_Numb_auxos.png", plot = q2,
#        width = 7, height = 6)
# 
# ##visualization auxotrophies per sample
# auxos_onesample2 <- ggplot(count2_auxos_reimagine, aes(PatNo, y=V1, group = location)) +
#   geom_line(aes(color=location), size =0.2) +
#   geom_point(aes(color=location), size =1) +
#   theme(legend.position = "right",
#         legend.justification = 	0.5,
#         axis.text.x = element_text(color = "black", angle = 90, vjust = 0.2, size = 10),
#         axis.text.y = element_text(color = "black", size = 10)) +
#   theme(axis.text.x = element_blank()) +
#   theme(axis.ticks.x = element_blank() ) +
#   theme(axis.title.y = element_text(size=12)) +
#   theme(axis.title.x = element_text(size=12)) +
#   theme(panel.background = element_blank()) +
#   labs(y="Abundance-weighted average of auxotrophies per MAG",x = "Samples") +
#   theme(axis.line.x = element_line(color="black", size = 0.4),
#         axis.line.y = element_line(color = "black", size =0.4)) +
#   scale_color_manual(values = c("Duodenum" = "#009E73", "Ileum" = "#F0E442", "Jejunum" ="#0072B2"))
# 
# ggsave("Output/Group2_auxos_onesample.pdf", plot = auxos_onesample2,
#        width = 7, height = 5)
# 
#######         analyzing GROUP 3:= Dudodenum & Ileum  & Jejunum & FD     ######
group3 <- reimagine[reimagine$Group3 == TRUE, ]
group3


##filter vor samples from all three locations
group3_d <- group3[location == "Duodenum", ]
group3_j <- group3[location == "Jejunum", ]
group3_i <- group3[location == "Farthest distance", ]
group3_s <- group3[location == "Colon", ]
v3 <- intersect(intersect(group3_d$PatNo, group3_j$PatNo),intersect(group3_i$PatNo, group3$PatNo))
group3 <- group3[PatNo %in% v3, ]


##add Auxotrophy information

sub <- unique(group3$sample)
relAA <- unique(Auxotrophy_2$Compound)

p <- list()
k <- 1


for (subi in sub) {
  print(subi)
  for (AAi in relAA) {
    x <- group3[sample == subi]
    y <- Auxotrophy_2[Compound == AAi]
    z <- merge(x,y, by.x = "Genomes", by.y = "Genomes")
    t <- z[Prototrophy == 0]
    p[[k]] <- t
    k <- k +1
  }
}

u3 <- rbindlist(p)

sumfreq3 <- aggregate(u3$Freq, by=list(sample=u3$sample, AA=u3$Compound, location =u3$location, PatNo = u3$PatNo), FUN=sum)
sumfreq3 <- aggregate(u3$Freq, by=list(sample=u3$sample, AA=u3$Compound, location =u3$location, PatNo = u3$PatNo), FUN=sum)
sumfreq3 <- as.data.table(sumfreq3)

sumfreq3[AA %in% c("Val","Met","Leu","Ile","Trp","Phe","Lys","His","Thr"), is.essential := "essential (human)"]
sumfreq3[is.na(is.essential), is.essential := "not essential (human)"]
# 
# ##kruskal-wallis-tets
# AA <- unique(sumfreq3$AA)
# k <- 1
# dtt <- list()
# 
# for(i in AA) {
#   Krusk_test <- sumfreq3[AA == i]
#   krusk_test <- kruskal.test(x ~ location, data = Krusk_test)
#   dt <- data.table(p.value = krusk_test$p.value,
#                    AA = i)
#   dtt[[k]] <- dt
#   k <- k+1
# }
# 
# krusk3 <- rbindlist(dtt)
# krusk3$padjust = p.adjust(krusk3$p.value, method = "fdr")
# krusk3[padjust < 0.05, sign.label1 := "P < 0.05"]
# 
# 
# #boxplot
# ü3 <- ggplot(sumfreq3[AA != "Gly"], aes(AA, x, fill = location)) +
#   geom_boxplot(outlier.shape = NA) +
#   #geom_jitter(alpha = 0.05, width = 0.2, color = "black") +
#   ylab("Relative abundance of auxotrophies [%]")+
#   xlab("Amino acids") +
#   theme(axis.line = element_line(size=0.2, colour = "black")) +
#   theme(panel.background = element_rect(fill="white", colour= "white")) +
#   theme(axis.title.y = element_text(colour = "black", size = 10, margin = margin(0,10,0,0)))+
#   theme(axis.title.x = element_text(colour = "black", size = 10, margin = margin(10,0,0,0))) +
#   theme(axis.text.x = element_text(size=8, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
#   theme(axis.text.y = element_text(size =8, colour = "black")) +
#   theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
#   #guides(fill = guide_legend(title = "Essentiality")) +
#   scale_fill_manual(values = c("#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
#   #theme(legend.position = "none") +
#   theme(legend.text = element_text(size=8)) +
#   theme(legend.title = element_text(size =10, face = "bold"))+
#   facet_grid(.~ is.essential, scales = "free_x", space = "free_x") +
#   theme(strip.text.x = element_text(size=10))
#   #stat_compare_means(method = "kruskal", label.y = 100, label = "p.signif", hide.ns=TRUE)
# ü3
# 
# ggsave("Output/Group3.pdf", plot = ü3,
#        width = 7, height = 5.5)
# 
# ##abundance-weighted number of auxotrophies per genome
# ##load scripts
# 
# source("Scripts/Auxotable_not_melted.R")
# 
# Auxotrophy[is.na(count), count:= 0]
# 
# numb3_auxos_reimagine <- merge(group3, Auxotrophy, by.x="Genomes", by.y = "Genomes")
# count3_auxos_reimagine <- numb3_auxos_reimagine[ ,sum(count*Freq), by = list(sample,location, PatNo)]
# count3_auxos_reimagine <- data.frame(count3_auxos_reimagine)
# krusk_test <- kruskal.test(V1 ~ location, data = count3_auxos_reimagine)
# dunn_test <- pairwise.wilcox.test(count3_auxos_reimagine$V1, count3_auxos_reimagine$location)
# 
# 
# 
# ##visualization
# q3 <- ggplot(count3_auxos_reimagine,aes(location, V1, fill = location)) +
#   geom_boxplot() +
#   ylab("Abundance-weighted average of auxotrophies per MAG")+
#   xlab("") +
#   theme(axis.line = element_line(size=0.2, colour = "black")) +
#   theme(panel.background = element_rect(fill="white", colour= "white")) +
#   theme(axis.title.y = element_text(colour = "black", size = 12, margin = margin(0,10,0,0)))+
#   theme(axis.title.x = element_text(colour = "black", size = 12, margin = margin(10,0,0,0))) +
#   theme(axis.text.x = element_text(size=12, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
#   theme(axis.text.y = element_text(size =10, colour = "black")) +
#   theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
#   #guides(fill = guide_legend(title = "Essentiality")) +
#   scale_fill_manual(values = c("#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
#   theme(legend.position = "none") +
#   stat_cor(method = "kruskal")
# q3
# 
# #q3 <- q3 +stat_compare_means(method = "wilcoxan.test", label.y = 12, label = "p.signif")
# #my_comparisons <- list( c("Duodenum", "Ileum"), c("Duodenum", "Jejunum"), c("Duodenum", "Stool"),
#                        # c("Ileum", "Jejunum"), c("Ileum", "Stool"), c("Jejunum","Stool") )
# #q3 +   stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
#   #stat_compare_means(label.y = 20) 
# #q3 +   stat_pvalue_manual(r, label = "p.adj", y.position = c(8,9,10,11,12,13))
# 
# 
# 
# #### visualization individually ###
#   
# q3.1 <- ggplot(count3_auxos_reimagine,aes(location, V1, fill = location)) +
#     geom_boxplot() +
#     ylab("Abundance-weighted average of auxotrophies per MAG")+
#     xlab("") +
#     theme(axis.line = element_line(size=0.2, colour = "black")) +
#     theme(panel.background = element_rect(fill="white", colour= "white")) +
#     theme(axis.title.y = element_text(colour = "black", size = 12, margin = margin(0,10,0,0)))+
#     theme(axis.title.x = element_text(colour = "black", size = 12, margin = margin(10,0,0,0))) +
#     theme(axis.text.x = element_text(size=12, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
#     theme(axis.text.y = element_text(size =10, colour = "black")) +
#     theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
#     #guides(fill = guide_legend(title = "Essentiality")) +
#     scale_fill_manual(values = c("#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
#     theme(legend.position = "none") 
# q3.1
#   
#   
# ggsave("Output/Group3_Numb_auxos.pdf", plot = q3,
#        width = 7, height = 6)
# 
# 
# w <- pairwise.wilcox.test(count3_auxos_reimagine$V1,count3_auxos_reimagine$location, p.adjust = "fdr")
# 
# ##visualization auxotrophies per sample
# auxos_indiv <- ggplot(count3_auxos_reimagine, aes(location, y=V1, group = PatNo)) +
#   geom_line(aes(color=location), size =0.2) +
#   geom_point(aes(color=location), size =1) +
#   theme(legend.position = "none",
#         legend.justification = 	0.5,
#         axis.text.x = element_text(color = "black", angle = 0, vjust = 0.2, size = 10),
#         axis.text.y = element_text(color = "black", size = 10)) +
#   theme(axis.ticks.x = element_blank() ) +
#   theme(axis.title.y = element_text(size=12)) +
#   theme(axis.title.x = element_text(size=12)) +
#   theme(panel.background = element_blank()) +
#   labs(y="Abundance-weighted average of auxotrophies per MAG",x = "Samples") +
#   theme(axis.line.x = element_line(color="black", size = 0.4),
#         axis.line.y = element_line(color = "black", size =0.4)) +
#   scale_color_manual(values = c("Duodenum" = "#F0E442", "Farthest distance" = "#0072B2", "Jejunum" ="#D55E00", "Colon" ="#CC79A7" ))
# auxos_indiv
# 
# ggsave("Output/Group3_auxos_onesample.pdf", plot = auxos_indiv,
#        width = 7, height = 5)
# 
# #### figure for every amino acid
# 
# ###Leucine
# Leu <- sumfreq3[AA == "Leu", ]
# Leu
# Leu_stat <- Leu %>% pairwise_wilcox_test(x ~ location, adjmethod = "fdr")
# 
# 
# # Leu1 <- ggplot(Leu, aes(x = factor(location, levels = c("Duodenum", "Jejunum", "Ileum", "Stool")), x*100)) +
# #   geom_boxplot(aes(fill = location), outlier.shape = NA)  +
# #   ylab("Relative abundance of auxotrophies [%]")+
# #   xlab("") +
# #   geom_jitter(alpha = 0.2, size =1) +
# #   theme(axis.line = element_line(size=0.2, colour = "black")) +
# #   theme(panel.background = element_rect(fill="white", colour= "white")) +
# #   theme(axis.title.y = element_text(colour = "black", size = 7, margin = margin(0,10,0,0)))+
# #   theme(axis.title.x = element_text(colour = "black", size = 7, margin = margin(10,0,0,0))) +
# #   theme(axis.text.x = element_text(size=7, colour = "black", hjust = 1,  angle = 45, margin = margin(10,0,0,0))) +
# #   theme(axis.text.y = element_text(size =7, colour = "black")) +
# #   theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm")) +
# #   guides(fill = guide_legend(title = "Location")) +
# #   scale_fill_manual(values = c("#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
# #   theme(legend.position = "none") +
# #   theme(legend.text = element_text(size=7)) +
# #   theme(legend.title = element_text(size =8, face = "bold"))+
# #   facet_grid(.~ "Leucine", scales = "free_x", space = "free_x") +
# #   theme(strip.text.x = element_text(size=7)) +
# #   stat_pvalue_manual(Leu_stat, y.position = c(80), step.increase = 0.05, label = "p.adj.signif", hide.ns = TRUE, size =4)  +
# #   scale_alpha(guide = 'none')
# 
# 
# 
# Leu$PatNo <- factor(Leu$PatNo)
# Leu1 <- ggplot(Leu, aes(location, y=x, group = PatNo, fill = factor(PatNo))) +
#   geom_line(aes(colour=PatNo))+
#   geom_point(aes(color=PatNo), size =1) +
#   theme(legend.position = "none",
#         legend.justification = 	0.5,
#         axis.text.x = element_text(color = "black", angle = 0, vjust = 0.2, size = 7),
#         axis.text.y = element_text(color = "black", size = 7)) +
#   theme(axis.ticks.x = element_blank() ) +
#   theme(axis.title.y = element_text(size=7)) +
#   theme(axis.title.x = element_text(size=7)) +
#   theme(panel.background = element_blank()) +
#   labs(y="Abundance-weighted average of auxotrophies per MAG",x = "") +
#   theme(axis.line.x = element_line(color="black", size = 0.4),
#         axis.line.y = element_line(color = "black", size =0.4)) +
#   scale_colour_manual(values = c("131" = "#88CCEE", "142" = "#CC6677", "158" = "#DDCC77", "180" = "#117733", "193" = "#332288", "245" = "#AA4499", 
#                                 "273" = "#44AA99")) +
#   facet_grid(.~ "Leucine", scales = "free_x", space = "free_x") 
# 
# Leu1
# 
# ###Isoleucine
# Ile <- sumfreq3[AA == "Ile", ]
# Ile
# 
# #Ile_stat <- Ile %>% pairwise_wilcox_test(x ~ location,adjmethod = "fdr",)
# 
# Ile$PatNo <- factor(Ile$PatNo)
# Ile1 <- ggplot(Ile, aes(location, y=x, group = PatNo, fill = factor(PatNo))) +
#   geom_line(aes(colour=PatNo))+
#   geom_point(aes(color=PatNo), size =1) +
#   theme(legend.position = "none",
#         legend.justification = 	0.5,
#         axis.text.x = element_text(color = "black", angle = 0, vjust = 0.2, size = 7),
#         axis.text.y = element_text(color = "black", size = 7)) +
#   theme(axis.ticks.x = element_blank() ) +
#   theme(axis.title.y = element_text(size=7)) +
#   theme(axis.title.x = element_text(size=7)) +
#   theme(panel.background = element_blank()) +
#   labs(y="Abundance-weighted average of auxotrophies per MAG",x = "") +
#   theme(axis.line.x = element_line(color="black", size = 0.4),
#         axis.line.y = element_line(color = "black", size =0.4)) +
#   scale_colour_manual(values = c("131" = "#88CCEE", "142" = "#CC6677", "158" = "#DDCC77", "180" = "#117733", "193" = "#332288", "245" = "#AA4499", 
#                                  "273" = "#44AA99")) +
#   facet_grid(.~ "Isoleucine", scales = "free_x", space = "free_x")
# Ile1
# 
# ###Valine
# Val <- sumfreq3[AA == "Val", ]
# 
# Val$PatNo <- factor(Val$PatNo)
# Val1 <- ggplot(Val, aes(location, y=x, group = PatNo, fill = factor(PatNo))) +
#   geom_line(aes(colour=PatNo))+
#   geom_point(aes(color=PatNo), size =1) +
#   theme(legend.position = "none",
#         legend.justification = 	0.5,
#         axis.text.x = element_text(color = "black", angle = 0, vjust = 0.2, size = 7),
#         axis.text.y = element_text(color = "black", size = 7)) +
#   theme(axis.ticks.x = element_blank() ) +
#   theme(axis.title.y = element_text(size=7)) +
#   theme(axis.title.x = element_text(size=7)) +
#   theme(panel.background = element_blank()) +
#   labs(y="Abundance-weighted average of auxotrophies per MAG",x = "") +
#   theme(axis.line.x = element_line(color="black", size = 0.4),
#         axis.line.y = element_line(color = "black", size =0.4)) +
#   scale_colour_manual(values = c("131" = "#88CCEE", "142" = "#CC6677", "158" = "#DDCC77", "180" = "#117733", "193" = "#332288", "245" = "#AA4499", 
#                                  "273" = "#44AA99")) +
#   facet_grid(.~ "Valine", scales = "free_x", space = "free_x")
# Val1
# 
# ###Threonine
# Thr <- sumfreq3[AA == "Thr", ]
# Thr
# 
# 
# Thr$PatNo <- factor(Thr$PatNo)
# Thr1 <- ggplot(Thr, aes(location, y=x, group = PatNo, fill = factor(PatNo))) +
#   geom_line(aes(colour=PatNo))+
#   geom_point(aes(color=PatNo), size =1) +
#   theme(legend.position = "none",
#         legend.justification = 	0.5,
#         axis.text.x = element_text(color = "black", angle = 0, vjust = 0.2, size = 7),
#         axis.text.y = element_text(color = "black", size = 7)) +
#   theme(axis.ticks.x = element_blank() ) +
#   theme(axis.title.y = element_text(size=7)) +
#   theme(axis.title.x = element_text(size=7)) +
#   theme(panel.background = element_blank()) +
#   labs(y="Abundance-weighted average of auxotrophies per MAG",x = "") +
#   theme(axis.line.x = element_line(color="black", size = 0.4),
#         axis.line.y = element_line(color = "black", size =0.4)) +
#   scale_colour_manual(values = c("131" = "#88CCEE", "142" = "#CC6677", "158" = "#DDCC77", "180" = "#117733", "193" = "#332288", "245" = "#AA4499", 
#                                  "273" = "#44AA99")) +
#   facet_grid(.~ "Threonine", scales = "free_x", space = "free_x")
# Thr1
# 
# ###Chorismate
# Chor <- sumfreq3[AA == "Chor", ]
# Chor
# 
# 
# Chor$PatNo <- factor(Chor$PatNo)
# Chor1 <- ggplot(Chor, aes(location, y=x, group = PatNo, fill = factor(PatNo))) +
#   geom_line(aes(colour=PatNo))+
#   geom_point(aes(color=PatNo), size =1) +
#   theme(legend.position = "none",
#         legend.justification = 	0.5,
#         axis.text.x = element_text(color = "black", angle = 0, vjust = 0.2, size = 7),
#         axis.text.y = element_text(color = "black", size = 7)) +
#   theme(axis.ticks.x = element_blank() ) +
#   theme(axis.title.y = element_text(size=7)) +
#   theme(axis.title.x = element_text(size=7)) +
#   theme(panel.background = element_blank()) +
#   labs(y="Abundance-weighted average of auxotrophies per MAG",x = "") +
#   theme(axis.line.x = element_line(color="black", size = 0.4),
#         axis.line.y = element_line(color = "black", size =0.4)) +
#   scale_colour_manual(values = c("131" = "#88CCEE", "142" = "#CC6677", "158" = "#DDCC77", "180" = "#117733", "193" = "#332288", "245" = "#AA4499", 
#                                  "273" = "#44AA99")) +
#   facet_grid(.~ "Chorismate", scales = "free_x", space = "free_x")
# Chor1
# 
# ###Tryptophan
# Trp <- sumfreq3[AA == "Trp", ]
# Trp
# 
# 
# Trp$PatNo <- factor(Trp$PatNo)
# Trp1 <- ggplot(Trp, aes(location, y=x, group = PatNo, fill = factor(PatNo))) +
#   geom_line(aes(colour=PatNo))+
#   geom_point(aes(color=PatNo), size =1) +
#   theme(legend.position = "none",
#         legend.justification = 	0.5,
#         axis.text.x = element_text(color = "black", angle = 0, vjust = 0.2, size = 7),
#         axis.text.y = element_text(color = "black", size = 7)) +
#   theme(axis.ticks.x = element_blank() ) +
#   theme(axis.title.y = element_text(size=7)) +
#   theme(axis.title.x = element_text(size=7)) +
#   theme(panel.background = element_blank()) +
#   labs(y="Abundance-weighted average of auxotrophies per MAG",x = "") +
#   theme(axis.line.x = element_line(color="black", size = 0.4),
#         axis.line.y = element_line(color = "black", size =0.4)) +
#   scale_colour_manual(values = c("131" = "#88CCEE", "142" = "#CC6677", "158" = "#DDCC77", "180" = "#117733", "193" = "#332288", "245" = "#AA4499", 
#                                  "273" = "#44AA99")) +
#   facet_grid(.~ "Tryptophan", scales = "free_x", space = "free_x")
# Trp1
# 
# ###Serine
# Ser <- sumfreq3[AA == "Ser", ]
# Ser
# 
# 
# Ser$PatNo <- factor(Ser$PatNo)
# Ser1 <- ggplot(Ser, aes(location, y=x, group = PatNo, fill = factor(PatNo))) +
#   geom_line(aes(colour=PatNo))+
#   geom_point(aes(color=PatNo), size =1) +
#   theme(legend.position = "none",
#         legend.justification = 	0.5,
#         axis.text.x = element_text(color = "black", angle = 0, vjust = 0.2, size = 7),
#         axis.text.y = element_text(color = "black", size = 7)) +
#   theme(axis.ticks.x = element_blank() ) +
#   theme(axis.title.y = element_text(size=7)) +
#   theme(axis.title.x = element_text(size=7)) +
#   theme(panel.background = element_blank()) +
#   labs(y="Abundance-weighted average of auxotrophies per MAG",x = "") +
#   theme(axis.line.x = element_line(color="black", size = 0.4),
#         axis.line.y = element_line(color = "black", size =0.4)) +
#   scale_colour_manual(values = c("131" = "#88CCEE", "142" = "#CC6677", "158" = "#DDCC77", "180" = "#117733", "193" = "#332288", "245" = "#AA4499", 
#                                  "273" = "#44AA99")) +
#   facet_grid(.~ "Serine", scales = "free_x", space = "free_x")
# Ser1
# 
# ###Cysteine
# Cys <- sumfreq3[AA == "Cys"]
# Cys
# 
# 
# Cys$PatNo <- factor(Cys$PatNo)
# Cys1 <- ggplot(Cys, aes(location, y=x, group = PatNo, fill = factor(PatNo))) +
#   geom_line(aes(colour=PatNo))+
#   geom_point(aes(color=PatNo), size =1) +
#   theme(legend.position = "none",
#         legend.justification = 	0.5,
#         axis.text.x = element_text(color = "black", angle = 0, vjust = 0.2, size = 7),
#         axis.text.y = element_text(color = "black", size = 7)) +
#   theme(axis.ticks.x = element_blank() ) +
#   theme(axis.title.y = element_text(size=7)) +
#   theme(axis.title.x = element_text(size=7)) +
#   theme(panel.background = element_blank()) +
#   labs(y="Abundance-weighted average of auxotrophies per MAG",x = "") +
#   theme(axis.line.x = element_line(color="black", size = 0.4),
#         axis.line.y = element_line(color = "black", size =0.4)) +
#   scale_colour_manual(values = c("131" = "#88CCEE", "142" = "#CC6677", "158" = "#DDCC77", "180" = "#117733", "193" = "#332288", "245" = "#AA4499", 
#                                  "273" = "#44AA99")) +
#   facet_grid(.~ "Cysteine", scales = "free_x", space = "free_x")
# Cys1
# 
# ###Tyr
# Tyr <- sumfreq3[AA == "Tyr", ]
# Tyr
# 
# 
# Tyr$PatNo <- factor(Tyr$PatNo)
# Tyr1 <- ggplot(Tyr, aes(location, y=x, group = PatNo, fill = factor(PatNo))) +
#   geom_line(aes(colour=PatNo))+
#   geom_point(aes(color=PatNo), size =1) +
#   theme(legend.position = "none",
#         legend.justification = 	0.5,
#         axis.text.x = element_text(color = "black", angle = 0, vjust = 0.2, size = 7),
#         axis.text.y = element_text(color = "black", size = 7)) +
#   theme(axis.ticks.x = element_blank() ) +
#   theme(axis.title.y = element_text(size=7)) +
#   theme(axis.title.x = element_text(size=7)) +
#   theme(panel.background = element_blank()) +
#   labs(y="Abundance-weighted average of auxotrophies per MAG",x = "") +
#   theme(axis.line.x = element_line(color="black", size = 0.4),
#         axis.line.y = element_line(color = "black", size =0.4)) +
#   scale_colour_manual(values = c("131" = "#88CCEE", "142" = "#CC6677", "158" = "#DDCC77", "180" = "#117733", "193" = "#332288", "245" = "#AA4499", 
#                                  "273" = "#44AA99")) +
#   facet_grid(.~ "Tyrosine", scales = "free_x", space = "free_x")
# Tyr1
# 
# 
# ###Phe
# Phe <- sumfreq3[AA == "Phe", ]
# Phe
# 
# 
# Phe$PatNo <- factor(Phe$PatNo)
# Phe1 <- ggplot(Phe, aes(location, y=x, group = PatNo, fill = factor(PatNo))) +
#   geom_line(aes(colour=PatNo))+
#   geom_point(aes(color=PatNo), size =1) +
#   theme(legend.position = "none",
#         legend.justification = 	0.5,
#         axis.text.x = element_text(color = "black", angle = 0, vjust = 0.2, size = 7),
#         axis.text.y = element_text(color = "black", size = 7)) +
#   theme(axis.ticks.x = element_blank() ) +
#   theme(axis.title.y = element_text(size=7)) +
#   theme(axis.title.x = element_text(size=7)) +
#   theme(panel.background = element_blank()) +
#   labs(y="Abundance-weighted average of auxotrophies per MAG",x = "") +
#   theme(axis.line.x = element_line(color="black", size = 0.4),
#         axis.line.y = element_line(color = "black", size =0.4)) +
#   scale_colour_manual(values = c("131" = "#88CCEE", "142" = "#CC6677", "158" = "#DDCC77", "180" = "#117733", "193" = "#332288", "245" = "#AA4499", 
#                                  "273" = "#44AA99")) +
#   facet_grid(.~ "Phenylalanine", scales = "free_x", space = "free_x")
# Phe1
# 
# ###Met
# Met <- sumfreq3[AA == "Met", ]
# Met
# 
# 
# Met$PatNo <- factor(Met$PatNo)
# Met1 <- ggplot(Met, aes(location, y=x, group = PatNo, fill = factor(PatNo))) +
#   geom_line(aes(colour=PatNo))+
#   geom_point(aes(color=PatNo), size =1) +
#   theme(legend.position = "none",
#         legend.justification = 	0.5,
#         axis.text.x = element_text(color = "black", angle = 0, vjust = 0.2, size = 7),
#         axis.text.y = element_text(color = "black", size = 7)) +
#   theme(axis.ticks.x = element_blank() ) +
#   theme(axis.title.y = element_text(size=7)) +
#   theme(axis.title.x = element_text(size=7)) +
#   theme(panel.background = element_blank()) +
#   labs(y="Abundance-weighted average of auxotrophies per MAG",x = "") +
#   theme(axis.line.x = element_line(color="black", size = 0.4),
#         axis.line.y = element_line(color = "black", size =0.4)) +
#   scale_colour_manual(values = c("131" = "#88CCEE", "142" = "#CC6677", "158" = "#DDCC77", "180" = "#117733", "193" = "#332288", "245" = "#AA4499", 
#                                  "273" = "#44AA99")) +
#   facet_grid(.~ "Methionine", scales = "free_x", space = "free_x")
# Met1
# 
# ###Lys
# Lys <- sumfreq3[AA == "Lys", ]
# Lys
# 
# Lys$PatNo <- factor(Lys$PatNo)
# Lys1 <- ggplot(Lys, aes(location, y=x, group = PatNo, fill = factor(PatNo))) +
#   geom_line(aes(colour=PatNo))+
#   geom_point(aes(color=PatNo), size =1) +
#   theme(legend.position = "none",
#         legend.justification = 	0.5,
#         axis.text.x = element_text(color = "black", angle = 0, vjust = 0.2, size = 7),
#         axis.text.y = element_text(color = "black", size = 7)) +
#   theme(axis.ticks.x = element_blank() ) +
#   theme(axis.title.y = element_text(size=7)) +
#   theme(axis.title.x = element_text(size=7)) +
#   theme(panel.background = element_blank()) +
#   labs(y="Abundance-weighted average of auxotrophies per MAG",x = "") +
#   theme(axis.line.x = element_line(color="black", size = 0.4),
#         axis.line.y = element_line(color = "black", size =0.4)) +
#   scale_colour_manual(values = c("131" = "#88CCEE", "142" = "#CC6677", "158" = "#DDCC77", "180" = "#117733", "193" = "#332288", "245" = "#AA4499", 
#                                  "273" = "#44AA99")) +
#   facet_grid(.~ "Lysine", scales = "free_x", space = "free_x")
# Lys1
# 
# ##!!!
# ####ANMERKUNG: In der Abb. selber fehlen für 4 Patienten die Auxotrophie-Frquenzen, 
# ############### was damit zusammenhängt, dass die Protophie-Frequenz bei >98,8% ist 
# ##############  bisher noch nicht die Auxotorphie-Frequenz händisch hinzugefügt 
# 
# ###Arg
# Arg <- sumfreq3[AA == "Arg",]
# Arg
# 
# Arg$PatNo <- factor(Arg$PatNo)
# Arg1 <- ggplot(Arg, aes(location, y=x, group = PatNo, fill = factor(PatNo))) +
#   geom_line(aes(colour=PatNo))+
#   geom_point(aes(color=PatNo), size =1) +
#   theme(legend.position = "none",
#         legend.justification = 	0.5,
#         axis.text.x = element_text(color = "black", angle = 0, vjust = 0.2, size = 7),
#         axis.text.y = element_text(color = "black", size = 7)) +
#   theme(axis.ticks.x = element_blank() ) +
#   theme(axis.title.y = element_text(size=7)) +
#   theme(axis.title.x = element_text(size=7)) +
#   theme(panel.background = element_blank()) +
#   labs(y="Abundance-weighted average of auxotrophies per MAG",x = "") +
#   theme(axis.line.x = element_line(color="black", size = 0.4),
#         axis.line.y = element_line(color = "black", size =0.4)) +
#   scale_colour_manual(values = c("131" = "#88CCEE", "142" = "#CC6677", "158" = "#DDCC77", "180" = "#117733", "193" = "#332288", "245" = "#AA4499", 
#                                  "273" = "#44AA99")) +
#   facet_grid(.~ "Arginine", scales = "free_x", space = "free_x")
# Arg1 
# 
# ###Pro
# Pro <- sumfreq3[AA == "Pro", ]
# Pro
# 
# Pro$PatNo <- factor(Pro$PatNo)
# Pro1 <- ggplot(Pro, aes(location, y=x, group = PatNo, fill = factor(PatNo))) +
#   geom_line(aes(colour=PatNo))+
#   geom_point(aes(color=PatNo), size =1) +
#   theme(legend.position = "none",
#         legend.justification = 	0.5,
#         axis.text.x = element_text(color = "black", angle = 0, vjust = 0.2, size = 7),
#         axis.text.y = element_text(color = "black", size = 7)) +
#   theme(axis.ticks.x = element_blank() ) +
#   theme(axis.title.y = element_text(size=7)) +
#   theme(axis.title.x = element_text(size=7)) +
#   theme(panel.background = element_blank()) +
#   labs(y="Abundance-weighted average of auxotrophies per MAG",x = "") +
#   theme(axis.line.x = element_line(color="black", size = 0.4),
#         axis.line.y = element_line(color = "black", size =0.4)) +
#   scale_colour_manual(values = c("131" = "#88CCEE", "142" = "#CC6677", "158" = "#DDCC77", "180" = "#117733", "193" = "#332288", "245" = "#AA4499", 
#                                  "273" = "#44AA99")) +
#   facet_grid(.~ "Proline", scales = "free_x", space = "free_x")
# Pro1
# 
# ###His
# His <- sumfreq3[AA == "His",]
# His
# 
# His$PatNo <- factor(His$PatNo)
# His1 <- ggplot(His, aes(location, y=x, group = PatNo, fill = factor(PatNo))) +
#   geom_line(aes(colour=PatNo))+
#   geom_point(aes(color=PatNo), size =1) +
#   theme(legend.position = "none",
#         legend.justification = 	0.5,
#         axis.text.x = element_text(color = "black", angle = 0, vjust = 0.2, size = 7),
#         axis.text.y = element_text(color = "black", size = 7)) +
#   theme(axis.ticks.x = element_blank() ) +
#   theme(axis.title.y = element_text(size=7)) +
#   theme(axis.title.x = element_text(size=7)) +
#   theme(panel.background = element_blank()) +
#   labs(y="Abundance-weighted average of auxotrophies per MAG",x = "") +
#   theme(axis.line.x = element_line(color="black", size = 0.4),
#         axis.line.y = element_line(color = "black", size =0.4)) +
#   scale_colour_manual(values = c("131" = "#88CCEE", "142" = "#CC6677", "158" = "#DDCC77", "180" = "#117733", "193" = "#332288", "245" = "#AA4499", 
#                                  "273" = "#44AA99")) +
#   facet_grid(.~ "Histidine", scales = "free_x", space = "free_x")
# 
# His1
# 
# #### combine all figures together
# All_AA <- ggarrange(Leu1,Ile1, Val1, Thr1, Chor1, Trp1, Ser1,Cys1, Tyr1, Phe1, Met1, Lys1, Arg1, Pro1, His1,
#                     ncol = 3, nrow = 5)
# 
# ggsave("Output/AA_Group3.pdf", plot = All_AA,
#        width = 15, height = 17)

########    fold change figure for all three groups #############


###log2 FC for group 1 & group 2
#getting number of individuals for each group 
#group1
n1 <- length(v)
n2 <- length(v2)
n3 <- length(v3)
##loop for fold change group 1
relAA <- unique(sumfreq$AA)

FC.log <- list()
k <- 1
for(i in relAA) {
  FC_i <- sumfreq[AA == i]
  mean_D <- mean(FC_i$x[FC_i$location == "Duodenum"])
  mean_S <- mean(FC_i$x[FC_i$location == "Colon"])
  FC <- log10(mean_D/mean_S)
  wilcox_all <- wilcox_test(x ~ location, data = FC_i, p.adjust.method = "fdr")
  table <- data.table(AA = i,
                      FC.log = FC,
                      p.value = wilcox_all$p,
                      location = "Duodenum/Colon",
                      group = "Group 1\n(n=39)")
  FC.log[[k]] <- table
  k <- k +1
}

FC_group1 <- rbindlist(FC.log)
FC_group1

FC_group1[p.value < 0.05, sign.label1 := "Padj < 0.05"]

## loop for fold change of group 2
relAA <- unique(sumfreq2$AA)
FC.log2 <- list()
k <- 1
for(i in relAA) {
  FC_i2 <- sumfreq2[AA == i]
  mean_D <- mean(FC_i2$x[FC_i2$location == "Duodenum"])
  mean_J <- mean(FC_i2$x[FC_i2$location == "Jejunum"])
  mean_I <- mean(FC_i2$x[FC_i2$location == "Farthest distance"])
  FC_DI <- log10(mean_D/mean_I)
  FC_DJ <- log10(mean_D/mean_J)
  FC_IJ <- log10(mean_I/mean_J)
  wilcox_all <- pairwise.wilcox.test(FC_i2$x, FC_i2$location, p.adjust.method = "fdr")
  table <- data.table(AA = i,
                      FC.log = c(FC_DI,FC_DJ,FC_IJ),
                      p.value = c(wilcox_all$p.value[1, "Duodenum" ],wilcox_all$p.value[2, "Duodenum" ], wilcox_all$p.value[2, "Farthest distance" ]),
                      location = c("Duodenum/Farthest distance", "Duodenum/Jejunum", "Farthest distance/Jejunum"),
                      group = "Group 2\n(n= 14)")
  FC.log2[[k]] <- table
  k <- k +1
}

wilcox_all$p.value
FC_group2 <- rbindlist(FC.log2)
FC_group2

FC_group2[p.value < 0.05, sign.label1 := "Padj < 0.05"]


#loop for fold change of group 3
relAA <- unique(sumfreq3$AA)
FC.log3 <- list()
k <- 1
for(i in relAA) {
  FC_i3 <- sumfreq3[AA == i]
  mean_D <- mean(FC_i3$x[FC_i3$location == "Duodenum"])
  mean_J <- mean(FC_i3$x[FC_i3$location == "Jejunum"])
  mean_I <- mean(FC_i3$x[FC_i3$location == "Farthest distance"])
  mean_S <- mean(FC_i3$x[FC_i3$location == "Colon"])
  FC_DS <- log10(mean_D/mean_S)
  FC_JS <- log10(mean_J/mean_S)
  FC_IS <- log10(mean_I/mean_S)
  wilcox_all <- pairwise.wilcox.test(FC_i3$x, FC_i3$location, p.adjust.method = "fdr")
  table <- data.table(AA = i,
                      FC.log = c(FC_DS,FC_JS,FC_IS),
                      p.value = c(wilcox_all$p.value[1, "Colon" ],wilcox_all$p.value[3, "Colon" ],wilcox_all$p.value[2, "Colon" ]),
                      location = c("Duodenum/Colon", "Jejunum/Colon", "Farthest distance/Colon"),
                      group = "Group 3\n(n= 7)")
  FC.log3[[k]] <- table
  k <- k +1
}

FC_group3 <- rbindlist(FC.log3)
FC_group3

FC_group3[p.value < 0.05, sign.label1 := "Padj < 0.05"]
library(forcats)
##merge files together from group1 and group 2
FC_all <- rbind(FC_group1, FC_group3)

# create a vector with letters in the desired order
x <- c("Duodenum/Colon", "Duodenum/Jejunum", "Duodenum/Farthest distance", "Farthest distance/Jejunum","Duodenum/Colon","Jejunum/Colon","Farthest distance/Colon")



FC_all$location <- factor(FC_all$location, levels = c("Duodenum/Colon","Jejunum/Colon","Farthest distance/Colon"))

#df3 <- df2[order(df2$Group.1), ]

##visualization
FC_heatmap <- ggplot(FC_all[AA != c("Gly")], aes( x=AA, y = factor(location, levels = c("Farthest distance/Colon", "Jejunum/Colon", "Duodenum/Colon")), fill = FC.log)) +
  geom_tile() +
  geom_point(aes(shape = sign.label1), size = 0.5) +
  scale_fill_gradient2(high = "#ca0020", mid = "white", low = "#0571b0") +
  scale_shape_manual(values = 8, na.translate = FALSE) +
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(x = "Auxotrophy", y = "", shape = "",
       fill = expression(log[10]~'(Fold Change)')) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.justification = 	1,
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1, size =10),
        axis.text.y = element_text(color = "black", size =10),
        axis.title.x = element_text(size=11),
        legend.title = element_text(size =11),
        legend.text = element_text(size=8),
        strip.text.y = element_text(size = 10)) +
  facet_grid(group ~. , scales = "free_y", space = "free_y") +
  theme(strip.text.y = element_text(size=8))



FC_heatmap 


# ggsave("Output/Heatmap_FC.pdf", plot = FC_heatmap,
#        width = 6, height = 5)

###only fold change for small intestine

FC_all$location <- factor(FC_all$location, levels = c("Duodenum/Jejunum","Duodenum/Farthest distance","Farthest distance/Jejunum"))

FC_heatmap_2 <- ggplot(FC_group2[AA != c("Gly")], aes( x=AA, y = factor(location, levels = c("Farthest distance/Jejunum","Duodenum/Farthest distance","Duodenum/Jejunum")), fill = FC.log)) +
  geom_tile() +
  geom_point(aes(shape = sign.label1), size = 0.5) +
  scale_fill_gradient2(high = "#ca0020", mid = "white", low = "#0571b0") +
  scale_shape_manual(values = 8, na.translate = FALSE) +
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(x = "Auxotrophy", y = "", shape = "",
       fill = expression(log[10]~'(Fold Change)')) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.justification = 	1,
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1, size =10),
        axis.text.y = element_text(color = "black", size =10),
        axis.title.x = element_text(size=11),
        legend.title = element_text(size =11),
        legend.text = element_text(size=8),
        strip.text.y = element_text(size = 10)) +
  facet_grid(group ~. , scales = "free_y", space = "free_y") +
  theme(strip.text.y = element_text(size=8))
FC_heatmap_2

# #ggsave("Output/Heatmap_FC_small_intestine.pdf", plot = FC_heatmap_2,
#        width = 6, height = 3)