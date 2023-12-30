#### predict auxotrophies in SIBO ####

##load scripts
source("Scripts/Mapping_ASV2HRGM_SIBO.R")

source("Scripts/Predict_auxotrophies.R")

source("Scripts/Auxotable_melted_merged.R")

##add Auxotrophy information
metadata
sub <- unique(SIBO_study$sample)
relAA <- unique(Auxotrophy_2$Compound)

p <- list()
k <- 1


for (subi in sub) {
  print(subi) 
  for (AAi in relAA) {
    x <- SIBO_study[sample == subi]
    y <- Auxotrophy_2[Compound == AAi]
    z <- merge(x,y, by.x = "Genomes", by.y = "Genomes")
    t <- z[Prototrophy == 0]
    p[[k]] <- t
    k <- k +1
  }
}

u <- rbindlist(p) 

sumfreq <- aggregate(u$Freq, by=list(sample= u$sample, AA= u$Compound, SIBO =u$V41), FUN=sum)
sumfreq <- as.data.table(sumfreq)

sumfreq[AA %in% c("Val","Met","Leu","Ile","Trp","Phe","Lys","His","Thr"), is.essential := "essential (human)"]
sumfreq[is.na(is.essential), is.essential := "not essential (human)"]


##wilcox test to add pvalues manually
AA <- unique(sumfreq$AA)
k <- 1
dtt <- list()

for(i in AA) {
  wilcox_data <- sumfreq[AA == i]
  wilcox_test <- wilcox.test(x ~ SIBO, data = wilcox_data, exact = FALSE)
  dt <- data.table(p.value = wilcox_test$p.value,
                   AA = i)
  dtt[[k]] <- dt
  k <- k+1
}


wilcox <- rbindlist(dtt)
wilcox$padjust = p.adjust(wilcox$p.value, method = "fdr")
wilcox[padjust < 0.05, sign.label1 := "*"]
wilcox

###add pvalues manually

sumfreq$sign.label1 <- NA
sumfreq$sign.label1 <- as.character(sumfreq$sign.label1)

sumfreq[sumfreq$sample == "SRR12533986" & sumfreq$AA == "His", "sign.label1"] <- "*"
sumfreq[sumfreq$sample == "SRR12533986" & sumfreq$AA == "Arg", "sign.label1"] <- "*"
sumfreq[sumfreq$sample == "SRR12533986" & sumfreq$AA == "Pro", "sign.label1"] <- "*"
sumfreq[sumfreq$sample == "SRR12533986" & sumfreq$AA == "Ser", "sign.label1"] <- "*"

ü <- ggplot(sumfreq[AA != "Gly"], aes(AA, x*100, fill = SIBO)) +
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
  guides(fill = guide_legend(title = "State")) +
  scale_fill_manual(values = c("#D55E00", "#F0E422"),labels=c('non-SIBO', 'SIBO')) +
  theme(legend.text = element_text(size=8)) +
  theme(legend.title = element_text(size =11))+
  facet_grid(.~ is.essential, scales = "free_x", space = "free_x") +
  theme(strip.text.x = element_text(size=11)) +
  geom_text(aes(label = sign.label1), position = position_dodge(width = 0.1),vjust = c(-10), size = 6) 
ü
ggsave("Output/SIBO_auxos.pdf", plot = ü,
       width = 9, height = 6)

# ###get number of samples SIBO vs non-SIBO
# samples_real <- unique(u$sample)
# SIBO_real <- SIBO[SIBO$sample %in% samples_real,]
# nrow(SIBO_real[SIBO_real$V41 == "non-SIBO"])
# nrow(SIBO_real[SIBO_real$V41 == "SIBO-replicate"])

