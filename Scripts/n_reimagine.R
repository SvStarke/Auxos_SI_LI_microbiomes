#------------------------------------------------------------------------------#
# Main script for the REIMAGINE cohort analysis
#------------------------------------------------------------------------------#

# Auxotrophy distribution by GI tract location

# get realtive abundance of auxotrophies
tmpaux <- auxos[rownames(abun_reim_rel),]
tmpaux <- t(abun_reim_rel) %*% tmpaux

dtaux <- data.table(as.table(tmpaux))
setnames(dtaux, c("sample","aa","freq"))

# merge in meta data from reimagine cohort
dtaux <- merge(db_meta_reim, dtaux, by = "sample")
dtaux[, is.essential := ifelse(aa %in% param$essentialAA,
                               "essential (human)",
                               "not essential (human)")]
dtaux$location <- factor(dtaux$location, levels = names(param$git.colors))

# this is a correlation analysis between a continous variable (Auxotr. freq) and
# a ordinal categorical variable (GI tract location). Since we have only 4
# categories Kendall's Tau is recommenden (https://journals.sagepub.com/doi/10.1177/8756479308317006) 
statgit <- list()
for(aai in unique(dtaux[!(aa %in% c("Gly","Ala","Glu","Asp")), aa])) {
  #print(aai)
  stat <- cor.test(formula = ~ as.integer(location) + freq,
                   data = dtaux[aa == aai],
                   method = "kendall")
  statgit[[aai]] <- data.table(aa = aai,
                               p = stat$p.value,
                               z = stat$statistic,
                               tau = stat$estimate)
}
statgit <- rbindlist(statgit)
statgit[, padj := p.adjust(p, method = "fdr")]
statgit[, is.essential := ifelse(aa %in% param$essentialAA,
                                 "essential (human)",
                                 "not essential (human)")]
statgit[, label := format(round(padj, digits = 3), scientific=F)]
statgit[, label := ifelse(label == "0.000","<i>p</i> < 0.001",paste0("<i>p</i> = ",label))]
statgit[, label := paste0(label,"<br>","&tau; = ",round(tau, digits = 3))]
statgit[padj < 0.05, label := paste0("<b>",label,"</b>")]
statgit[, location := "SI-FD"]

gitlabels <- dtaux[!duplicated(sample),.N,by = location]
gitlabels$location <- factor(gitlabels$location, levels = names(param$git.colors))
gitlabels <- gitlabels[order(location)]
gitlabels <- gitlabels[,paste0(location,"\n(n=",N,")")]

p_reimAux1 <- ggplot(dtaux[!(aa %in% c("Glu","Ala","Gly","Asp"))], aes(x = aa,
                  y = freq, fill = location)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 0.8) +
  geom_richtext(data = statgit[!(aa %in% c("Glu","Ala","Gly","Asp"))],
            aes(aa, label = label), y = 1.075, fill = "lightgrey", label.color = "black",
            size = 3.15, label.r = unit(0, "lines")) +
  facet_grid(.~is.essential, scales = "free", space = "free") +
  scale_y_continuous(labels = scales::percent, limits = c(0,1.15), breaks = seq(0,1,by=0.25)) +
  scale_fill_manual(values = param$git.colors, labels = gitlabels) +
  theme_bw() +
  labs(x = "Amino acid", y = "Rel. abundance of auxotrophs",
       fill = "Location") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

export_plot(p_reimAux1, "Output/p_reimAuxos", width = 12, height = 4.8)
fwrite(dtaux[!(aa %in% c("Glu","Ala","Gly","Asp")), .(sample,location, aa, freq,is.essential)],
       "Output/p_reimAuxos.tsv", sep = "\t", quote = FALSE)

# average-w. auxos
nrAuxos <- apply(auxos[rownames(abun_reim_rel),],1,sum)
weightedAuxos <- t(abun_reim_rel) %*% nrAuxos
weightedAuxos <- data.table(sample = rownames(weightedAuxos), avgAuxos = weightedAuxos[,1])
weightedAuxos <- merge(db_meta_reim, weightedAuxos, by = "sample")
weightedAuxos$location <- factor(weightedAuxos$location, levels = names(param$git.colors))

my_comparisons <- lapply(apply(combn(names(param$git.colors),2),2,list),unlist)

p_reimAux2 <- ggplot(weightedAuxos, aes(location, avgAuxos, fill = location)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 21, width = 0.15) +
  scale_fill_manual(values = param$git.colors) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none") +
  # stat_compare_means(method = "wilcox", label.x = 1.2, label.y = 10.5,
  #                    aes(label = sprintf("p = %5.4f", as.numeric(..p.format..)))) +
  labs(y = "Weighted-average of auxotrophies\nper genotype") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",
                     tip.length = 0, method='wilcox.test')+
  # stat_compare_means(label.y = 20) +
  labs(x = "Location")
p_reimAux2

export_plot(p_reimAux2, "Output/p_reimAuxosAvg", width = 4.2, height = 4.8)
fwrite(weightedAuxos[,.(sample, location, avgAuxos)],
       "Output/p_reimAuxosAvg.tsv", sep = "\t", quote = FALSE)


# PEPTIDASES

pepd_profile <- t(abun_reim_rel) %*% t(mero_mat_bin[,rownames(abun_reim_rel)])
# remove zero-only pepdidases
pepd_profile <- pepd_profile[,colSums(pepd_profile) > 0]

# PCoA
pepDist <- vegdist(pepd_profile, method = "bray")
dt_pcoa <- cmdscale(pepDist)
dt_pcoa <- data.table(sample = rownames(dt_pcoa),
                      `PCoA dimension 1` = dt_pcoa[,1],
                      `PCoA dimension 2` = dt_pcoa[,2])
dt_pcoa <- merge(dt_pcoa, db_meta_reim, by = "sample", sort = FALSE)
dt_pcoa$location <- factor(dt_pcoa$location, levels = names(param$git.colors))

permanova <- adonis2(pepDist ~ location, data = dt_pcoa)

annotations <- data.frame(
  xpos = Inf,
  ypos =  Inf,
  annotateText = paste0("<i>p</i> = ", permanova$`Pr(>F)`[1],"<br>",
                        "<i>R<sup>2</sup></i> = ", round(permanova$R2[1], digits = 3)),
  hjustvar = 1.05,
  vjustvar = 1.25,
  location = "SI-Duodenum")

p_reimPep1 <- ggplot(dt_pcoa, aes(`PCoA dimension 1`,`PCoA dimension 2`, fill = location)) +
  stat_ellipse(aes(col = location)) +
  geom_point(shape = 21) +
  geom_richtext(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText),
                fill = NA, label.color = NA) +
  scale_fill_manual(values = param$git.colors) +
  scale_colour_manual(values = param$git.colors) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black")) +
  labs(fill = "Location", col = "Location")

export_plot(p_reimPep1, "Output/p_reimPepPCoA", width = 4.9, height = 3.7)
fwrite(dt_pcoa[,.(sample, location, `PCoA dimension 1`, `PCoA dimension 2`)], "Output/p_reimPepPCoA.tsv", sep = "\t", quote = FALSE)

# statistics on relative abundances of pepdidases (compare only SI-Duodenum with LI (Stool))
# here: Paired stats
dt_pepd <- data.table(pepd_profile)
dt_pepd$sample <- rownames(pepd_profile)
dt_pepd <- melt(dt_pepd, id.vars = "sample", value.name = "freq", variable.name = "code")
dt_pepd <- merge(dt_pepd, db_meta_reim, by = "sample")
dt_pepd <- dt_pepd[location %in% c("SI-Duodenum","LI (Stool)")]
patPairs <- dt_pepd[!duplicated(sample)][,.N,by= PatNo][N == 2, PatNo]
dt_pepd <- dt_pepd[PatNo %in% patPairs][order(location, PatNo)]
statpep <- list()
for(pepi in unique(dt_pepd$code)) {
  
  # Filter to those peptidases, which have >=75% prevalence in at least one group (Duodenum or Stool)
  if((dt_pepd[code == pepi & location == "LI (Stool)", sum(freq>0)/.N]) >= (3/4) |
     (dt_pepd[code == pepi & location == "SI-Duodenum", sum(freq>0)/.N]) >= (3/4)) {
    stat <- wilcox.test(dt_pepd[code == pepi & location == "LI (Stool)", freq],
                        dt_pepd[code == pepi & location == "SI-Duodenum", freq],
                        paired = TRUE)
    effs <- data.table(rstatix::wilcox_effsize(data = dt_pepd[code == pepi], formula = freq ~ location, paired = TRUE))
    effs$p <- stat$p.value
    effs$LImedian <- dt_pepd[code == pepi & location == "LI (Stool)", median(freq)]
    effs$SImedian <- dt_pepd[code == pepi & location == "SI-Duodenum", median(freq)]
    effs$LImean <- dt_pepd[code == pepi & location == "LI (Stool)", mean(freq)]
    effs$SImean <- dt_pepd[code == pepi & location == "SI-Duodenum", mean(freq)]
    statpep[[pepi]] <- effs
  } else {
    statpep[[pepi]] <- data.table(p = NA_real_,
                                  `.y.` = "freq",
                                  group1 = "LI (Stool)",
                                  group2 = "SI-Duodenum",
                                  effsize = NA_real_,
                                  n1 = NA_integer_, n2 = NA_integer_,
                                  magnitude = NA_character_)
  }
}
statpep <- rbindlist(statpep, idcol = "code", fill = TRUE)
statpep[!is.na(p), padj := p.adjust(p, method = "fdr")]
#statpep[order(padj)]
statpep[, family := substr(code,1,1)]
statpep[LImedian < SImedian, effsize := -effsize]
statpep[LImedian == 0 & SImedian == 0 & LImean < SImean, effsize := -effsize]
statpep[, xmin := min(0,effsize), by = code]
statpep[, xmax := max(0,effsize), by = code]

p_reimPep2 <- ggplot(statpep, aes(effsize, code, col = family, fill = padj < 0.05)) +
  geom_pointrange(aes(xmin = xmin, xmax = xmax),
                  shape = 21) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_grid(family~., space = "free",scales = "free") +
  scale_fill_manual(values = c("white", "black"), na.translate = F) +
  scale_color_manual(values = param$mero.familyColors,
                     labels = param$mero.families) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black"),
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.background = element_blank(),
        legend.title = element_markdown()) +
  labs(x = "Signed effect size\n(LI vs. SI-Duodenum)",
       fill = "<i>p</i><sub>adj</sub> < 0.05",
       color = "Chemical mechanism<br>of catalysis")

# point plot for pepdidase presence/absence
pepd_profile_bin <- pepd_profile > 0
dt_pepbin <- data.table(pepd_profile_bin)
dt_pepbin$sample <- rownames(pepd_profile_bin)
dt_pepbin <- melt(dt_pepbin, id.vars = "sample", value.name = "isThere", variable.name = "code")
dt_pepbin[, family := substr(code,1,1)]
dt_pepbin <- merge(dt_pepbin, db_meta_reim, by = "sample")
dt_pepbin$location <- factor(dt_pepbin$location, levels = names(param$git.colors))

p_reimPep3 <- ggplot(dt_pepbin, aes(sample, code, col = family, shape = isThere)) +
  geom_point(size = 0.5) +
  facet_grid(family~location, space = "free",scales = "free") +
  scale_shape_manual(values=c(NA,19)) +
  scale_color_manual(values = param$mero.familyColors,
                     labels = param$mero.families) +
  theme_bw() +
  theme(axis.text.y = element_text(color = "black"),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_text(face = "bold"),
        strip.background = element_blank(),
        legend.position = "none") +
  labs(shape = "Rel. abundance > 0",
       color = "Chemical mechanism\nof catalysis", 
       x = "Samples", y = "Peptidase species")


p_reimPepA <- egg::ggarrange(p_reimPep3, p_reimPep2, nrow = 1, draw = FALSE,
                             widths = c(0.65, 0.35), labels = c("a","b"),
                             label.args = list(gp = grid::gpar(font = 2,
                                                              cex = 1.2)))

export_plot(p_reimPepA, "Output/n_reimPepD", width = 12, height = 8.3)
fwrite(dt_pepbin, "Output/n_reimPepD_a.tsv", sep = "\t", quote = FALSE)
fwrite(statpep, "Output/n_reimPepD_b.tsv", sep = "\t", quote = FALSE)

statpep_reim <- copy(statpep)

#
# Heatmap
#

dbm <- copy(db_meta_reim)
setkey(dbm, "sample")

# range [0-1]
rel_codes <- dt_pepbin[, sum(isThere/.N), by = .(location,code)][V1 >= 0.75, unique(code)]
hm_data <- pepd_profile[,rel_codes]
rel_spls <- rownames(hm_data)
hm_data <- apply(hm_data, 2, function(x) (x-min(x)) / (max(x) - min(x)))
#hm_data <- apply(hm_data, 2, sqrt)
hm_data <- t(hm_data)
colnames(hm_data) <- rel_spls

dbm <- dbm[colnames(hm_data)]

# cluster codes
codes_clust <- hclust(dist(hm_data, method = "euclidean"), method = "average")

# cluster samples
spls_order <- c()
for(loci in names(param$git.colors)) {
  tmpmat <- hm_data[,dbm[location == loci, sample]]
  clustmp <- hclust(dist(t(tmpmat), method = "euclidean"))
  spls_order <- c(spls_order, clustmp$labels[clustmp$order])
}

dt_hm <- data.table(as.table(hm_data))
setnames(dt_hm, c("code","sample","scaled.freq"))
dt_hm <- merge(dt_hm, dbm, by = "sample")
dt_hm[, family := substr(code,1,1)]
dt_hm$location <- factor(dt_hm$location, levels = names(param$git.colors))
dt_hm$code <- factor(dt_hm$code, levels = codes_clust$labels[codes_clust$order])
dt_hm$sample <- factor(dt_hm$sample, levels = spls_order)

p_reimHMa <- ggplot(dt_hm, aes(sample, code, fill = scaled.freq, col = scaled.freq)) +
  geom_tile() +
  scale_fill_viridis_c(breaks = c(0,0.25,0.5,0.75,1),
                       labels = c("0.00\n(min)","0.25","0.50","0.75","1.00\n(max)"),
                       option = "inferno") +
  scale_color_viridis_c(breaks = c(0,0.25,0.5,0.75,1),
                        labels = c("0.00\n(min)","0.25","0.50","0.75","1.00\n(max)"),
                        option = "inferno") +
  facet_grid(.~location, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(axis.text = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom",
        plot.margin = unit(c(2,0,2,2), "pt")) +
  labs(fill = "Row range-scaled\nfrequency", color = "Row range-scaled\nfrequency",
       y = "Peptidase species",
       x = "Sample")

dt_fam <- data.table(code = unique(dt_hm$code))
dt_fam[, x := "A"]
dt_fam[, family := substr(code,1,1)]
dt_fam$code <- factor(dt_fam$code, levels = codes_clust$labels[codes_clust$order])

p_reimHMb <- ggplot(dt_fam, aes(x = x, y = code, fill = family)) +
  geom_tile() +
  scale_fill_manual(values = param$mero.familyColors,
                    labels = param$mero.families) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(2,2,2,0), "pt")) +
  labs(fill = "Chemical mechanism\nof catalysis")

p_reimHM <- egg::ggarrange(p_reimHMa, p_reimHMb,
                           widths = c(0.975,0.025), draw = FALSE)

export_plot(p_reimHM, "Output/p_reimHM", width = 12, height = 8)
fwrite(dt_hm[,.(sample, code, scaled.freq, location, family)],
       "Output/p_reimHM.tsv", sep = "\t", quote = FALSE)

# avg number of peptidases per genotype
avgpep <- mero[,.N, by = model]
avgpep <- merge(avgpep, data.table(model = rownames(abun_reim_rel)),
                all.y = TRUE)
avgpep[is.na(N), N := 0]
setkey(avgpep, "model")

avgpep <- t(abun_reim_rel) %*% avgpep[rownames(abun_reim_rel),N]
avgpep <- data.table(sample = rownames(avgpep), avgNrPep = avgpep[,1])

avgpep <- merge(avgpep, db_meta_reim, by = "sample")
avgpep$location <- factor(avgpep$location, levels = names(param$git.colors))

p_reimAvgPep <- ggplot(avgpep, aes(location, avgNrPep, fill = location)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 21, width = 0.15) +
  scale_fill_manual(values = param$git.colors) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none") +
  labs(y = "Weighted-average of pepdidase genes\nper genotype") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",
                     tip.length = 0, method='wilcox.test') +
  labs(x = "Location")

export_plot(p_reimAvgPep, "Output/p_reimPepdidasesAvg", width = 4.2, height = 4.8)
fwrite(avgpep[,.(sample, location, avgNrPep)],
       "Output/p_reimPepdidasesAvg.tsv", sep = "\t", quote = FALSE)

#
# Last step here: Are there correlations between auxotrophies and pepdidases?
# Detail: We'll do this within the Duodenum and Stool group (both are the largest)
#

AuxPepCor <- list()
k <- 1
for(loci in c("SI-Duodenum","LI (Stool)")) {
  spls <- db_meta_reim[location == loci, sample]
  for(aai in colnames(tmpaux)) {
    if(!(aai %in% c("Ala","Glu","Gly","Asp"))) {
      dttmpa <- dtaux[aa == aai & sample %in% spls,
                      .(sample, auxfreq = freq)]
      for(pepi in unique(statpep$code)) {
        dttmp <- merge(dt_pepd[sample %in% spls & code == pepi,
                               .(sample,pepfreq = freq)],
                       dttmpa,
                       by = "sample") 
        if(dttmp[,sum(pepfreq>0)/.N] > 3/4) {
          stat <- cor.test(dttmp$pepfreq, dttmp$auxfreq, method = "spearman")
          AuxPepCor[[k]] <- data.table(aa = aai, code = pepi,
                                       p = stat$p.value, rho = stat$estimate,
                                       location = loci)
          k <- k + 1
        }
      }
    }
  }
}
AuxPepCor <- rbindlist(AuxPepCor)
AuxPepCor[,padj := p.adjust(p, method = "fdr")]
AuxPepCor[, family := substr(code, 1,1)]
AuxPepCor$location <- factor(AuxPepCor$location, levels = c("SI-Duodenum", 
                                                            "LI (Stool)"))
AuxPepCor[, plab := ">= 0.05"]
# AuxPepCor[padj < 0.1, plab := "< 0.1"]
AuxPepCor[padj < 0.05, plab := "< 0.05"]
AuxPepCor[padj < 0.01, plab := "< 0.01"]
AuxPepCor[padj < 0.001, plab := "< 0.001"]
AuxPepCor$plab <- factor(AuxPepCor$plab, levels = c(">= 0.05","< 0.05", "< 0.01", "< 0.001"))
aakeep <- AuxPepCor[padj < 0.001, unique(aa)]
pepkeep <- AuxPepCor[padj < 0.001, unique(code)]

p_reimAuxPep <- ggplot(AuxPepCor[aa %in% aakeep & code %in% pepkeep],
                       aes(aa, code, fill = rho, shape = plab)) +
  geom_tile() +
  geom_point() +
  scale_fill_gradient2(low = "#b2182b", high = "#2166ac") +
  scale_shape_manual(values = c(NA,4,1,19)) +
  scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
  facet_grid(.~location) +
  theme_bw() +
  theme(axis.text = element_text(colour = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_markdown(),
        legend.position = "right") +
  labs(x = "Amino acid auxotrophy", y = "Pepdidase species",
       shape = "<i>p</i><sub>adj</sub>",
       fill = "Spearman's &rho;")

export_plot(p_reimAuxPep, "Output/p_reimAuxPep", width = 10, height = 8.6)
fwrite(AuxPepCor, "Output/p_reimAuxPep.tsv", sep = "\t", quote = FALSE)

