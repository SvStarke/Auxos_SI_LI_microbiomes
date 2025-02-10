#------------------------------------------------------------------------------#
# Main script for the SIBO cohort analysis
#------------------------------------------------------------------------------#

# AUXOTROPHIES

# get realtive abundance of auxotrophies
tmpaux <- auxos[rownames(abun_sibo_rel),]
tmpaux <- t(abun_sibo_rel) %*% tmpaux
auxtab_sibo <- tmpaux # for later SIBO x Reimagine analysis

dtaux <- data.table(as.table(tmpaux))
setnames(dtaux, c("sample","aa","freq"))

# merge in meta data from SIBO cohort
dtaux <- merge(db_meta_sibo, dtaux, by = "sample")
dtaux[, is.essential := ifelse(aa %in% param$essentialAA,
                               "essential (human)",
                               "not essential (human)")]

# stats
statdt <- list()
for(aai in unique(dtaux$aa)) {
  stat <- wilcox.test(freq ~ SIBO, data = dtaux[aa == aai])
  medf <- dtaux[aa == aai, median(freq),by=SIBO]
  statdt[[aai]] <- data.table(aa = aai,
                              median.nonSIBO = medf[SIBO == "non-SIBO", V1],
                              median.SIBO = medf[SIBO == "SIBO", V1],
                              pval = stat$p.value,
                              W = stat$statistic)
}
statdt <- rbindlist(statdt)
statdt[!is.nan(pval), padj := p.adjust(pval, method = "fdr")]
statdt[, plab := pval_labs(padj)]
statdt$SIBO = "SIBO"
statdt[, is.essential := ifelse(aa %in% param$essentialAA,
                                 "essential (human)",
                                 "not essential (human)")]

# plotting
p_siboAux1 <- ggplot(dtaux[!(aa %in% c("Glu","Ala","Gly","Asp"))], aes(aa, freq, fill = SIBO)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 0.8) +
  geom_text(data = statdt[!(aa %in% c("Glu","Ala","Gly","Asp"))],
            aes(aa, label = plab), y = 1) +
  facet_grid(.~is.essential, scales = "free", space = "free") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = param$sibo.colors) +
  theme_bw() +
  labs(x = "Amino acid", y = "Rel. abundance of auxotrophs") +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))


# Abundance weighted average of auxotrophies per genotype

nrAuxos <- apply(auxos[rownames(abun_sibo_rel),],1,sum)
weightedAuxos <- t(abun_sibo_rel) %*% nrAuxos
weightedAuxos <- data.table(sample = rownames(weightedAuxos), avgAuxos = weightedAuxos[,1])
weightedAuxos <- merge(db_meta_sibo, weightedAuxos, by = "sample")

p_siboAux2 <- ggplot(weightedAuxos, aes(SIBO, avgAuxos, fill = SIBO)) +
  geom_boxplot(outlier.shape = 21) +
  scale_fill_manual(values = param$sibo.colors) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none") +
  stat_compare_means(method = "wilcox", label.x = 1.2, label.y = 10.5,
                     aes(label = sprintf("p = %5.4f", as.numeric(..p.format..)))) +
  labs(y = "Weighted-average of auxotrophies\nper genotype")


# PEPTIDASES

pepd_profile <- t(abun_sibo_rel) %*% t(mero_mat_bin[,rownames(abun_sibo_rel)])
# remove zero-only pepdidases
pepd_profile <- pepd_profile[,colSums(pepd_profile) > 0]

# PCoA
pepDist <- vegdist(pepd_profile, method = "bray")
dt_pcoa <- cmdscale(pepDist)
dt_pcoa <- data.table(sample = rownames(dt_pcoa),
                      `PCoA dimension 1` = dt_pcoa[,1],
                      `PCoA dimension 2` = dt_pcoa[,2])
dt_pcoa <- merge(dt_pcoa, db_meta_sibo, by = "sample", sort = FALSE)

permanova <- adonis2(pepDist ~ SIBO, data = dt_pcoa)

annotations <- data.frame(
  xpos = Inf,
  ypos =  Inf,
  annotateText = paste0("<i>p</i> = ", permanova$`Pr(>F)`[1],"<br>",
                        "<i>R<sup>2</sup></i> = ", round(permanova$R2[1], digits = 3)),
  hjustvar = 1.05,
  vjustvar = 1.25,
  SIBO = "SIBO")

p_siboPep1 <- ggplot(dt_pcoa, aes(`PCoA dimension 1`,`PCoA dimension 2`, fill = SIBO)) +
  stat_ellipse(aes(col = SIBO)) +
  geom_point(shape = 21) +
  geom_richtext(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText),
                fill = NA, label.color = NA) +
  scale_fill_manual(values = param$sibo.colors) +
  scale_colour_manual(values = param$sibo.colors) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none")

p_sibo1 <- egg::ggarrange(p_siboAux1, p_siboAux2, p_siboPep1,
                          widths = c(0.575,0.125,0.3),
                          labels = c("a","b","c"),
                          label.args = list(gp = grid::gpar(font = 2,
                                                            cex = 1.2)))


export_plot(p_sibo1, "Output/p_sibo1", width = 12, height = 4.8)
fwrite(dtaux[!(aa %in% c("Glu","Ala","Gly","Asp"))],
       "Output/p_sibo1_a.tsv", sep = "\t", quote = FALSE)
fwrite(weightedAuxos, "Output/p_sibo1_b.tsv", sep = "\t", quote = FALSE)
fwrite(dt_pcoa, "Output/p_sibo1_c.tsv", sep = "\t", quote = FALSE)


# statistics on relativae abundances of pepdidases
dt_pepd <- data.table(pepd_profile)
dt_pepd$sample <- rownames(pepd_profile)
dt_pepd <- melt(dt_pepd, id.vars = "sample", value.name = "freq", variable.name = "code")
dt_pepd <- merge(dt_pepd, db_meta_sibo, by = "sample")

statpep <- list()
for(pepi in unique(dt_pepd$code)) {
  if((dt_pepd[code == pepi & SIBO == "SIBO", sum(freq>0)/.N]) >= (3/4) |
     (dt_pepd[code == pepi & SIBO == "non-SIBO", sum(freq>0)/.N]) >= (3/4)) {
    stat <- wilcox.test(freq ~ SIBO, data = dt_pepd[code == pepi])
    effs <- data.table(rstatix::wilcox_effsize(data = dt_pepd[code == pepi], formula = freq ~ SIBO))
    effs$p <- stat$p.value
    effs$SIBOmedian <- dt_pepd[code == pepi & SIBO == "SIBO", median(freq)]
    effs$nSIBOmedian <- dt_pepd[code == pepi & SIBO == "non-SIBO", median(freq)]
    effs$SIBOmean <- dt_pepd[code == pepi & SIBO == "SIBO", mean(freq)]
    effs$nSIBOmean <- dt_pepd[code == pepi & SIBO == "non-SIBO", mean(freq)]
    statpep[[pepi]] <- effs
  } else {
    statpep[[pepi]] <- data.table(p = NA_real_,
                                  `.y.` = "freq",
                                  group1 = "non-SIBO",
                                  group2 = "SIBO",
                                  effsize = NA_real_,
                                  n1 = NA_integer_, n2 = NA_integer_,
                                  conf.low = NA_real_, conf.high = NA_real_,
                                  magnitude = NA_character_)
  }
}
statpep <- rbindlist(statpep, idcol = "code", fill = TRUE)
statpep[!is.na(p), padj := p.adjust(p, method = "fdr")]
#statpep[order(padj)]
statpep[, family := substr(code,1,1)]
statpep[SIBOmedian < nSIBOmedian, effsize := -effsize]
statpep[SIBOmedian == 0 & nSIBOmedian == 0 & SIBOmean < nSIBOmean, effsize := -effsize]
statpep[, xmin := min(0,effsize), by = code]
statpep[, xmax := max(0,effsize), by = code]

p_siboPep2 <- ggplot(statpep, aes(effsize, code, col = family, fill = padj < 0.05)) +
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
  labs(x = "Signed effect size\n(SIBO vs. non-SIBO)",
       fill = "<i>p</i><sub>adj</sub> < 0.05",
       color = "Chemical mechanism<br>of catalysis")

# point plot for pepdidase presence/absence
pepd_profile_bin <- pepd_profile > 0
dt_pepbin <- data.table(pepd_profile_bin)
dt_pepbin$sample <- rownames(pepd_profile_bin)
dt_pepbin <- melt(dt_pepbin, id.vars = "sample", value.name = "isThere", variable.name = "code")
dt_pepbin[, family := substr(code,1,1)]
dt_pepbin <- merge(dt_pepbin, db_meta_sibo, by = "sample")

p_siboPep3 <- ggplot(dt_pepbin, aes(sample, code, col = family, shape = isThere)) +
  geom_point(size = 0.5) +
  facet_grid(family~SIBO, space = "free",scales = "free") +
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

p_siboPepA <-egg::ggarrange(p_siboPep3, p_siboPep2, nrow = 1, draw = FALSE,
                            widths = c(0.65, 0.35), labels = c("a","b"),
                            label.args = list(gp = grid::gpar(font = 2,
                                                              cex = 1.2)))

export_plot(p_siboPepA, "Output/n_siboPepD", width = 12, height = 7)
fwrite(dt_pepbin, "Output/n_siboPepD_a.tsv", sep = "\t", quote = FALSE)
fwrite(statpep, "Output/n_siboPepD_b.tsv", sep = "\t", quote = FALSE)

statpep_sibo <- copy(statpep)

#
# Heatmap
#

dbm <- copy(db_meta_sibo)
setkey(dbm, "sample")

# range [0-1]
rel_codes <- dt_pepbin[, sum(isThere/.N), by = .(SIBO,code)][V1 >= 0.75, unique(code)]
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
for(siboi in names(param$sibo.colors)) {
  tmpmat <- hm_data[,dbm[SIBO == siboi, sample]]
  clustmp <- hclust(dist(t(tmpmat), method = "euclidean"))
  spls_order <- c(spls_order, clustmp$labels[clustmp$order])
}

dt_hm <- data.table(as.table(hm_data))
setnames(dt_hm, c("code","sample","scaled.freq"))
dt_hm <- merge(dt_hm, dbm, by = "sample")
dt_hm[, family := substr(code,1,1)]
dt_hm$code <- factor(dt_hm$code, levels = codes_clust$labels[codes_clust$order])
dt_hm$sample <- factor(dt_hm$sample, levels = spls_order)

p_siboHMa <- ggplot(dt_hm, aes(sample, code, fill = scaled.freq, col = scaled.freq)) +
  geom_tile() +
  scale_fill_viridis_c(breaks = c(0,0.25,0.5,0.75,1),
                       labels = c("0.00\n(min)","0.25","0.50","0.75","1.00\n(max)"),
                       option = "inferno") +
  scale_color_viridis_c(breaks = c(0,0.25,0.5,0.75,1),
                        labels = c("0.00\n(min)","0.25","0.50","0.75","1.00\n(max)"),
                        option = "inferno") +
  facet_grid(.~SIBO, scales = "free_x", space = "free_x") +
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

p_siboHMb <- ggplot(dt_fam, aes(x = x, y = code, fill = family)) +
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

p_siboHM <- egg::ggarrange(p_siboHMa, p_siboHMb,
                           widths = c(0.975,0.025), draw = FALSE)

export_plot(p_siboHM, "Output/p_siboHM", width = 12, height = 6.2)
fwrite(dt_hm[,.(sample, code, scaled.freq, SIBO, family)],
       "Output/p_siboHM.tsv", sep = "\t", quote = FALSE)


# avg number of peptidases per genotype
avgpep <- mero[,.N, by = model]
avgpep <- merge(avgpep, data.table(model = rownames(abun_sibo_rel)),
                all.y = TRUE)
avgpep[is.na(N), N := 0]
setkey(avgpep, "model")

avgpep <- t(abun_sibo_rel) %*% avgpep[rownames(abun_sibo_rel),N]
avgpep <- data.table(sample = rownames(avgpep), avgNrPep = avgpep[,1])

avgpep <- merge(avgpep, db_meta_sibo, by = "sample")

p_siboAvgPep <- ggplot(avgpep, aes(SIBO, avgNrPep, fill = SIBO)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 21, width = 0.15) +
  scale_fill_manual(values = param$sibo.colors) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none") +
  labs(y = "Weighted-average of pepdidase genes\nper genotype") +
  stat_compare_means(comparisons = list(c("SIBO","non-SIBO")), label = "p.signif",
                     tip.length = 0, method='wilcox.test') +
  labs(x = "SIBO")

export_plot(p_reimAvgPep, "Output/p_siboPepdidasesAvg", width = 4.2, height = 4.8)
fwrite(avgpep[,.(sample, SIBO, avgNrPep)],
       "Output/p_siboPepdidasesAvg.tsv", sep = "\t", quote = FALSE)


# compare stats from Reimagine and SIBO
statpep_compare <- merge(
  statpep_sibo[,.(code, family, padj_sibo = padj, effsize_sibo = effsize)],
  statpep_reim[,.(code, family, padj_reim = padj, effsize_reim = effsize)],
  by = c("code","family"))
statpep_compare <- statpep_compare[padj_sibo < 0.05 & !is.na(effsize_reim) & abs(effsize_reim) > 0.1]

comptab <- statpep_compare[, table(a = effsize_sibo > 0, b = effsize_reim > 0)]
statcomp <- fisher.test(comptab)

compdt <- data.table(comptab)
setnames(compdt, c("SIBO","LI","N"))
compdt[, SIBO := ifelse(SIBO,"Increase","Decrease")]
compdt[, LI := ifelse(LI,"Colon","Duodenum")]

p_sibo_reim_compare <- ggplot(compdt, aes(SIBO, LI, size = N, label = as.character(N))) +
  geom_point() +
  scale_size_area() +
  geom_text(col = "black", size = 4, nudge_x = 0.25, nudge_y = 0.25) +
  labs(x = "SIBO association", y = "Primary location") +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))

export_plot(p_sibo_reim_compare, "Output/p_sibo_reim_comparison",
            width = 2.2, height = 1.7)
