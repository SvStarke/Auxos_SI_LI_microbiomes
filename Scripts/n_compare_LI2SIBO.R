#-------------------------------------------------------------------------------
# Peptidases
#-------------------------------------------------------------------------------
pepd_profile_sibo <- t(abun_sibo_rel) %*% t(mero_mat_bin[,rownames(abun_sibo_rel)])
pepd_profile_reim <- t(abun_reim_rel) %*% t(mero_mat_bin[,rownames(abun_reim_rel)])

pepd_profile_both <- rbind(pepd_profile_sibo, pepd_profile_reim)
# remove zero-only pepdidases
pepd_profile_both <- pepd_profile_both[,colSums(pepd_profile_both) > 0]

# PCoA
pepDist <- vegdist(pepd_profile_both, method = "bray")
pca <- prcomp(pepd_profile_both)
dt_pcoa <- cmdscale(pepDist)
dt_pcoa <- data.table(sample = rownames(dt_pcoa),
                      `PCoA dimension 1` = dt_pcoa[,1],
                      `PCoA dimension 2` = dt_pcoa[,2],
                      PC1 = pca$x[,1],
                      PC2 = pca$x[,2])
dt_pcoa[, study := ifelse(sample %in% rownames(pepd_profile_sibo), "SIBO", "Reimagine")]
dt_pcoa <- merge(dt_pcoa, db_meta_reim, by = "sample", sort = FALSE, all.x = TRUE)
dt_pcoa <- merge(dt_pcoa, db_meta_sibo, by = "sample", sort = FALSE, all.x = TRUE)

dt_pcoa <- dt_pcoa[(study == "Reimagine" & location.main == "Large intestine") | (study == "SIBO")]
dt_pcoa[, tmpgrp := SIBO]
dt_pcoa[is.na(tmpgrp), tmpgrp := "Large intestine"]
dt_pcoa[tmpgrp == "SIBO", tmpgrp := "Small intestine (SIBO)"]
dt_pcoa[tmpgrp == "non-SIBO", tmpgrp := "Small intestine (non-SIBO)"]


ggplot(dt_pcoa, aes(`PCoA dimension 1`,`PCoA dimension 2`, fill = tmpgrp)) +
  geom_point(shape = 21) +
  labs(fill = "Microbiome") +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank())

dtbc <- data.table(as.table(as.matrix(pepDist)))
setnames(dtbc, c("sampleA","sampleB","BC"))
dtbc[, studyA:= ifelse(sampleA %in% rownames(pepd_profile_sibo), "SIBO", "Reimagine")]
dtbc[, studyB:= ifelse(sampleB %in% rownames(pepd_profile_sibo), "SIBO", "Reimagine")]
dtbc <- dtbc[studyA == "Reimagine" & studyB != "Reimagine"]
dtbc <- dtbc[sampleA %in% db_meta_reim[location == "LI (Stool)", sample]]
dtbc[, SIBO := ifelse(sampleB %in% db_meta_sibo[SIBO == "SIBO", sample], "SIBO","non-SIBO")]
dtbc <- dtbc[, .(mBC = median(BC)), by = .(sampleB, SIBO)]

wilcox.test(dtbc[SIBO == "SIBO", mBC],
            dtbc[SIBO == "non-SIBO", mBC])

p_lisibo_pep <- ggplot(dtbc, aes(SIBO, mBC, fill = SIBO)) +
  geom_boxplot(outlier.shape = 21) +
  labs(y = "Dissimilarity to large intestine microbiomes",
       title = "Peptidase gene profiles") +
  stat_compare_means() +
  theme_bw() +
  scale_fill_manual(values = param$sibo.colors) +
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())

#-------------------------------------------------------------------------------
# Auxoterophies
#-------------------------------------------------------------------------------
tmpaux <- auxos[rownames(abun_reim_rel),]
auxprof_reim <- t(abun_reim_rel) %*% tmpaux
tmpaux <- auxos[rownames(abun_sibo_rel),]
auxprof_sibo <- t(abun_sibo_rel) %*% tmpaux

auxprof_both <- rbind(auxprof_sibo, auxprof_reim)
auxprof_both <- auxprof_both[,colnames(auxprof_both) %notin% c("Ala","Glu","Gly","Asp")]
auxprof_both <- auxprof_both / rowSums(auxprof_both)

# PCoA
auxDist <- vegdist(auxprof_both, method = "bray")
pca <- prcomp(auxprof_both)
dt_pcoa <- cmdscale(auxDist)
dt_pcoa <- data.table(sample = rownames(dt_pcoa),
                      `PCoA dimension 1` = dt_pcoa[,1],
                      `PCoA dimension 2` = dt_pcoa[,2],
                      PC1 = pca$x[,1],
                      PC2 = pca$x[,2])
dt_pcoa[, study := ifelse(sample %in% rownames(auxprof_sibo), "SIBO", "Reimagine")]
dt_pcoa <- merge(dt_pcoa, db_meta_reim, by = "sample", sort = FALSE, all.x = TRUE)
dt_pcoa <- merge(dt_pcoa, db_meta_sibo, by = "sample", sort = FALSE, all.x = TRUE)

dt_pcoa <- dt_pcoa[(study == "Reimagine" & location.main == "Large intestine") | (study == "SIBO")]
dt_pcoa[, tmpgrp := SIBO]
dt_pcoa[is.na(tmpgrp), tmpgrp := "Large intestine"]
dt_pcoa[tmpgrp == "SIBO", tmpgrp := "Small intestine (SIBO)"]
dt_pcoa[tmpgrp == "non-SIBO", tmpgrp := "Small intestine (non-SIBO)"]

ggplot(dt_pcoa, aes(`PCoA dimension 1`,`PCoA dimension 2`, fill = tmpgrp)) +
  geom_point(shape = 21) +
  labs(fill = "Microbiome") +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank())

ggplot(dt_pcoa, aes(`PC1`,`PC2`, fill = tmpgrp)) +
  geom_point(shape = 21) +
  labs(fill = "Microbiome") +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank())


dtbc <- data.table(as.table(as.matrix(auxDist)))
setnames(dtbc, c("sampleA","sampleB","BC"))
dtbc[, studyA:= ifelse(sampleA %in% rownames(pepd_profile_sibo), "SIBO", "Reimagine")]
dtbc[, studyB:= ifelse(sampleB %in% rownames(pepd_profile_sibo), "SIBO", "Reimagine")]
dtbc <- dtbc[studyA == "Reimagine" & studyB != "Reimagine"]
dtbc <- dtbc[sampleA %in% db_meta_reim[location == "LI (Stool)", sample]]
dtbc[, SIBO := ifelse(sampleB %in% db_meta_sibo[SIBO == "SIBO", sample], "SIBO","non-SIBO")]
dtbc <- dtbc[, .(mBC = median(BC)), by = .(sampleB, SIBO)]

wilcox.test(dtbc[SIBO == "SIBO", mBC],
            dtbc[SIBO == "non-SIBO", mBC])


p_lisibo_aux <- ggplot(dtbc, aes(SIBO, mBC, fill = SIBO)) +
  geom_boxplot(outlier.shape = 21) +
  labs(y = "Dissimilarity to large intestine microbiomes",
       title = "Auxotrophy profiles") +
  stat_compare_means() +
  theme_bw() +
  scale_fill_manual(values = param$sibo.colors) +
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "none")


p_lisibo <- egg::ggarrange(p_lisibo_aux, p_lisibo_pep, ncol = 2,
                           labels = c("a","b"),
                           label.args = list(gp=grid::gpar(font=2)))

export_plot(p_lisibo, "Output/p_LI2SIBO", width = 6, height = 3.8)
