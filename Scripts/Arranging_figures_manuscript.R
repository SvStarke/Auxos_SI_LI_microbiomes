###arrange figures for publication 

####figure about auxotrophy frequencies along the GI-tract

# source("Scripts/abundancies_Auxos_Gi-tract.R")
# 
# fig_auxos_GI <- ggarrange(ü, FC_heatmap,
#                       labels = c("a", "b"),
#                       ncol = 2, nrow = 1, widths = c(1,1))
# fig_auxos_GI
# 
# ggsave("Output/Fig_auxos_Gi_tract.pdf", plot = fig_auxos_GI,
#        width = 15, height = 6)
# ggsave("Output/Fig_auxos_Gi_tract.png", plot = fig_auxos_GI,
#        width = 15, height = 6)


# ###figure about proteolytic activity along the GI-tract
# 
# source("Scripts/abundancies_Auxos_Gi-tract.R")
# 
# source("Scripts/Proteolytic_activity.R")
# 
# fig_peptidases_GI <- ggarrange(bray_curtis_plot, per_pep_loc,
#                           labels = c("b", "c"),
#                           ncol = 2, nrow = 1, widths = c(1.5,1))
# 
# fig_peptidases_GI2 <- ggarrange(t_new, fig_peptidases_GI,
#                                labels = c("a"),
#                                ncol = 1, nrow = 2, heights = c(1.3,1))
# fig_peptidases_GI2
# 
# ggsave("Output/Fig_peptidases_Gi_tract.pdf", plot = fig_peptidases_GI2,
#        width = 16, height = 14)
# ggsave("Output/Fig_peptidases_Gi_tract.png", plot = fig_peptidases_GI2,
#        width = 16, height = 14)

###Alternative to combine figure 1 and 2

source("Scripts/abundancies_Auxos_Gi-tract.R")

source("Scripts/Proteolytic_activity.R")


fi_pep_auxos_tmp <- ggarrange(ü, bray_curtis_plot,
                          labels = c("a", "b"),
                           ncol = 2, nrow = 1, widths = c(1.7,1))

fi_pep_auxos <- ggarrange(fi_pep_auxos_tmp, t_new,
                                labels = c("","c"),
                                ncol = 1, nrow = 2, heights = c(1,1))

ggsave("Output/Fig_pept_auxos_Combo.pdf", plot = fi_pep_auxos,
       width = 17.5, height = 14)


###figures about SIBO
source("Scripts/SIBO_auxos.R")

source("Scripts/SIBO_Proteolytic_activity.R")

fig_SIBO_tmp <- ggarrange(ü,bray_curtis_plot,
                      labels = c("a", "b"),
                      ncol = 2, nrow = 1, widths = c(2,1.6))

fig_SIBO <- ggarrange(fig_SIBO_tmp,SIBO_pep410,
                      labels = c("", "c"),
                      nrow = 2)
fig_SIBO


ggsave("Output/Fig_comb_SIBO.pdf", plot = fig_SIBO,
       width = 14, height = 11)


##figures proteolytic activity along the GI-tract

fi_proteo_gi <- ggarrange(p, fi_all,
                          labels = c("a"),
                          ncol=1, nrow=2)
fi_proteo_gi 


ggsave("Output/Fig_proteo_GI.pdf", plot = fi_proteo_gi,
       width = 10, height = 12)


##figure about auxos and peptidases in SIBO
fi_auxo_pep_SIBO <- ggarrange(abun_auxo_SIBO,abun_pept_SIBO,
                              labels = c("a","b"),
                              ncol = 2, nrow=1, widths = c(1,1))



ggsave("Output/Fig_auxo_pep_SIBO.pdf", plot = fi_auxo_pep_SIBO,
       width =7, height = 4)



##supplementary figure Fold change
source("Scripts/abundancies_Auxos_Gi-tract.R")

FC_heatmap_all <- ggarrange(FC_heatmap, FC_heatmap_2,
                            labels = c("a","b"),
                            ncol = 2, nrow = 1)

ggsave("Output/Fig_FC_heatmap_all.pdf", plot = FC_heatmap_all,
       width = 14, height = 6)

