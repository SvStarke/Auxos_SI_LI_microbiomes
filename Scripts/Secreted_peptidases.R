prediction_SP <- read.delim("Data/prediction_results_Reimagine.txt")

prediction_SP$Mernum <- substr(prediction_SP$X..ID, 1, 10)

##peptidases found

test <- fread("Data/all.m8")
test

colnames(test) <- c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
                    "evalue", "bitscore")
test

test$qaccver <- substr(test$qaccver, 1, 16)

reimagine_merops <- reimagine

reimagine_merops <- merge(reimagine, test, by.x= "Genomes", by.y= "qaccver")
domain <-fread("Data/domain.csv")
domain

domain_results <- merge(domain, test, by.x="mernum", by.y="saccver")

#View(domain_results)
##delete all inhibitors

domain_results <- domain_results[type != "inhibitor",]

domain_results_p <- merge(reimagine, domain_results, by.x="Genomes", by.y = "qaccver")
domain_results_p$family <- substr(domain_results_p$code, 1,1)



secpep <- merge(domain_results_p, prediction_SP, by.x= "mernum", by.y = "Mernum")

secpep <- secpep[secpep$Prediction != "OTHER"]
secpep$location[secpep$location == "Ileum"] <- "FD"
secpep$location[secpep$location == "Stool"] <- "Colon"
secpep$location_sort <- factor(secpep$location, levels=c("Duodenum", "FD", "Colon"))

secrp <- ggplot(secpep, aes(sample, code)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.text.y = element_text(size = 11, colour = "black"),
        axis.ticks = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "white"),
        strip.text.x = element_text(colour = "black", size =11),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black")) +
  labs(x = "Samples", y = "Peptidases")  + 
  facet_grid(.~ location_sort,  scales = "free", space = "free")
secrp


ggsave("Output/Secr_pep.pdf", plot = secrp,
       width = 15, height = 5)
# secrp_new <- ggplot(secpep, aes (location, code, fill = code)) +
#   geom_col(aes(colour = code, fill = code)) +
#   guides(colour = guide_legend(title = "Peptidases")) +
#   theme(axis.text.y = element_blank(),
#         axis.ticks.y=element_blank(),
#         panel.grid =  element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         axis.text.x= element_text(colour = "black", size= 10)) +
#   scale_colour_manual(values = c("#F748A5", "#2271B2", "#359B73"),
#                       labels = c("M23.005" = "M23.005", "M23.007" = "M23.007", "M57.001")) +
#   scale_fill_manual(values = c("#F748A5", "#2271B2", "#359B73"),
#                       labels = c("M23.005" = "M23.005", "M23.007" = "M23.007", "M57.001")) +
#   xlab("Location") +
#   ylab("Secreted Peptidases") +
#   guides(colour = "none", fill=guide_legend(title="Peptidases")) 
#   
# secrp_new
# 
# 
# ggsave("Output/Secr_pep_bar_plot.pdf", plot = secrp_new,
#        width = 6, height = 3)


