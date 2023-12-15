########################################
############## PACKAGES ################
########################################

library(ggplot2) 
library(ggpubr)

########################################
############ READ DATA #################
########################################

final_matrix <- read.csv("final_matrix.csv", header=T, encoding='UTF-8')

########################################
############# FIGURE 2 #################
########################################

### Biodiversity-positive activities ###

biodiv_positive_colnames <- c("Pond", "Unmowpatch", 
                              "Pesti_No", "PollSup_Nhab", 
                              "PollSup_Ahab", "PollSup_Food", 
                              "PollSup_Wate", "NMowM_Join")

positive_data <- data.frame(egyes = colSums(final_matrix[, biodiv_positive_colnames]), 
                            nulla = nrow(final_matrix) - colSums(final_matrix[, biodiv_positive_colnames]),
                            percentage_of_one = colSums(final_matrix[, biodiv_positive_colnames]/(nrow(final_matrix))*100))

positive_data$names <- c("Have pond", "Leave unmown\npatches", "Do not use\npesticides",
                         "Natural habitats", "Artificial habitats", 
                         "Flower sources", "Water sources", 
                         "Jointed to\nNoMowMay\ncampaign")

positives <- ggplot(positive_data, aes(x=reorder(names, percentage_of_one),
                                       y=percentage_of_one, fill = percentage_of_one)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  ylab("Presence (%)")+
  xlab("")+
  ylim(0, 100)+
  scale_fill_continuous(low="aquamarine1", high="cyan4") +
  theme_classic()+
  theme(legend.position = "none")

### Biodiversity-negative activities ###
biodiv_negative_colnames <- c("Mow_Voft", "Ferti_Synt", "Pesti_Synt",
                              "Herbic", "Garea_Undis")

negative_data <- data.frame(egyes = colSums(final_matrix[, biodiv_negative_colnames]), 
                            nulla = nrow(final_matrix) - colSums(final_matrix[, biodiv_negative_colnames]),
                            percentage_of_one = colSums(final_matrix[, biodiv_negative_colnames]/(nrow(final_matrix))*100))

negative_data$names <- c("Mowing\nseveral times\na month", "Use synthetic\nfertilisers", 
                         "Use synthetic\npesticides",
                         "Use herbicides", "Lack of\nundisturbed\narea*")

negative_data$egyes[5] <- negative_data$nulla[5] #Since to the question whether respondents had undisturbed areas
#the ‘no’ answer was considered as a biodiversity-negative practice 
negative_data$percentage_of_one[5] <- 100-negative_data$percentage_of_one[5]

netagives <- ggplot(negative_data, aes(x=reorder(names, percentage_of_one), 
                                       y=percentage_of_one, 
                                       fill = percentage_of_one)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  ylab("Presence (%)")+
  xlab("")+
  ylim(0, 100)+
  scale_fill_continuous(low="darksalmon", high="darkred") +
  theme_classic()+
  theme(legend.position = "none")

### Visulaization ###

pos_and_neg <- ggarrange(positives + theme(text = element_text(size=15)), 
                         netagives + theme(text = element_text(size=15)), 
                         labels = c("A", "B"),
                         ncol = 2, nrow = 1)

ggsave("Figure_2.png", pos_and_neg, width = 10, height = 5)