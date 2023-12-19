########################################
############## PACKAGES ################
########################################

library(BiodiversityR)
library(ggord)
library(ggpubr)

########################################
############ READ DATA #################
########################################

filtered_df <- read.csv("filtered_df.csv", header=T, encoding='UTF-8')
final_matrix <- read.csv("final_matrix.csv", header=T, encoding='UTF-8')

########################################
############# FIGURE 1 #################
########################################

dfsocio <- filtered_df[, 2:6]
rownames(dfsocio) <- filtered_df$ID
colnames(dfsocio) <- c("Live_in", "Gender", "Age", "Education", "Child")

dfsocio[, 1:ncol(dfsocio)]<-lapply(dfsocio[,1:ncol(dfsocio)], as.factor)
levels(dfsocio$Child)<-c("No", "Yes")

levels(dfsocio$Age)<-c("36_55", "Over_55", "Under_36")

colnames(final_matrix)
dfnonsocio <- final_matrix[, -c(1:11, 58)]

colnames(dfnonsocio)
colnames(dfnonsocio) <- c("`Garden size large`", "`Garden size small`", "`Garden size medium`", "`Garden experience long`",    
                          "`Garden experience middle`", "`Garden experience newly`",  "`Garden type house`",   "`Garden type orchard`",    
                          "`Garden type kitchen`", "`Garden type flower`", "`Garden type other`", "`Garden type vineyard`",    
                          "`Mowing very often`", "`Mowing often`", "`Mowing rarely`", "`Mowing never`",
                          "`Growing crops for own use`", "`Growing crops no`", "`Growing crops both`",   
                          "`Biodiv network maybe`", "`Biodiv network no`", "`Biodiv network yes`", 
                          "`NMMc not know`", "`NMMc did not joint`", "`NMMc jointed`",   
                          "`Getting seeds from shop`", "`Getting seeds in person`",      
                          "`Getting seeds via collection`", "`Fertilisers yes nonsynthetic`",    
                          "`Fertilisers no`", "`Fertilisers yes synthetic`", "`Pesticides no`",      
                          "`Pestices yes eco bio`", "`Pesticides homemade`", "`Pesticides yes synthetic`",   
                          "`Support pollinators no`", "`Support pollinators natural habitats`", "`Support pollinators flower sources`",  
                          "`Support pollinators water sources`", "`Support pollinators artificial habitats`", 
                          "`Informaton from Internet`",  "`Information personal links`", "`Informaton tradicional media`",     
                          "`Observing several pollinator`", "`Observing many pollinators`",  "`Observing few pollinators`",   
                          "`Have a pond`", "`Having vegetables`",  "`Having herbs`",  "`Having fruit trees`",   
                          "`Having grapes`",  "`Having ornamental plants`",  "`Having evergreens`",  "`Having lawn`",    
                          "`Having undisturbed area`",  "`Motivation beautiful garden`", "`Motivation nature conservation`",     
                          "`Motivation making money`", "`Motivation production`",  "`Garden perception pastime`", 
                          "`Knowledge about plants`",  "`Knowledge about birds`", "`Knowledge about insects`",  
                          "`Leave unmown patches`",  "`Planting ornamental plants`", "`Use herbicides`", 
                          "`Opinion pollinator friendly`")  
colnames(dfsocio)

########################################
########## db-RDA analysis #############
########################################

RDA_mod<-dbrda(dfnonsocio~Live_in+Gender+Age+Education+Child, dfsocio,
               distance = "jaccard")

RDA_mod_MDS <-metaMDS(dfnonsocio)

sppscores(RDA_mod)<-dfnonsocio
RDA_mod$CCA$v
summary(RDA_mod)

RDA_mod$CCA$tot.chi
global_r2<-RsquareAdj(RDA_mod)  # Total variance explained by the RDA

anova.cca(RDA_mod, permutations = 999)
gc()
anova.cca(RDA_mod, by = "axis")   # Test which axis are significant
anova.cca(RDA_mod, by = "terms")  # Test which terms are significant

# adding species points to dbrda, which does not have them
# https://stackoverflow.com/questions/46531969/vegan-dbrda-species-scores-are-empty-despite-community-matrix-provided
source("specscores_dbrda.R")
sp.scores <- specscores.dbrda(RDA_mod,dfnonsocio)
spp_scores<-as.data.frame(scores(sp.scores,display="species"))
spp_scores<-spp_scores[spp_scores$dbRDA1 < -0.2 | spp_scores$dbRDA1 > 0.2 &
                         spp_scores$dbRDA2 < -0.2 | spp_scores$dbRDA2 > 0.2,]

########################################
########### Visualization ##############
########################################

windows()
dbRDA_plot_soc_all <-ggord(RDA_mod, dfsocio$Live_in,
                           polylntyp = "dashed", poly = F, 
                           alpha_el=0.2, veclsz = 1, veccol = "red", 
                           arrow=0.3, size=0.5,
                           addpts = spp_scores,
                           addsize = 3,
                           labcol = "red",
                           grp_title = "Live in", coord_fix = F,
                           ptslab = T, repel = T)+
  theme(text = element_text(size=20))

ggsave("dbRDA_plot_soc_all.svg", dbRDA_plot_soc_all, width = 15, height = 10, device = "svg")
dev.off()