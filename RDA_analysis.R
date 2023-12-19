########################################
############## PACKAGES ################
########################################

library(dplyr)
library(cooccur)
library(igraph) 
library(vegan) 
library(ggord)

########################################
############# FUNCTIONS ################
########################################

### Function for creating network ###

make_net <- function(dat, probability)
{
  
  # dat = co_df
  # probability = 0.6
  
  dat<-dat$results[dat$results$p_lt<0.05 | dat$results$p_gt<0.05, ]
  dat$sign<-ifelse(dat$p_lt<0.05, -1, 1)
  # association strength is calculated as the difference between the observed and expected co-occurrence *sign (to keep it always positive)
  dat$weight<-(dat$obs_cooccur-dat$exp_cooccur)*dat$sign
  
  dat<- dat[dat$prob_cooccur>=probability, ]
  
  if(nrow(dat)>0){
    nodes <- data.frame(id= unique(c(dat$sp1, dat$sp2)), name = unique(c(dat$sp1_name, dat$sp2_name)),
                        color = "#606482")
    
    edges <- data.frame(from = dat$sp1, to = dat$sp2, obs_cooccur = dat$obs_cooccur,
                        prob_cooccur = dat$prob_cooccur, exp_cooccur = dat$exp_cooccur,
                        weight = dat$weight, sign = dat$sign,
                        color = ifelse(dat$p_lt <= 0.05, "blue", adjustcolor("coral2", alpha.f = .3)))
    
    
    
    igr_net<- graph_from_data_frame(d=edges, vertices=nodes, directed=F)
    igr_net}
}

########################################

### all network generation into a function ###

network_gen<-function(dat, cutoffs = c(0.05, 0.95), prob_cutoff = 0.6){
  
  # dat = filt_dat
  cutoff_low<-cutoffs[1]
  cutoff_high<-cutoffs[2]
  
  #####################################################
  ### EXCLUDE cols with too few or many responses  ####
  
  final_matrix_f = dat[,colSums(dat) > nrow(dat)*cutoff_low]
  cat(paste("Lower cutoff removed", 
            paste(colnames(dat[,!colSums(dat) > nrow(dat)*cutoff_low]), collapse = ", "),
            "columns.\n"))
  
  final_matrix_ff = final_matrix_f[,colSums(final_matrix_f) < nrow(final_matrix_f)*cutoff_high]
  cat(paste("Higher cutoff removed", 
            paste(colnames(final_matrix_f[,!colSums(final_matrix_f) < 
                                            nrow(final_matrix_f)*cutoff_high]), collapse = ", "),
            "columns.\n"))
  
  #####################################################
  
  co <- cooccur(t(final_matrix_ff), spp_names = TRUE)
  
  net<-make_net(co, prob_cutoff)
  
  if(!is.null(net))
  {
    # detects modules and adds the full module output to the network and 
    # membership (which node belongs to which module) to the node information
    net$mod<-cluster_louvain(net)
    V(net)$mod<-cluster_louvain(net)$membership
    
    
    # repeat the same process but module detection will be run on positive edges only
    
    # removing negative edges
    g=delete.edges(net, which(E(net)$sign ==-1))
    
    net$mod_poz<-cluster_louvain(g)
    V(net)$mod_poz<-cluster_louvain(g)$membership
    net}
}

########################################

########################################
############ READ DATA #################
########################################

filtered_df <- read.csv("filtered_df.csv", header=T, encoding='UTF-8')
final_matrix <- read.csv("final_matrix.csv", header=T, encoding='UTF-8')

new_colnames_for_filterd_df <- c("NUTS3", "Live_in", "Sex", "Age", 
                                 "Education", "Children", "Garden_size", "Pond", "G_experience", 
                                 "Garden_type", "Garea_Vege", "Garea_Herb", "Garea_Frui", 
                                 "Garea_Grap", "Garea_Orna", "Garea_Egree", "Garea_Lawn", "Garea_Undis",
                                 "Habi_Beau",  "Habi_Cons", "Habi_Mone", "Habi_Produ", 
                                 "Moti_Pastime", "Knowledge_Plant", "Knowledge_Bird", "Knowledge_Insect",  
                                 "Mowing", "Unmowpatch", "Plat_Oplants", "Get_Seeds_Plants",  
                                 "Fertiliser", "Herbic", "Growing",  "Pesticide", "Observing_polli", 
                                 "Support_Polli", "Polli_Friendly", "Biodiw_Network", "NOMowMay",
                                 "Information", "ID")

########################################
############# FIGURE 5 #################
########################################

potential_links<-as.data.frame(combn(new_colnames_for_filterd_df[-c(1:6,41)], 2))
potential_link_names<-c(apply(potential_links, 2, function(x){paste(x[1], x[2], sep="+")}))

expl_dat_df<-filtered_df[, c(2:6)]
colnames(expl_dat_df)<-c("Live_in", "Gender", "Age", "Education", "Children")

expl_dat_df[, 1:ncol(expl_dat_df)]<-lapply(expl_dat_df[,1:ncol(expl_dat_df)], as.factor)
levels(expl_dat_df$Children)<-c("No", "Yes")   
levels(expl_dat_df$Age)<-c("36_55", "Over_55", "Below_36")

sapply(expl_dat_df, table)                

combinations <- expand.grid(Live_in = levels(expl_dat_df$Live_in), 
                            Gender = levels(expl_dat_df$Gender), 
                            Age = levels(expl_dat_df$Age),
                            Education = levels(expl_dat_df$Education),
                            Children = levels(expl_dat_df$Children))
colnames(final_matrix)
expl_dat_mat<-cbind(expl_dat_df, final_matrix[,-c(1:11, 58)])

# Please note, due to the randomization 
# processes combinations maybe different

all_sub_matrices<-sapply(1:nrow(combinations), 
                         function(x){
                           
                           print(x)
                           filt_dat<-expl_dat_mat[expl_dat_mat$Live_in==combinations[x, "Live_in"] &
                                                    expl_dat_mat$Gender==combinations[x, "Gender"] &
                                                    expl_dat_mat$Age==combinations[x, "Age"] &
                                                    expl_dat_mat$Education==combinations[x, "Education"] &
                                                    expl_dat_mat$Children==combinations[x, "Children"], 
                                                  6: ncol(expl_dat_mat)]
                           
                           if(nrow(filt_dat)>=20){
                             
                             
                             filt_net<-network_gen(filt_dat)
                             if(!is.null(filt_net))
                             {
                               links<-as_long_data_frame(filt_net)
                               
                               link_presence_vector<-sapply(1:ncol(potential_links), 
                                                            function(k){
                                                              # k = 10
                                                              # k = 359
                                                              
                                                              from <- potential_links[1, k]
                                                              to <- potential_links[2, k]
                                                              
                                                              # from = "Gtype_House"
                                                              # to = "Garea_Orna"
                                                              
                                                              out<-0
                                                              if(from %in% links$from_name & 
                                                                 to %in% links$to_name) 
                                                              {out =  as.numeric(links[links$from_name == from & links$to_name == to, "sign"])
                                                              out<-ifelse(length(out)==0, 0, out)
                                                              }
                                                              if(to %in% links$from_name & 
                                                                 from %in% links$to_name )
                                                              {out =  as.numeric(links[links$from_name == to & links$to_name == from, "sign"])
                                                              out<-ifelse(length(out)==0, 0, out)
                                                              }
                                                              out[1]}
                               )
                               link_presence_vector
                               
                             } else
                             {cat("Not enough data")}
                             
                             
                           }})

filtered_sub_matrices<-all_sub_matrices[!sapply(all_sub_matrices, is.null)]

all_networks_mat<-do.call("rbind", filtered_sub_matrices)
all_networks_mat<-as.data.frame(all_networks_mat)

colnames(all_networks_mat)<-potential_link_names

filt<-apply(all_networks_mat, 2, function(m) all(m==0))
networks_mat<-all_networks_mat[, !filt]

env_mat<-combinations[!sapply(all_sub_matrices, is.null), ]
RDA_mod<-rda(networks_mat~Live_in+Gender+Age+Education+Children, env_mat)

summary(RDA_mod)

RDA_mod$CCA$tot.chi
global_r2<-RsquareAdj(RDA_mod)    # Total variance explained by the RDA

anova.cca(RDA_mod, permutations = 999)
anova.cca(RDA_mod, by = "axis")   # Test which axis are significant
anova.cca(RDA_mod, by = "terms")  # Test which terms are significant

limits1_low<-min(RDA_mod$CCA$v[,"RDA1"])*0.5
limits1_high<-max(RDA_mod$CCA$v[,"RDA1"])*0.5
limits2_low<-min(RDA_mod$CCA$v[,"RDA2"])*0.5
limits2_high<-max(RDA_mod$CCA$v[,"RDA2"])*0.5

filt_points<-(RDA_mod$CCA$v[,"RDA1"]<limits1_low | 
                RDA_mod$CCA$v[,"RDA1"]>limits1_high) & 
  (RDA_mod$CCA$v[,"RDA2"]<limits2_low | 
     RDA_mod$CCA$v[,"RDA2"]>limits2_high) 


RDA_mod$CCA$v<-RDA_mod$CCA$v[filt_points,]
RDA_mod$CA$v<-RDA_mod$CA$v[filt_points,]

########################################
########### Visualization ##############
########################################

windows()
ggord(RDA_mod, env_mat$Education,
      polylntyp = "dashed", poly = F, 
      alpha_el=0.2, veclsz = 1, veccol = "red", 
      arrow=0.3, size=2,
      addsize = 3,
      labcol = "red",
      grp_title = "Education", coord_fix = F,
      ptslab = T, repel = T)
