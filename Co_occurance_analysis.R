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

final_matrix <- read.csv("final_matrix.csv", header=T, encoding='UTF-8')

########################################
########### CO-OCCURANCE ##############
########################################

### Network without socio-demoraphic variables

colnames(final_matrix)
final_matrix_wo_sd_vars<-final_matrix[, -c(1:11, 58)] #delet socio-d.

net_wo_sd_vars<-network_gen(final_matrix_wo_sd_vars)

E(net_wo_sd_vars)$color[E(net_wo_sd_vars)$color != "blue"] <- "darksalmon"
E(net_wo_sd_vars)$color[E(net_wo_sd_vars)$color == "blue"] <- "cornflowerblue"

windows()
par(mfrow=c(1,2))
#pdf("network_with_soc.pdf")
plot(net_wo_sd_vars$mod, net_wo_sd_vars, #mod!
     layout = layout.circle(net_wo_sd_vars, 
                            order = order(V(net_wo_sd_vars)$mod)),
     vertex.label.color="black",
     edge.color = E(net_wo_sd_vars)$color,
     vertex.label.dist = 1, vertex.label.cex = .7,
     vertex.label.degree = -pi/2)
#dev.off()

#pdf("fig_network_with_soc_pos.pdf")
plot(net_wo_sd_vars$mod_poz, net_wo_sd_vars, #mod_poz!
     layout = layout.circle(net_wo_sd_vars, 
                            order = order(V(net_wo_sd_vars)$mod_poz)),
     vertex.label.color="black",
     mark.col = adjustcolor(c("chocolate1", "deepskyblue", "chartreuse3", 
                              "darkgoldenrod2", "cornflowerblue"), alpha.f = 0.3),
     mark.border = adjustcolor(c("chocolate1", "deepskyblue", "chartreuse3", 
                                 "darkgoldenrod2", "cornflowerblue"), alpha.f = 1),
     edge.color = E(net_wo_sd_vars)$color,
     vertex.label.dist = 1, vertex.label.cex = .7,
     vertex.label.degree = -pi/2)
#dev.off()

# network - only negative edges
neg_net<-delete_edges(net_wo_sd_vars, E(net_wo_sd_vars)[sign==1])

# network - only positive edges
poz_net<-delete_edges(net_wo_sd_vars, E(net_wo_sd_vars)[sign== -1])

# number of connections
sort(degree(neg_net))
sort(degree(net_wo_sd_vars))
sort(degree(poz_net))

# betweenness
sort(betweenness(neg_net))
sort(betweenness(net_wo_sd_vars))
sort(betweenness(poz_net))

