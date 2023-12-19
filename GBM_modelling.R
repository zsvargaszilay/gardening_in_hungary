########################################
############## PACKAGES ################
########################################

library(dplyr)
library(scales)
library(rsample)      # data splitting 
library(gbm)          # basic implementation

########################################
############ READ DATA #################
########################################

mymatrix <- read.csv("matrix_for_gbm.csv", header=T, encoding='UTF-8')

########################################
############# FIGURE 3 #################
########################################

### Excluding highly subjective questions ### 

mymatrix <- mymatrix[, -c(44:45)]

colnames(mymatrix) <- c("Get_shop", "Get_pers", "Get_colle",
                        "Ferti_Synt", "Pesti_No", "Pesti_Synt",
                        "PollSup_Nhab", "PollSup_Food", "PollSup_Wate",
                        "PollSup_Ahab", "Info_Web", "Info_Pers", "Info_Trad",
                        "Live_in", "Gender", "Age", "Education", "Child", "Garden_size",
                        "Pond", "G_time", "G_type", "Garea_Vege", "Garea_Herb",
                        "Garea_Frui", "Garea_Grap", "Garea_Orna", "Garea_egre",
                        "Garea_Lawn", "Garea_Undi", "Habi_Beau", "Habi_Cons", "Habi_Money", "Habi_Prod",
                        "Motivation", "Glife_Pl", "Glife_Bi", "Glife_In", "Mowing_Freq", "Unmowp", "Oplants", "Herbicid",
                        "Growing", "SeePoll_Med", "SeePoll_Lot", "SeePoll_Few")

### biodiversity-positive activites ###

mymatrix_for_positives <- mymatrix[ , -c(4, 6, 30, 39, 42)] # delete biodiv-positive activities

# Merging: Food support + Water support = PollSup_Sources
mymatrix_for_positives$PollSup_Sources <- apply(mymatrix_for_positives[, 6:7], 1, function(x){
  ifelse(sum(x) > 0, 1, 0)
})
mymatrix_for_positives <- mymatrix_for_positives[ , -c(6:7)]

# create a vector
mymatrix_for_positives_sum <- apply(mymatrix_for_positives[, c(4:6, 16, 34, 40)], 1, 
                                    function(x){
                                      sum(as.numeric(x))
                                    })

### biodiversity-negative activites ###

mymatrix_for_negatives <- mymatrix[ , -c(5, 7:10, 20, 40)] # delete biodiv-negative activities

### create a vector ###
mymatrix_for_negatives_sum <- apply(mymatrix_for_negatives[, c(4:5, 24, 33, 35)], 1, 
                                    function(x){
                                      sum(as.numeric(x))
                                    })

########################################

### scaling ###
mymatrix_for_positives_sum <- rescale(mymatrix_for_positives_sum, to = c(0,1))
mymatrix_for_negatives_sum <- rescale(mymatrix_for_negatives_sum, to = c(0,1))

mymatrix_for_negatives_sum <- mymatrix_for_negatives_sum*-1

pos_and_neg_vector <- mymatrix_for_positives_sum + mymatrix_for_negatives_sum

### delet the coloumns of biodiversity-positive and - negative activities ###

deleted_cols <- c(c(4, 6, 30, 39, 42), c(5, 7:10, 20, 40))
colnames(mymatrix[ , deleted_cols])
#> colnames(mymatrix[ , deleted_cols])
#[1] "Ferti_Synt"   "Pesti_Synt"   "Garea_Undi"   "Mowing_Freq"  "Herbicid"    
#[6] "Pesti_No"     "PollSup_Nhab" "PollSup_Food" "PollSup_Wate" "PollSup_Ahab"
#[11] "Pond"         "Unmowp"

new_mymatrix <- mymatrix[ , -deleted_cols]

### Add the biodiversity friendliness score (BDF score) to the matrix

new_mymatrix$sum <- pos_and_neg_vector 

### Transform all of variables to factor (except for cols of BDF score) ###

new_mymatrix[1:34] <- lapply(new_mymatrix[1:34], factor)

### Filtering: variables in which response agreement was over 95% were removed ### 

mymatrix_for_filt <- apply(new_mymatrix[, 1:34], 2, function(x){
  largest_level <- rev(sort(table(x)))[1] # reverse (rev) the sequence (sort)
  !largest_level > length(x)*0.95 
})

mymatrix_for_filt <- c(mymatrix_for_filt, TRUE)

new_mymatrix_final <- new_mymatrix[,mymatrix_for_filt]

sapply(new_mymatrix_final, levels)

########################################
########### GBM MODELLING ##############
########################################

# Create training (70%) and test (30%) sets

set.seed(123) # for reproducibility
ames_split <- initial_split(new_mymatrix_final, prop = .7)
ames_train <- training(ames_split)
ames_test  <- testing(ames_split)

set.seed(123)

# train GBM model
gbm.fit <- gbm(
  formula = sum ~ .,
  distribution = "gaussian",
  data = ames_train,
  n.trees = 10000,
  interaction.depth = 1,
  shrinkage = 0.001,
  cv.folds = 5,
  n.cores = NULL,
  verbose = FALSE
)  

print(gbm.fit)
#The best cross-validation iteration was 9997.
#There were 31 predictors of which 31 had non-zero influence.

# get MSE and compute RMSE
sqrt(min(gbm.fit$cv.error))
# [1] 0.3320361

# plot loss function as a result of n trees added to the ensemble
windows()
gbm.perf(gbm.fit, method = "cv")
# 9997 tree

# create hyperparameter grid
hyper_grid <- expand.grid(
  shrinkage = c(.01, .1, .3),
  interaction.depth = c(1, 3, 5),
  n.minobsinnode = c(5, 10, 15),
  bag.fraction = c(.65, .8, 1), 
  optimal_trees = 9997,      
  min_RMSE = 0.3320361
)

# total number of combinations
nrow(hyper_grid)
# [1] 81

# randomize data --> emiatt eltérő lehet az eredmény
random_index <- sample(1:nrow(ames_train), nrow(ames_train))
random_ames_train <- ames_train[random_index, ]

# grid search 
for(i in 1:nrow(hyper_grid)) {
  
  # reproducibility
  set.seed(123)
  
  # train model
  gbm.tune <- gbm(
    formula = sum ~ .,
    distribution = "gaussian",
    data = random_ames_train,
    n.trees = 9997,
    interaction.depth = hyper_grid$interaction.depth[i],
    shrinkage = hyper_grid$shrinkage[i],
    n.minobsinnode = hyper_grid$n.minobsinnode[i],
    bag.fraction = hyper_grid$bag.fraction[i],
    train.fraction = .75,
    n.cores = NULL, # will use all cores by default
    verbose = FALSE
  )
  
  # add min training error and trees to grid
  hyper_grid$optimal_trees[i] <- which.min(gbm.tune$valid.error)
  hyper_grid$min_RMSE[i] <- sqrt(min(gbm.tune$valid.error))
}

hyper_grid %>% 
  dplyr::arrange(min_RMSE) %>%
  head(10)

###
# modify hyperparameter grid
hyper_grid <- expand.grid(
  shrinkage = c(.01, .05, .1),
  interaction.depth = c(3, 5, 7),
  n.minobsinnode = c(5, 7, 10),
  bag.fraction = c(.65, .8, 1), 
  optimal_trees = 85,               # a place to dump results
  min_RMSE = 0.3352214              # a place to dump results
)

# total number of combinations
nrow(hyper_grid)
# shrinkage interaction.depth n.minobsinnode bag.fraction optimal_trees  min_RMSE
#1        0.1                 3              5         0.80            85 0.3352214

set.seed(123)

# train GBM model
gbm.fit.final <- gbm(
  formula = sum ~ .,
  distribution = "gaussian",
  data = ames_train,
  n.trees = 85,
  interaction.depth = 3,
  shrinkage = 0.1,
  n.minobsinnode = 5,
  bag.fraction = .80, 
  train.fraction = 1,
  n.cores = NULL,
  verbose = FALSE
)  

########################################
####### LOAD our final model ###########
########################################
# Due to the randomization processes models maybe different #
# If you need the exact model we used, load this file #

load("gbm_fit_final.RDA")

########################################
###########   SHAP   ###################
########################################

library(shapviz)
library(kernelshap)
source("sv_importance_G.R") # edited function from https://cran.r-project.org/web/packages/shapviz/shapviz.pdf
library(ggpirate)

windows()
sapply(ames_test, class)

shap_values <- kernelshap(gbm.fit.final, 
                          X = new_mymatrix_final[,-ncol(new_mymatrix_final)], 
                          bg_X = ames_test)

sv <- shapviz(shap_values)
sv_importance(sv, kind = "no")

sapply(sv$X, levels)

windows()
sv2<-sv
colnames(sv)
sv<-sv2
new_colnames <- c("Getting seeds/plants in person*", "Getting seeds/plants via collection*", "Informaton from Internet*", 
                  "Information through personal link*", "Informaton through tradicional media*", "Type of residence", 
                  "Gender", "Age",  "Education level*", "Having children*", "Garden size", 
                  "Duration of gardening experience", "Garden type", "Having vegetables*",  "Having herbs*",  
                  "Having fruit trees*",  "Having grapes*", "Having ornamental*", "Having evergreens*",  "Having lawn*", 
                  "Motivation: beautiful garden*", "Motivation: nature conservation*", 
                  "Motivation: production*", "Garden perception*", "Knowledge about plants*",
                  "Knowledge about birds*", "Knowledge about insects*", "Planting ornamental plants*", 
                  "Growing crops", "Observing several pollinators groups*", "Observing many pollinator groups*")
colnames(sv$X)[1:31]<-new_colnames
colnames(sv$S)[1:31]<-new_colnames

sv_fig = sv_importance_G(sv, show_numbers = TRUE, kind = "both",
                         alpha = 0.4, fill = "lightgray", 
                         color_bar_title = "Factor levels",
                         max_display = 10)+ #  the first ten variables
                          #max_display = 31) #  all variables
  scale_color_brewer(palette = "Set2", 
                     labels = c("I", "II", "III", "IV", "V", "VI", "VII"))+
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5)))+
  theme_minimal()


ggsave("sv_fig.pdf", sv_fig + theme(text = element_text(size=20)),
       width = 18, height = 10, dpi = 300)
ggsave("sv_fig.png", sv_fig + theme(text = element_text(size=20)),
       width = 18, height = 10, dpi = 300)
ggsave("sv_fig_full_variables.pdf", sv_fig_full_variables + 
         theme(text = element_text(size=20)), width = 18, height = 15, dpi = 300)
ggsave("sv_fig_full_variables.png", sv_fig_full_variables +
         theme(text = element_text(size=20)), width = 18, height = 15, dpi = 300)




