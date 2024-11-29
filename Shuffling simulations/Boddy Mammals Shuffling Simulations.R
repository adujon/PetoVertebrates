rm(list=ls()) #remove all objects into memmory
library(readxl)
library(ape)
library(nlme)
library(geiger)
library(ggeffects)
library(caper)
library(phylolm)
library(phytools)
library(MuMIn)
library(ggplot2)
library(rr2)
library(DHARMa)
library(dplyr)
library(tidyr)
library(ggpubr)
library(grid)
library(brms)
library(parallel)
library(tidybayes)
library(bayesplot)

################################################################################

#setup import and export folders
export_dir <- "C:/Users/antoi/My Drive/Postdoc/Zoo data/Peto/models"
import_dir <- "C:/Users/antoi/My Drive/Postdoc/Zoo data/Peto/data"

##########################
### import phylogeny data ###
##########################

#import the phylogenetic tree for mammals
setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Peto/data")
phylo_all <- read.nexus("MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre") 

#rewrite the tip labels to match the species names in the data file
current_tip_labels <- phylo_all$tip.label
new_tip_labels <- c()

for(i in 1:length(current_tip_labels)){
  
  vec <- strsplit(current_tip_labels[i], "[_]")[[1]][1:3]
  
  if(vec[1] == ""){
    new_tip_labels[i] <- paste(vec[2], vec[3], sep = "_")
  } else {
    new_tip_labels[i] <- paste(vec[1], vec[2], sep = "_")
  }
  
  
}

df_label <- data.frame(current_tip_labels, new_tip_labels)
phylo_all$tip.label<-df_label[[2]][match(phylo_all$tip.label, df_label[[1]])]

##########################
### import cancer data ###
##########################

data <-  read.csv("Boddy cancer mammals.csv")
data <- data.frame(data) #convert to dataframe
colnames(data)[2] <- "Species" #rename species column
colnames(data)[6] <- "NNecropsies" #rename species column
colnames(data)[7] <- "NNeoplasia" #rename neoplasia column
colnames(data)[8] <- "NMalignant" #rename malignant column

row.names(data) <- data$Species #change row names to species names
data$logMass <- log10(data$adult_mass_kg) #log transform the body mass
data$logMaxLongevity <- log10(data$max_lifespan_yr * 12) #log transform the life expectancy

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#################### prepare the data for brms models ########################################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## fix the name for the donkeys Equus africanus asinus to match names in the csv
sel_donkey <- which(phylo_all$tip.label == "Equus_africanus")
phylo_all$tip.label[sel_donkey] <- "Equus_asinus"

## fix the name for the  Pygmy marmoset to match names in the csv
sel_donkey <- which(phylo_all$tip.label == "Cebuella_pygmaea")
phylo_all$tip.label[sel_donkey] <- "Callithrix_pygmaea"

#check the number of species left
phylo_brms <- phylo_all

#compute the var-covar matrix from the tree to use in brms
vcv_brms_OU <- vcv(corMartins(1, phylo_brms, fixed=FALSE,form=~phylo_brms$tip.label)) #OU model

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
############ Modelling the effect of logMass on neoplasia and longevity  #####################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

#define the number of simulations
n_simul <- 50

#######################
# neoplasia - body mass
#######################

#create an empty dataframe to store the results of the simulations
simul_df <- data.frame(Data = NULL, TumorType = NULL, RiskFactor = NULL, Simulation = NULL, Estimate = NULL, L95CI = NULL, U95CI = NULL)


for(i in 1:n_simul){
  
  #print current simulation
  print(paste("Running simulation", i))
  
  #copy the original data
  data_simul <- data
  
  #compute prevalences to keep them constant accross species
  prevalences <- data_simul$NNeoplasia/data_simul$NNecropsies
  
  #shuffle the number of necropsies accross species
  data_simul$NNecropsies <- sample(data_simul$NNecropsies, replace = FALSE)
  
  #calculate the number of animal with tumours based on the new number of necropsies
  data_simul$NNeoplasia <- round(data_simul$NNecropsies *  prevalences)
  
  #fit the phylogenetic regression model using brms
  fit_logmass_neoplasia <- brm(NNeoplasia | trials(NNecropsies) ~ logMass + (1|gr(Species, cov = vcv_brms_OU)),
                               family =binomial(),
                               data = data_simul,
                               data2 = list(vcv_brms_OU = vcv_brms_OU), 
                               chains = 4,
                               cores = 24,
                               iter = 10000, warmup = 5000,
                               thin = 10,
                               seed = 1,
                               control=list(max_treedepth = 12, adapt_delta = 0.99))
  
  #get the effect sizes for the current model
  res <- summary(fit_logmass_neoplasia)$fixed
  
  #append the dataframe with the results
  res_df <- data.frame(Data = "Boddy et al.", TumorType = "Neoplasia", RiskFactor = "Mass", Simulation = i, Estimate = res[2,1], U95CI = res[2,3],  L95CI = res[2,4])
  simul_df <- rbind(simul_df, res_df)
  
}

setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Peto/scripts/Binomial models scripts/Shuffling simulations/Results")
write.csv(simul_df, "Boddy Neoplasia Mammal Mass.csv", row.names = F)

#######################
# malignancy - body mass
#######################

#create an empty dataframe to store the results of the simulations
simul_df <- data.frame(Data = NULL, TumorType = NULL, RiskFactor = NULL, Simulation = NULL, Estimate = NULL, L95CI = NULL, U95CI = NULL)

for(i in 1:n_simul){
  
  #print current simulation
  print(paste("Running simulation for proportion", i))
  
  #copy the original data
  data_simul <- data
  
  #compute prevalences to keep them constant accross species
  prevalences <- data_simul$NMalignant/data_simul$NNecropsies
  
  #shuffle the number of necropsies accross species
  data_simul$NNecropsies <- sample(data_simul$NNecropsies, replace = FALSE)
  
  #calculat the number of animal with tumours based on the new number of necropsies
  data_simul$NMalignant <- round(data_simul$NNecropsies *  prevalences)
  
  #fit the phylogenetic regression model using brms
  fit_logmass_malignant<- brm(NMalignant | trials(NNecropsies) ~ logMass + (1|gr(Species, cov = vcv_brms_OU)),
                               family =binomial(),
                               data = data_simul,
                               data2 = list(vcv_brms_OU = vcv_brms_OU), 
                               chains = 4,
                               cores = 24,
                               iter = 10000, warmup = 5000,
                               thin = 10,
                               seed = 1,
                               control=list(max_treedepth = 12, adapt_delta = 0.99))
  
  #get the effect sizes for the current model
  res <- summary(fit_logmass_malignant)$fixed
  
  #append the dataframe with the results
  res_df <- data.frame(Data = "Boddy et al.", TumorType = "Malignant", RiskFactor = "Mass", Simulation = i, Estimate = res[2,1], U95CI = res[2,3],  L95CI = res[2,4])
  simul_df <- rbind(simul_df, res_df)
  
}

setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Peto/scripts/Binomial models scripts/Shuffling simulations/Results")
write.csv(simul_df, "Boddy Malignant Mammal Mass.csv", row.names = F)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
############ Modelling the effect of logmaxLogevity on neoplasia and malignancy  #############
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##


#######################
# neoplasia - longevity
#######################

#create an empty dataframe to store the results of the simulations
simul_df <- data.frame(Data = NULL, TumorType = NULL, RiskFactor = NULL, Simulation = NULL, Estimate = NULL, L95CI = NULL, U95CI = NULL)


for(i in 1:n_simul){
  
  #print current simulation
  print(paste("Running simulation", i))
  
  #copy the original data
  data_simul <- data
  
  #compute prevalences to keep them constant accross species
  prevalences <- data_simul$NNeoplasia/data_simul$NNecropsies
  
  #shuffle the number of necropsies accross species
  data_simul$NNecropsies <- sample(data_simul$NNecropsies, replace = FALSE)
  
  #calculat the number of animal with tumours based on the new number of necropsies
  data_simul$NNeoplasia <- round(data_simul$NNecropsies *  prevalences)
  
  #fit the phylogenetic regression model using brms
  fit_logMaxLongevity_neoplasia <- brm(NNeoplasia | trials(NNecropsies) ~ logMaxLongevity + (1|gr(Species, cov = vcv_brms_OU)),
                               family =binomial(),
                               data = data_simul,
                               data2 = list(vcv_brms_OU = vcv_brms_OU), 
                               chains = 4,
                               cores = 24,
                               iter = 10000, warmup = 5000,
                               thin = 10,
                               seed = 1,
                               control=list(max_treedepth = 12, adapt_delta = 0.99))
  
  #get the effect sizes for the current model
  res <- summary(fit_logMaxLongevity_neoplasia)$fixed
  
  #append the dataframe with the results
  res_df <- data.frame(Data = "Boddy et al.", TumorType = "Neoplasia", RiskFactor = "Longevity", Simulation = i, Estimate = res[2,1], U95CI = res[2,3],  L95CI = res[2,4])
  simul_df <- rbind(simul_df, res_df)
  
}

setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Peto/scripts/Binomial models scripts/Shuffling simulations/Results")
write.csv(simul_df, "Boddy Neoplasm Mammal Longevity.csv", row.names = F)

#######################
# malignancy - longevity
#######################


#create an empty dataframe to store the results of the simulations
simul_df <- data.frame(Data = NULL, TumorType = NULL, RiskFactor = NULL, Simulation = NULL, Estimate = NULL, L95CI = NULL, U95CI = NULL)


for(i in 1:n_simul){
  
  #print current simulation
  print(paste("Running simulation for proportion", i))
  
  #copy the original data
  data_simul <- data
  
  #compute prevalences to keep them constant accross species
  prevalences <- data_simul$NMalignant/data_simul$NNecropsies
  
  #shuffle the number of necropsies accross species
  data_simul$NNecropsies <- sample(data_simul$NNecropsies, replace = FALSE)
  
  #calculat the number of animal with tumours based on the new number of necropsies
  data_simul$NMalignant <- round(data_simul$NNecropsies *  prevalences)
  
  #fit the phylogenetic regression model using brms
  fit_logMaxLongevity_malignant <- brm(NMalignant | trials(NNecropsies) ~ logMaxLongevity + (1|gr(Species, cov = vcv_brms_OU)),
                                       family =binomial(),
                                       data = data_simul,
                                       data2 = list(vcv_brms_OU = vcv_brms_OU), 
                                       chains = 4,
                                       cores = 24,
                                       iter = 10000, warmup = 5000,
                                       thin = 10,
                                       seed = 1,
                                       control=list(max_treedepth = 12, adapt_delta = 0.99))
  
  #get the effect sizes for the current model
  res <- summary(fit_logMaxLongevity_malignant)$fixed
  
  #append the dataframe with the results
  res_df <- data.frame(Data = "Boddy et al.", TumorType = "Malignant", RiskFactor = "Longevity", Simulation = i, Estimate = res[2,1], U95CI = res[2,3],  L95CI = res[2,4])
  simul_df <- rbind(simul_df, res_df)
  
}

setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Peto/scripts/Binomial models scripts/Shuffling simulations/Results")
write.csv(simul_df, "Boddy Malignant Mammal Longevity.csv", row.names = F)
