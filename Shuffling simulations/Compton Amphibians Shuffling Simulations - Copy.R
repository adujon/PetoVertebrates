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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import phylogeny data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#import the phylogenetic tree for squamates
setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Peto/data")
phylo_all <- read.tree("amph_shl_new_Consensus_7238.tre") 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import cancer data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data <-  read.csv("Compton et al amphibians.csv")
data <- data.frame(data) #convert to dataframe
colnames(data)[11] <- "NNecropsies"  
colnames(data)[12] <- "NNeoplasia"  
colnames(data)[16] <- "NMalignant"  

#change the dataframe rown names to species names
row.names(data) <- data$Species

#convert the body mass to kg and maximum longevity to months
data$adult_mass_kg <- data$Body_Mass_Grams / 1000
data$max_longevity_months <- data$Maximum_Longevity_Years * 12

data$logMass <- log10(data$adult_mass_kg) #log transform the body mass
data$logMaxLongevity <- log10(data$max_longevity_months) #log transform the  Max Longevity in months

#remove species with no mass or longevity data
sel_rm <- c()
for(i in 1:nrow(data)){
  
  if(is.na(data$adult_mass_kg[i])&is.na(data$logMaxLongevity[i])){
    
    sel_rm <- c( sel_rm, i)
  }
  
}

data <- data[-sel_rm,]

#remove species not in tree
data <- data[-which(data$Species == "Lithobates_pipiens"),]
data <- data[-which(data$Species == "Lithobates_catesbeianus"),]

##############################################################################################
#################### prepare the data for brms models ########################################
##############################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Prepare the data for brms models ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#check which species are both in the main dataframe and the phylogeny tree
obj_brms<- name.check(phylo_all,data)
obj_brms


## fix the name for the Mississippi gopher frog
sel_ran <- which(phylo_all$tip.label == "Rana_sevosa")
phylo_all$tip.label[sel_ran] <- "Lithobates_sevosus"


#recheck the conflicts have been revolved
obj_brms<- name.check(phylo_all,data)
obj_brms

#drop the species from the tree for which we have no cancer data
phylo_brms <- drop.tip(phylo_all, obj_brms$tree_not_data)


#compute the var-covar matrix from the tree to use in brms
vcv_brms_OU <- vcv(corMartins(1, phylo_brms, fixed=FALSE,form=~phylo_brms$tip.label)) #OU model
vcv_brms_Pagel <- vcv(corPagel(value = 1, phy = phylo_brms, fixed = FALSE, form=~phylo_brms$tip.label)) #Pagel model

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#######################  Modelling the effect of logMass on neoplasia and malignancy  ########
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
  
  #calculat the number of animal with tumours based on the new number of necropsies
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
  res_df <- data.frame(Data = "Compton et al.", TumorType = "Neoplasia", RiskFactor = "Mass", Simulation = i, Estimate = res[2,1], U95CI = res[2,3],  L95CI = res[2,4])
  simul_df <- rbind(simul_df, res_df)
  
}

setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Peto/scripts/Binomial models scripts/Shuffling simulations/Results")
write.csv(simul_df, "Compton Neoplasia Amphibians Mass.csv", row.names = F)

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
  res_df <- data.frame(Data = "Compton et al.", TumorType = "Malignant", RiskFactor = "Mass", Simulation = i, Estimate = res[2,1], U95CI = res[2,3],  L95CI = res[2,4])
  simul_df <- rbind(simul_df, res_df)
  
}

setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Peto/scripts/Binomial models scripts/Shuffling simulations/Results")
write.csv(simul_df, "Compton Malignant Amphibians Mass.csv", row.names = F)

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
  res_df <- data.frame(Data = "Compton et al.", TumorType = "Neoplasia", RiskFactor = "Longevity", Simulation = i, Estimate = res[2,1], U95CI = res[2,3],  L95CI = res[2,4])
  simul_df <- rbind(simul_df, res_df)
  
}

setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Peto/scripts/Binomial models scripts/Shuffling simulations/Results")
write.csv(simul_df, "Compton Neoplasm Amphibians Longevity.csv", row.names = F)

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
  res_df <- data.frame(Data = "Compton et al.", TumorType = "Malignant", RiskFactor = "Longevity", Simulation = i, Estimate = res[2,1], U95CI = res[2,3],  L95CI = res[2,4])
  simul_df <- rbind(simul_df, res_df)
  
}

setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Peto/scripts/Binomial models scripts/Shuffling simulations/Results")
write.csv(simul_df, "Compton Malignant Amphibians Longevity.csv", row.names = F)
