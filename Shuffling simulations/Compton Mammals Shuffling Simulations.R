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

#import the phylogenetic tree for mammals
setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Peto/data")
phylo_all <- read.nexus("MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre") 

#rewrite the tip labels to match the species names in the datafile
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

data <-  read.csv("Compton et al Vertebrates.csv")
data <- data.frame(data) #convert to dataframe
colnames(data)[11] <- "NNecropsies"  
colnames(data)[12] <- "NNeoplasia"  
colnames(data)[16] <- "NMalignant"  

#subset the datset to only keep mammals
data <- subset(data, data$Class == "Mammalia")

#change the dataframe rown names to species names
row.names(data) <- data$Species

#replace the -1 in the mass and maxlongevity columns by NA
data$adult_weight.g[data$adult_weight.g == -1] <- NA
data$max_longevity.months.[data$max_longevity.months. == -1] <- NA

#convert the body mass to kg
data$adult_mass_kg <- data$adult_weight.g / 1000

data$logMass <- log10(data$adult_mass_kg) #log transform the body mass
data$logNecropsies <- log10(data$NNecropsies) #log transform the number of necropsies
data$logMaxLongevity <- log10(data$max_longevity.months.) #log transform the  Max Longevity in months

##############################################################################################
#################### prepare the data for brms models ########################################
##############################################################################################

## fix the name for the red-necked wallaby to match names in the csv
sel_wall <- which(phylo_all$tip.label == "Macropus_rufogriseus")
phylo_all$tip.label[sel_wall] <- "Notamacropus_rufogriseus"

## fix the name for the parma wallaby to match names in the csv
sel_wall <- which(phylo_all$tip.label == "Macropus_parma")
phylo_all$tip.label[sel_wall] <- "Notamacropus_parma"

## fix the name for the Common wallaroo to match names in the csv
sel_wall <- which(phylo_all$tip.label == "Macropus_robustus")
phylo_all$tip.label[sel_wall] <- "Osphranter_robustus"

## fix the name for the Red kangaroo to match names in the csv
sel_kan <- which(phylo_all$tip.label == "Macropus_rufus")
phylo_all$tip.label[sel_kan] <- "Osphranter_rufus"

## fix the name for the Giant eland to match names in the csv
sel_gia <- which(phylo_all$tip.label == "Tragelaphus_derbianus")
phylo_all$tip.label[sel_gia] <- "Taurotragus_derbianus"

#check which species are both in the main dataframe and the phylogeny tree
obj_brms<- name.check(phylo_all,data)
obj_brms

#drop the species from the tree for which we have no cancer data
phylo_brms <- drop.tip(phylo_all, obj_brms$tree_not_data)

#compute the var-covar matrix from the tree to use in brms
vcv_brms_OU <- vcv(corMartins(1, phylo_brms, fixed=FALSE,form=~phylo_brms$tip.label)) #OU model

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
write.csv(simul_df, "Compton Neoplasia Mammal Mass.csv", row.names = F)

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
write.csv(simul_df, "Compton Malignant Mammal Mass.csv", row.names = F)

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
write.csv(simul_df, "Compton Neoplasm Mammal Longevity.csv", row.names = F)

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
write.csv(simul_df, "Compton Malignant Mammal Longevity.csv", row.names = F)
