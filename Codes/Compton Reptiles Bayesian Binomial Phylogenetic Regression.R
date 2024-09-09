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

#import the phylogenetic tree 
setwd(import_dir)
phylo_all <- read.tree("Compton_bird_phylogeny.tre") 

##########################
### import phylogeny data ###
##########################

#import the phylogenetic tree for mammals
setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Peto/data")
phylo_all <- read.tree("squam_shl_new.tre") 

##########################
### import cancer data ###
##########################

data <-  read.csv("Compton et al Vertebrates.csv")
data <- data.frame(data) #convert to dataframe
colnames(data)[11] <- "NNecropsies"  
colnames(data)[12] <- "NNeoplasia"  
colnames(data)[16] <- "NMalignant"  

#subset the datset to only keep mammals
data <- subset(data, data$Class == "Reptilia")

#change the dataframe rown names to species names
row.names(data) <- data$Species

#replace the -1 in the mass and maxlongevity columns by NA
data$adult_weight.g[data$adult_weight.g == -1] <- NA
data$max_longevity.months.[data$max_longevity.months. == -1] <- NA

#convert the body mass to kg
data$adult_mass_kg <- data$adult_weight.g / 1000

data$logMass <- log(data$adult_mass_kg) #log transform the body mass
data$logNecropsies <- log(data$NNecropsies) #log transform the number of necropsies
data$logMaxLongevity <- log(data$max_longevity.months.) #log transform the  Max Longevity in months

#remove species with no mass or longevity data
sel_rm <- c()
for(i in 1:nrow(data)){
  
  if(is.na(data$adult_mass_kg[i])&is.na(data$logMaxLongevity[i])){
    
    sel_rm <- c( sel_rm, i)
  }
  
}

data <- data[-sel_rm,]

#remove species not in tree
data <- data[-which(data$Species == "Chelus_fimbriatus"),]
data <-data[-which(data$Species == "Gopherus_agassizii"),]
data <-data[-which(data$Species == "Paleosuchus_palpebrosus"),]
data <-data[-which(data$Species == "Testudo_kleinmanni"),]

##############################################################################################
############################# visualize the data ########################################
##############################################################################################

#histogram of the body mass
hist(data$adult_mass_kg)
hist(data$logMass) #better when log transformed

#histogram of the number of necropsies
hist(data$NNecropsies)

#histovrame log longevity
hist(data$logMaxLongevity)

#compute the number of species with no individuals that died from malignant cancer
length(which(data$NNeoplasia == 0)) / nrow(data) 

#histogram of the number of individuals that died from malignant cancer
hist(data$NNeoplasia)

##############################################################################################
#################### prepare the data for brms models ########################################
##############################################################################################


#check which species are both in the main dataframe and the phylogeny tree
obj_brms<- name.check(phylo_all,data)
obj_brms

#drop the species from the tree for which we have no cancer data
phylo_brms <- drop.tip(phylo_all, obj_brms$tree_not_data)

#compute the var-covar matrix from the tree to use in brms
vcv_brms <- vcv(phylo_brms, model = "Brownian")


##############################################################################################
#######################  Modelling the effect of logMass on neoplasia ##########################
##############################################################################################

#fit the phylogenetic regression model using brms
fit_logmass_neoplasia <- brm(NNeoplasia | trials(NNecropsies) ~ logMass + (1|gr(Species, cov = vcv_brms)),
                             family =binomial(),
                             data = data,
                             data2 = list(vcv_brms = vcv_brms), 
                             chains = 4,
                             cores = 32,
                             iter = 10000, warmup = 5000,
                             thin = 10,
                             seed = 1,
                             control=list(max_treedepth = 12, adapt_delta = 0.99))

#inspect the chains for improper mixing
plot(fit_logmass_neoplasia)

#check if the iteration in the chains are not autocorrelated
#this will plot the autocorrelation function for each of the four chains
mcmc_acf(fit_logmass_neoplasia , pars = "b_Intercept", lags = 10)
mcmc_acf(fit_logmass_neoplasia , pars = "b_logMass", lags = 10)
mcmc_acf(fit_logmass_neoplasia , pars = "sd_Species__Intercept", lags = 10)

#check if the posterior distribution matches the observed distribution
pp_check(fit_logmass_neoplasia , ndraws = 100)

#plot the conditional effect of mass for this model
cond_effects <- conditional_effects(fit_logmass_neoplasia , "logMass")
plot(cond_effects)

#export the model
setwd(export_dir)
saveRDS(fit_logmass_neoplasia , "Compton_reptiles_neoplasm_logMass.rds")

##############################################################################################
#######################  Modelling the effect of logMass on malignancy #######################
##############################################################################################

#fit the phylogenetic regression model using brms
fit_logmass_malignant <- brm(NMalignant | trials(NNecropsies) ~ logMass + (1|gr(Species, cov = vcv_brms)),
                            family =binomial(),
                            data = data,
                            data2 = list(vcv_brms = vcv_brms), 
                            chains = 4,
                            cores = 32,
                            iter = 10000, warmup = 5000,
                            thin = 10,
                            seed = 1,
                            control=list(max_treedepth = 12, adapt_delta = 0.99))

#inspect the chains for improper mixing
plot(fit_logmass_malignant)

#check if the iteration in the chains are not autocorrelated
#this will plot the autocorrelation function for each of the four chains
mcmc_acf(fit_logmass_malignant, pars = "b_Intercept", lags = 10)
mcmc_acf(fit_logmass_malignant, pars = "b_logMass", lags = 10)
mcmc_acf(fit_logmass_malignant, pars = "sd_Species__Intercept", lags = 10)

#check if the posterior distribution matches the observed distribution
pp_check(fit_logmass_malignant, ndraws = 100)

#plot the conditional effect of mass for this model
cond_effects <- conditional_effects(fit_logmass_malignant, "logMass")
plot(cond_effects)

#export the model
setwd(export_dir)
saveRDS(fit_logmass_malignant, "Compton_reptiles_malignant_logMass.rds")


##############################################################################################
#######################  Modelling the effect of logMaxLongevity on neoplasia ##############
##############################################################################################

#fit the phylogenetic regression model using brms
fit_logMaxLongevity_neoplasia <- brm(NNeoplasia | trials(NNecropsies) ~ logMaxLongevity + (1|gr(Species, cov = vcv_brms)),
                             family =binomial(),
                             data = data,
                             data2 = list(vcv_brms = vcv_brms), 
                             chains = 4,
                             cores = 32,
                             iter = 10000, warmup = 5000,
                             thin = 10,
                             seed = 1,
                             control=list(max_treedepth = 12, adapt_delta = 0.99))

#inspect the chains for improper mixing
plot(fit_logMaxLongevity_neoplasia)

#check if the iteration in the chains are not autocorrelated
#this will plot the autocorrelation function for each of the four chains
mcmc_acf(fit_logMaxLongevity_neoplasia , pars = "b_Intercept", lags = 10)
mcmc_acf(fit_logMaxLongevity_neoplasia , pars = "b_logMaxLongevity", lags = 10)
mcmc_acf(fit_logMaxLongevity_neoplasia , pars = "sd_Species__Intercept", lags = 10)

#check if the posterior distribution matches the observed distribution
pp_check(fit_logMaxLongevity_neoplasia , ndraws = 100)

#plot the conditional effect of mass for this model
cond_effects <- conditional_effects(fit_logMaxLongevity_neoplasia , "logMaxLongevity")
plot(cond_effects)

#export the model
setwd(export_dir)
saveRDS(fit_logMaxLongevity_neoplasia , "Compton_reptiles_neoplasm_logMaxLongevity.rds")

##############################################################################################
#######################  Modelling the effect of logMaxLongevity on malignancy #######################
##############################################################################################

#fit the phylogenetic regression model using brms
fit_logMaxLongevity_malignant <- brm(NMalignant | trials(NNecropsies) ~ logMaxLongevity + (1|gr(Species, cov = vcv_brms)),
                             family =binomial(),
                             data = data,
                             data2 = list(vcv_brms = vcv_brms), 
                             chains = 4,
                             cores = 32,
                             iter = 10000, warmup = 5000,
                             thin = 10,
                             seed = 1,
                             control=list(max_treedepth = 12, adapt_delta = 0.99))

#inspect the chains for improper mixing
plot(fit_logMaxLongevity_malignant)

#check if the iteration in the chains are not autocorrelated
#this will plot the autocorrelation function for each of the four chains
mcmc_acf(fit_logMaxLongevity_malignant, pars = "b_Intercept", lags = 10)
mcmc_acf(fit_logMaxLongevity_malignant, pars = "b_logMaxLongevity", lags = 10)
mcmc_acf(fit_logMaxLongevity_malignant, pars = "sd_Species__Intercept", lags = 10)

#check if the posterior distribution matches the observed distribution
pp_check(fit_logMaxLongevity_malignant, ndraws = 100)

#plot the conditional effect of mass for this model
cond_effects <- conditional_effects(fit_logMaxLongevity_malignant, "logMaxLongevity")
plot(cond_effects)

#export the model
setwd(export_dir)
saveRDS(fit_logMaxLongevity_malignant, "Compton_reptiles_malignant_logMaxLongevity.rds")

