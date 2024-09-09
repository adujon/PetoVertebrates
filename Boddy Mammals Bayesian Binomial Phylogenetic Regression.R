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
data$logMass <- log(data$adult_mass_kg) #log transform the body mass
data$logMaxLongevity <- log(data$max_lifespan_yr * 12) #log transform the life expectancy

##############################################################################################
############################# visualize the data ########################################
##############################################################################################

#histogram of the body mass
hist(data$adult_mass_kg)
hist(data$logMass) #better when log transformed

#histogram of the number of necropsies
hist(data$NNecropsies)

#histovrame of log life expectancy
hist(data$logMaxLongevity)

#compute the number of species with no individuals that died from malignant cancer
length(which(data$NNeoplasia == 0)) / nrow(data) 

#histogram of the number of individuals that died from malignant cancer
hist(data$NNeoplasia)

##############################################################################################
#################### prepare the data for brms models ########################################
##############################################################################################

## fix the name for the donkeys Equus africanus asinus to match names in the csv
sel_donkey <- which(phylo_all$tip.label == "Equus_africanus")
phylo_all$tip.label[sel_donkey] <- "Equus_asinus"

## fix the name for the  Pygmy marmoset to match names in the csv
sel_donkey <- which(phylo_all$tip.label == "Cebuella_pygmaea")
phylo_all$tip.label[sel_donkey] <- "Callithrix_pygmaea"

#check the number of species left
phylo_brms <- phylo_all

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

#check if the iterations in the chains are not autocorrelated
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
saveRDS(fit_logmass_neoplasia , "Boddy_mammals_neoplasm_logMass.rds")

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
saveRDS(fit_logmass_malignant, "Boddy_mammals_malignant_logMass.rds")


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
saveRDS(fit_logMaxLongevity_neoplasia , "Boddy_mammals_neoplasm_logMaxLongevity.rds")

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
saveRDS(fit_logMaxLongevity_malignant, "Boddy_mammals_malignant_logMaxLongevity.rds")

