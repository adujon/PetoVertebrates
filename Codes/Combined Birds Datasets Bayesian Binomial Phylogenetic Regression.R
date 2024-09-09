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
phylo_all <- read.tree("Combined_data_bird_phylogeny.tre") 

##########################
### import cancer data ###
##########################

#import Kapsetaki data
data_kapsetaki <-  read.csv("Kapsetaki Cancer Birds.csv")
data_kapsetaki <- data.frame(data_kapsetaki) 
colnames(data_kapsetaki)[6] <- "NNecropsies" 

#replace a species names so it matches the same species in Compton's data
sel <- which(data_kapsetaki$Species == "Erythrura_gouldiae") 
data_kapsetaki$Species[sel] <- "Chloebia_gouldiae"

#compute the number of animal withs neoplasia and malignancy
data_kapsetaki$NNeoplasia <- as.integer(data_kapsetaki$NNecropsies * data_kapsetaki$NeoplasiaPrevalence)
data_kapsetaki$NMalignant <- as.integer(data_kapsetaki$NNecropsies * data_kapsetaki$MalignancyPrevalence)

data_kapsetaki$adult_mass_kg <- data_kapsetaki$adult_mass / 1000 #convert the body mass to kg
data_kapsetaki$logMass <- log(data_kapsetaki$adult_mass_kg) #log transform the body mass
data_kapsetaki$logMaxLongevity <- log(data_kapsetaki$max_longevity) #log transform the max longevity in days

#get only the data relevance for the analysis
data_kapsetaki <- data_kapsetaki[, c("Species", "logMass", "logMaxLongevity", "NNeoplasia", "NMalignant", "NNecropsies")]
data_kapsetaki$Study <- "Kapsetaki"

#import Compton's raw data
data_compton <-  read.csv("Compton et al Vertebrates.csv")
data_compton <- data.frame(data_compton) 
colnames(data_compton)[11] <- "NNecropsies" 
colnames(data_compton)[12] <- "NNeoplasia" 
colnames(data_compton)[16] <- "NMalignant" 

#subset the dataset to only keep birds
data_compton <- subset(data_compton, data_compton$Class == "Aves")

#change the data_frame row names to species names
row.names(data_compton) <- data_compton$Species #change row names to species names

#replace the -1 in the mass and max longevity columns by NA
data_compton$adult_weight.g[data_compton$adult_weight.g == -1] <- NA
data_compton$max_longevity.months.[data_compton$max_longevity.months. == -1] <- NA

data_compton$adult_mass_kg <- data_compton$adult_weight.g / 1000 #convert the body mass to kg
data_compton$logMass <- log(data_compton$adult_mass_kg) #log transform the body mass
data_compton$logMaxLongevity <- log(data_compton$max_longevity.months.) #log transform the max longevity in months

#keep only the data needed for the analysis
data_compton <- data_compton[, c("Species", "logMass", "logMaxLongevity", "NNeoplasia", "NMalignant", "NNecropsies")]
data_compton$Study <- "Compton"

#combine the two datasets
data <- rbind(data_kapsetaki, data_compton)
data$SpeciesRepetitions <- data$Species

##############################################################################################
############################# visualize the data #############################################
##############################################################################################

#histogram of the body mass
hist(data$logMass) #better when log transformed

#histogram of the number of necropsies
hist(data$NNecropsies)

#distribution looks ok outliers are not too big
hist(data$logMaxLongevity)

#compute the number of species with no individuals that died from malignant cancer
length(which(data$NNeoplasia == 0)) / nrow(data) 

#histogram of the number of individuals that died from malignant cancer
hist(data$NNeoplasia)

##############################################################################################
#################### prepare the data for brms models ########################################
##############################################################################################

#check the number of species left
phylo_brms <- phylo_all

#compute the var-covar matrix from the tree to use in brms
vcv_brms <- vcv(phylo_brms, model = "Brownian")

#add a tiny ammount to the diagonal of the var-covar matrix to help the model
# converge
diag(vcv_brms) <- diag(vcv_brms)+0.000000002 

##############################################################################################
#######################  Modelling the effect of logMass on neoplasia ##########################
##############################################################################################

#fit the phylogenetic regression model using brms
fit_logmass_neoplasia <- brm(NNeoplasia | trials(NNecropsies) ~ logMass + Study + (1|gr(Species, cov = vcv_brms)) + (1|SpeciesRepetitions),
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
saveRDS(fit_logmass_neoplasia , "Combined_data_birds_neoplasm_logMass.rds")

##############################################################################################
#######################  Modelling the effect of logMass on malignancy #######################
##############################################################################################

#fit the phylogenetic regression model using brms
fit_logmass_malignant <- brm(NMalignant | trials(NNecropsies) ~ logMass + Study + (1|gr(Species, cov = vcv_brms)) + (1|SpeciesRepetitions),
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

#check if the iterations in the chains are not autocorrelated
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
saveRDS(fit_logmass_malignant, "Combined_data_birds_malignant_logMass.rds")


##############################################################################################
#######################  Modelling the effect of logMaxLongevity on neoplasia ################
##############################################################################################

#fit the phylogenetic regression model using brms
fit_logMaxLongevity_neoplasia <- brm(NNeoplasia | trials(NNecropsies) ~ logMaxLongevity + Study + (1|gr(Species, cov = vcv_brms)) + (1|SpeciesRepetitions),
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
saveRDS(fit_logMaxLongevity_neoplasia , "Combined_data_birds_neoplasm_logMaxLongevity.rds")

##############################################################################################
#######################  Modelling the effect of logMaxLongevity on malignancy #######################
##############################################################################################

#fit the phylogenetic regression model using brms
fit_logMaxLongevity_malignant <- brm(NMalignant | trials(NNecropsies) ~ logMaxLongevity + Study + (1|gr(Species, cov = vcv_brms)) + (1|SpeciesRepetitions),
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
saveRDS(fit_logMaxLongevity_malignant, "Combined_data_birds_malignant_logMaxLongevity.rds")

